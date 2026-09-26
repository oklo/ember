#!/usr/bin/env python3
"""Test extra H/He3 composition nodes before generating full source planes.

All 102 independent chemical targets and native replies are reused. Direct
potential values at eight fixed material states isolate composition sampling
from temperature/density interpolation. Every new source group is retained.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import json
import math
from pathlib import Path
import subprocess
import numpy as np
from composition_potential_v2 import RGAS,PotentialSpline,source_mixing_value
from metal_eos_composition import mixture
from audit_smooth_eos_runtime import sha
from audit_diffusion_eos_two_compositions import normalized_error


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('work',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--jobs',type=int,default=8);a=p.parse_args()
    if not 1<=a.jobs<=8:raise ValueError('invalid workers')
    a.work.mkdir()
    reference_path=Path('docs/results/smooth_eos_chemical_domain_v1.json');ref=json.loads(reference_path.read_text())
    native_path=Path('/tmp/ember-smooth-eos-chemical-domain-run-v1/calculations/native.output')
    if sha(native_path)!=ref['artifacts_sha256'][str(native_path)]:raise ValueError('native comparisons changed')
    native=[json.loads(s) for s in native_path.read_text().splitlines()]
    family=Path('/tmp/ember-smooth-eos-physical-family-v1/freeeos300_gs98_z020.dat')
    base_x=list(map(float,family.read_text().splitlines()[1].split()[2:]));base_y=[0.,.06,.12]
    extra_x=[.00025,.00075];extra_y=[.03];materials=ref['materials']
    axes={i:[x for x in sorted(base_x+extra_x) if not m['dense'] or x<=.2] for i,m in enumerate(materials)}
    ys=sorted(base_y+extra_y);jobs={}
    for mi,xs in axes.items():
        for x in xs:
            for y in ys:jobs.setdefault((x,y),[]).append(mi)
    build=ref['source_builds'][0];receipt=build['receipt'];probe=Path(build['directory'])/'probe'
    if sha(probe)!=receipt['probe_sha256'] or sha(receipt['library'])!=receipt['library_sha256']:
        raise ValueError('source build changed')
    inputs={str(p):sha(p) for p in [Path(__file__),reference_path,native_path,family,
            Path('scripts/composition_potential_v2.py'),Path('scripts/metal_eos_composition.py'),probe]}
    spec=dict(scope=__doc__,materials=materials,hydrogen=base_x,helium3=base_y,added_hydrogen=extra_x,
              added_helium3=extra_y,source_groups=len(jobs),source_states=sum(len(m) for m in jobs.values()),
              source_build=build,criteria=ref['criteria'],inputs_sha256=inputs)
    (a.work/'specification.json').write_text(json.dumps(spec,indent=2)+'\n')
    def source(item):
        (x,y),mis=item;mix=mixture(x,y);scale=mix['source_mass_scale']
        request=' '.join(format(v,'.17g') for v in mix['eps'])+'\n3 223 -2\n'
        request+=''.join(f'{math.log(materials[mi]["rho"]*scale):.17g} {math.log(materials[mi]["T"]):.17g}\n' for mi in mis)
        r=subprocess.run([str(probe)],input=request,text=True,capture_output=True,timeout=20)
        values=[list(map(float,line.split())) for line in r.stdout.splitlines()]
        record=dict(X=x,Y3=y,material_indices=mis,input=request,returncode=r.returncode,
                    stdout=r.stdout,stderr=r.stderr,probe_sha256=receipt['probe_sha256'],states=[])
        if r.returncode or len(values)!=len(mis):return record
        for mi,v in zip(mis,values,strict=True):
            m=materials[mi];valid=len(v)==22 and all(map(math.isfinite,v)) and v[0]==0
            defect=None
            if valid:
                valid=abs(v[2]/(m['rho']*scale)-1)<1e-9 and abs(v[3]/m['T']-1)<1e-10 and v[4]>0 and v[10]!=0
            if valid:
                defect=max(abs(v[2]*v[9]/v[4]+v[8]-1),abs(m['T']*v[12]/v[10]-1),abs(v[2]*m['T']*v[11]/v[4]+v[8]))
                valid=defect<1e-7
            record['states'].append(dict(material=mi,valid=valid,first_law_defect=defect,
                                         phi=(v[5]/m['T']-v[6])*scale if valid else None))
        return record
    values={};failures=[]
    with (a.work/'source_groups.jsonl').open('x') as saved,ThreadPoolExecutor(a.jobs) as pool:
        futures={pool.submit(source,item):item[0] for item in jobs.items()}
        for future in as_completed(futures):
            try:r=future.result()
            except (subprocess.TimeoutExpired,ValueError) as error:
                r=dict(X=futures[future][0],Y3=futures[future][1],error=repr(error),states=[])
            saved.write(json.dumps(r,allow_nan=False)+'\n');saved.flush()
            if len(r['states'])!=len(jobs[r['X'],r['Y3']]) or not all(v['valid'] for v in r['states']):failures.append({k:v for k,v in r.items() if k not in ['input','stdout','stderr']})
            for v in r['states']:
                if v['valid']:values[v['material'],r['X'],r['Y3']]=v['phi']
    modes={'original':(False,False),'hydrogen_only':(True,False),'helium3_only':(False,True),'both':(True,True)}
    splines={}
    if not failures:
        for mi,all_x in axes.items():
            for name,(hx,hy) in modes.items():
                xx=all_x if hx else [x for x in all_x if x in base_x];yy=ys if hy else base_y
                grid=np.array([[values[mi,x,y]-source_mixing_value(x,y) for y in yy] for x in xx])
                splines[mi,name]=PotentialSpline(xx,yy,grid)
    comparisons=[]
    if not failures:
        for c,n in zip(ref['comparisons'],native,strict=True):
            if not c.get('source_reference_passes'):raise ValueError('unresolved independent chemical reference')
            gref=np.array(c['source_gradient_over_Rgas'])*RGAS;href=np.array(c['source_Hessian_over_Rgas'])*RGAS
            result={k:c[k] for k in ['material','material_name','X','Y3']};result['methods']={}
            for name in modes:
                _,g,h=splines[c['material'],name].jet(c['X'],c['Y3'])
                ge=float(np.max(np.abs(g-gref))/RGAS);he=normalized_error(h-href,href)
                result['methods'][name]=dict(gradient_error_over_Rgas=ge,normalized_Hessian_error=he,
                    passed=ge<ref['criteria']['gradient_over_Rgas'] and he<ref['criteria']['normalized_Hessian'])
                if name=='original':
                    result['original_direct_nodes_to_native_gradient_over_Rgas']=float(np.max(np.abs(g-np.array(n['values'][22:24])))/RGAS)
                    result['original_direct_nodes_to_native_Hessian_error']=normalized_error(h-np.array(n['values'][24:28]).reshape(2,2),href)
            comparisons.append(result)
    summaries={name:dict(passed=sum(c['methods'][name]['passed'] for c in comparisons),
                 maximum_gradient_error_over_Rgas=max((c['methods'][name]['gradient_error_over_Rgas'] for c in comparisons),default=None),
                 maximum_normalized_Hessian_error=max((c['methods'][name]['normalized_Hessian_error'] for c in comparisons),default=None)) for name in modes}
    passed=not failures and len(comparisons)==len(ref['comparisons']) and summaries['both']['passed']==len(comparisons)
    report=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),outcome='refinement_pilot_passed' if passed else 'failed',
                accepted_for_evolution=False,source_failures=failures,source_groups=len(jobs),source_states=spec['source_states'],
                criteria=ref['criteria'],method_summaries=summaries,comparisons=comparisons,
                inputs_sha256=inputs,artifacts_sha256={str(p):sha(p) for p in a.work.iterdir() if p.is_file()},
                remaining_work=['Generate only demonstrated necessary full composition additions under a bounded allocation',
                                'Implement and verify four-He3-plane interpolation','Repeat failed physical comparisons with existing direct references',
                                'Broader EOS acceptance and diffusion coupling'])
    with a.output.open('x') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'outcome':report['outcome'],'source_failures':len(failures),'methods':summaries}))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
