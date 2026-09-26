#!/usr/bin/env python3
"""Independent chemical-derivative checks in ionized, neutral and dense matter.

Two source quadratures and two composition spacings test the direct reference.
The known ideal source mixing is subtracted before finite differences, and
the full isotope mixing derivatives are restored analytically. New source
states are grouped by composition and saved immediately for later reuse.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path
import subprocess
import numpy as np
from composition_potential_v2 import RGAS,source_mixing_value,mixing_derivatives
from metal_eos_composition import mixture
from audit_smooth_eos_runtime import query,sha
from audit_diffusion_eos_two_compositions import normalized_error


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['probe','family','work','output']:p.add_argument(name,type=Path)
    p.add_argument('--jobs',type=int,default=8);a=p.parse_args()
    if not 1<=a.jobs<=8:raise ValueError('invalid worker count')
    a.work.mkdir()
    builds=[]
    for directory in ['/tmp/ember-freeeos-precision-builder-check-v1',
                      '/tmp/ember-freeeos-quadrature13-builder-check-v1']:
        path=Path(directory);r=json.loads((path/'build_receipt.json').read_text())
        if sha(path/'probe')!=r['probe_sha256'] or sha(r['library'])!=r['library_sha256']:
            raise ValueError('source build changed')
        builds.append(dict(directory=directory,receipt=r,receipt_sha256=sha(path/'build_receipt.json')))
    materials=[dict(name='cool_atmosphere',T=4157.,rho=7.2e-5,dense=False),
               dict(name='warm_dilute_atmosphere',T=7500.,rho=2.3e-6,dense=False),
               dict(name='hydrogen_ionization',T=15000.,rho=2.1e-5,dense=False),
               dict(name='helium_ionization_dense_envelope',T=27183.,rho=.037,dense=False),
               dict(name='helium_ionization_dilute_envelope',T=27183.,rho=3.7e-5,dense=False),
               dict(name='low_temperature_dense_extension',T=1123000.,rho=10365.007827247957,dense=True),
               dict(name='intermediate_dense_extension',T=7913000.,rho=193870.50813468645,dense=True),
               dict(name='hot_dense_extension',T=20710000.,rho=820861.9941714394,dense=True)]
    axis=np.array(list(map(float,a.family.read_text().splitlines()[1].split()[2:])))
    controls=[];tasks={};offsets={(0.,0.)}
    for s in [1.,.5]:offsets.update([(s,0),(-s,0),(0,s),(0,-s),(s,s),(s,-s),(-s,s),(-s,-s)])
    for mi,m in enumerate(materials):
        for x in ([.0005,.08125,.1625] if m['dense'] else [.0005,.08125,.1625,.55,.725]):
            cell=np.searchsorted(axis,x,side='right')-1;hx=float((axis[cell+1]-axis[cell])/16)
            for y in [.005,.03,.09]:
                hy=min(.002,y/4);control=dict(material=mi,X=x,Y3=y,hx=hx,hy=hy);controls.append(control)
                for bi in range(2):
                    for dx,dy in offsets:
                        key=(bi,x+dx*hx,y+dy*hy)
                        tasks.setdefault(key,set()).add(mi)
    criteria=dict(source_quadrature=1e-4,source_spacing=1e-4,
                  normalized_Hessian=.005,gradient_over_Rgas=1e-4)
    inputs={str(p.resolve()):sha(p) for p in [Path(__file__),a.probe,a.family,
            Path('scripts/composition_potential_v2.py'),Path('scripts/metal_eos_composition.py')]}
    specification=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),materials=materials,
                       controls=controls,source_builds=builds,source_options=[3,223,-2],
                       source_groups=len(tasks),source_states=sum(len(v) for v in tasks.values()),
                       criteria=criteria,inputs_sha256=inputs)
    (a.work/'specification.json').write_text(json.dumps(specification,indent=2)+'\n')
    def direct(job):
        (bi,x,y),mis=job;mis=sorted(mis);mix=mixture(x,y);scale=mix['source_mass_scale']
        request=' '.join(format(v,'.17g') for v in mix['eps'])+'\n3 223 -2\n'
        request+=''.join(f'{math.log(materials[mi]["rho"]*scale):.17g} {math.log(materials[mi]["T"]):.17g}\n' for mi in mis)
        command=[str(Path(builds[bi]['directory'])/'probe')]
        record=dict(build_index=bi,probe_sha256=builds[bi]['receipt']['probe_sha256'],X=x,Y3=y,
                    material_indices=mis,input=request,input_sha256=hashlib.sha256(request.encode()).hexdigest())
        try:
            r=subprocess.run(command,input=request,text=True,capture_output=True,timeout=20)
            rows=[list(map(float,line.split())) for line in r.stdout.splitlines()]
            record.update(returncode=r.returncode,stdout=r.stdout,stderr=r.stderr,
                          output_complete=r.returncode==0 and len(rows)==len(mis) and
                          all(len(v)==22 and all(map(math.isfinite,v)) for v in rows))
            records=[]
            if record['output_complete']:
                for mi,v in zip(mis,rows,strict=True):
                    m=materials[mi];P=v[4];T=m['T'];rho=m['rho'];reason=None
                    if v[0]!=0:reason='source_nonconvergence'
                    elif abs(v[2]/(rho*scale)-1)>1e-9 or abs(v[3]/T-1)>1e-10:reason='source_coordinate_mismatch'
                    elif P<=0 or v[10]==0:reason='source_nonpositive_pressure_or_zero_heat_response'
                    else:
                        defect=max(abs(v[2]*v[9]/P+v[8]-1),abs(T*v[12]/v[10]-1),abs(v[2]*T*v[11]/P+v[8]))
                        if defect>1e-7:reason='source_first_law_failure'
                    records.append(dict(material=mi,source_info=v[0],valid=reason is None,reason=reason,
                                        phi=(v[5]/T-v[6])*scale if reason is None else None))
            record['states']=records
        except (subprocess.TimeoutExpired,ValueError) as error:
            record.update(output_complete=False,error=repr(error))
            if isinstance(error,subprocess.TimeoutExpired):
                record['stdout']=(error.stdout or b'').decode(errors='replace') if isinstance(error.stdout,bytes) else error.stdout
                record['stderr']=(error.stderr or b'').decode(errors='replace') if isinstance(error.stderr,bytes) else error.stderr
        return record
    sources={};source_failures=[]
    with (a.work/'source_groups.jsonl').open('x') as saved,ThreadPoolExecutor(a.jobs) as pool:
        futures=[pool.submit(direct,t) for t in sorted(tasks.items())]
        for i,future in enumerate(as_completed(futures),1):
            r=future.result();saved.write(json.dumps(r,allow_nan=False)+'\n');saved.flush()
            if not r['output_complete']:source_failures.append({k:v for k,v in r.items() if k not in ['input','stdout','stderr']})
            for v in r.get('states',[]):
                key=(r['build_index'],v['material'],r['X'],r['Y3'])
                sources[key]=v
                if not v['valid']:source_failures.append(dict(key=key,**v))
            if i%50==0:print(json.dumps({'source_groups_saved':i,'source_failures':len(source_failures)}),flush=True)
    points=[[c['X'],c['Y3'],materials[c['material']]['T'],materials[c['material']]['rho'],1] for c in controls]
    native=query(a.probe,a.family,points,a.work,'native')
    def source_jet(c,bi):
        x,y,hx,hy=c['X'],c['Y3'],c['hx'],c['hy']
        def f(dx,dy):
            xx=x+dx*hx;yy=y+dy*hy;v=sources[bi,c['material'],xx,yy]
            if not v['valid']:raise ValueError(v['reason'])
            return np.longdouble(v['phi'])-np.longdouble(source_mixing_value(xx,yy))
        def difference(s):
            center=f(0,0);g=np.array([(f(s,0)-f(-s,0))/(2*s*hx),(f(0,s)-f(0,-s))/(2*s*hy)],dtype=float)
            xx=((f(s,0)-center)+(f(-s,0)-center))/(s*hx)**2
            yy=((f(0,s)-center)+(f(0,-s)-center))/(s*hy)**2
            xy=(f(s,s)-f(s,-s)-f(-s,s)+f(-s,-s))/(4*s*s*hx*hy)
            return g,np.array([[xx,xy],[xy,yy]],dtype=float)
        g0,h0=difference(1);g1,h1=difference(.5);gi,hi=mixing_derivatives(x,y)
        return gi+(4*g1-g0)/3,hi+(4*h1-h0)/3,hi+h1
    results=[]
    for c,native_row in zip(controls,native,strict=True):
        result=dict(c,material_name=materials[c['material']]['name'])
        try:
            if not native_row['ok']:raise ValueError(native_row['error'])
            g0,h0,_=source_jet(c,0);g,h,hfine=source_jet(c,1)
            precision=normalized_error(h0-h,h);spacing=normalized_error(hfine-h,h)
            v=np.array(native_row['values']);gradient=float(np.max(np.abs(v[22:24]-g))/RGAS)
            error=normalized_error(v[24:28].reshape(2,2)-h,h)
            result.update(source_gradient_over_Rgas=(g/RGAS).tolist(),source_Hessian_over_Rgas=(h/RGAS).tolist(),
                          source_quadrature_error=precision,source_spacing_error=spacing,
                          native_Hessian_error=error,native_gradient_error_over_Rgas=gradient,
                          source_reference_passes=precision<criteria['source_quadrature'] and spacing<criteria['source_spacing'],
                          native_comparison_passes=error<criteria['normalized_Hessian'] and gradient<criteria['gradient_over_Rgas'])
            result['passed']=result['source_reference_passes'] and result['native_comparison_passes']
        except (ValueError,KeyError) as error:result.update(passed=False,error=repr(error))
        results.append(result)
    passed=all(c['passed'] for c in results) and not source_failures
    report=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed' if passed else 'failed',
                accepted_for_evolution=False,source_groups=specification['source_groups'],source_states=specification['source_states'],
                source_failures=source_failures,materials=materials,comparisons=results,criteria=criteria,
                inputs_sha256=inputs,source_builds=builds,
                artifacts_sha256={str(p):sha(p) for p in a.work.iterdir() if p.is_file()},
                remaining_work=['Resolve any failed numerical reference or composition interpolation controls',
                                'Composition-boundary and source-mask physics','Consistent diffusion forces and evolution coupling'])
    with a.output.open('x') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'outcome':report['outcome'],'comparisons':len(results),'passed':sum(c['passed'] for c in results),
                      'source_failures':len(source_failures)}))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
