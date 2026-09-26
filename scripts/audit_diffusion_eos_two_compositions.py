#!/usr/bin/env python3
"""Compare H/He3 free-energy derivatives with independent FreeEOS queries.

Use existing Y3=0/.12 source planes and test adding a Y3=.06 plane at three
fixed material coordinates. Independent targets at Y3=.03/.09 supply both
diagonal and mixed composition derivatives. No full source plane is generated
and no production input is changed. All new direct evaluations are saved.
"""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess

import numpy as np

from composition_potential import (RGAS,PotentialSpline,mixing_value,
                                   source_mixing_value,mixing_derivatives)
from metal_eos_composition import mixture


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def normalized_error(delta,reference):
    e,q=np.linalg.eigh(reference)
    if np.min(e)<=0:raise ValueError('reference free-energy curvature is not positive')
    inv=(q/np.sqrt(e))@q.T
    return float(np.max(np.abs(np.linalg.eigvalsh(inv@delta@inv))))


def main():
    root=Path(__file__).resolve().parents[1]
    work=Path('/tmp/ember-diffusion-eos-two-compositions-run-v1')
    parent=root/'docs/results/diffusion_eos_composition_curvature_v2.json'
    prior=json.loads(parent.read_text())
    for p,h in prior['EOS_inputs_sha256'].items():assert digest(p)==h
    family=root/'data/eos/numerical_electron_base_v1/freeeos300_gs98_z020.dat'
    axis=np.array([float(v) for v in family.read_text().splitlines()[1].split()[2:]])
    with (root/'docs/reports/2026-09-11/evolution_latest.csv').open() as f:
        history=list(csv.DictReader(f))
    peak=max(history,key=lambda r:float(r['central_Y3']))
    c0=prior['controls'][0];c2=prior['controls'][2]
    material=[{'name':'late_center','T':c0['T'],'rho':c0['rho'],'X':c0['X']},
              {'name':'late_outer_nonconvective','T':c2['T'],'rho':c2['rho'],'X':c2['X']},
              {'name':'helium3_peak_center','T':float(peak['central_T_K']),
               'rho':float(peak['central_rho']),'X':.55}]
    old_stdout=Path('/tmp/ember-diffusion-eos-spline-run-v1/stdout.txt')
    old_report=root/'docs/results/diffusion_eos_hydrogen_spline_v1.json'
    old=json.loads(old_report.read_text());assert digest(old_stdout)==old['output_sha256'][str(old_stdout)]
    cached={tuple(row[:4]):row for row in np.loadtxt(old_stdout)}
    queries=[]
    for m in material:
        for x in axis:
            for y in [0.,.12]:
                key=(x,y,m['T'],m['rho'])
                if key not in cached:queries.append(key)
    inp=''.join(' '.join(format(v,'.17g') for v in k)+'\n' for k in queries)
    runtime_probe=Path('/tmp/ember-diffusion-eos-curvature-run-v1/probe')
    p=subprocess.run([str(runtime_probe),str(family)],input=inp,text=True,capture_output=True,timeout=90)
    (work/'runtime_input.txt').write_text(inp);(work/'runtime_stdout.txt').write_text(p.stdout)
    (work/'runtime_stderr.txt').write_text(p.stderr)
    if p.returncode:raise RuntimeError(p.stderr)
    runtime=np.array([list(map(float,l.split())) for l in p.stdout.splitlines()])
    assert runtime.shape==(len(queries),10) and np.all(np.isfinite(runtime))
    for key,row in zip(queries,runtime,strict=True):
        assert np.allclose(key,row[:4],rtol=1e-15,atol=0)
        cached[key]=row
    builds=prior['source_builds']
    for b in builds:
        assert digest(Path(b['directory'])/'probe')==b['receipt']['probe_sha256']
        assert digest(b['receipt']['library'])==b['receipt']['library_sha256']
    targets=[]
    for mi,m in enumerate(material):
        x=m['X'];j=np.searchsorted(axis,x,side='right')-1
        hx=float((axis[j+1]-axis[j])/16)
        for y in [.03,.09]:targets.append({'material':mi,'X':x,'Y3':y,'hx':hx,'hy':.002})
    offsets={(0.,0.)}
    for h in [1.,.5]:
        offsets.update([(h,0),(-h,0),(0,h),(0,-h),
                        (h,h),(h,-h),(-h,h),(-h,-h)])
    tasks=set()
    for mi,m in enumerate(material):
        for x in axis:tasks.add((0,mi,float(x),.06))
    for t in targets:
        for bi in range(2):
            for dx,dy in offsets:tasks.add((bi,t['material'],t['X']+dx*t['hx'],t['Y3']+dy*t['hy']))
    tasks=sorted(tasks)
    def direct(key):
        bi,mi,x,y=key;m=material[mi];comp=mixture(x,y);scale=comp['source_mass_scale']
        text=' '.join(format(v,'.17g') for v in comp['eps'])+'\n3 223 -2\n'
        text+=f'{math.log(m["rho"]*scale):.17g} {math.log(m["T"]):.17g}\n'
        p=subprocess.run([str(Path(builds[bi]['directory'])/'probe')],input=text,
                         text=True,capture_output=True,timeout=20)
        v=list(map(float,p.stdout.split()))
        row={'key':key,'input':text,'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr}
        if p.returncode or len(v)!=22 or v[0]!=0 or not all(map(math.isfinite,v)):
            raise RuntimeError(row)
        assert abs(v[2]/(m['rho']*scale)-1)<1e-9 and abs(v[3]/m['T']-1)<1e-10
        row.update(phi=(v[5]/m['T']-v[6])*scale,P=v[4],E=v[5]*scale)
        return row
    direct_records=[]
    with (work/'direct_source.jsonl').open('x') as f:
        with ThreadPoolExecutor(max_workers=4) as pool:
            for v in pool.map(direct,tasks):
                direct_records.append(v);f.write(json.dumps(v,allow_nan=False)+'\n');f.flush()
    sources={tuple(v['key']):v for v in direct_records}
    splines={}
    for mi,m in enumerate(material):
        ys=np.array([0.,.06,.12]);r=np.empty((len(axis),3))
        for i,x in enumerate(axis):
            for j,y in enumerate(ys):
                if y==.06:
                    prad=4*5.670374419e-5/2.99792458e10*m['T']**4/3
                    r[i,j]=(sources[0,mi,x,y]['phi']-source_mixing_value(x,y)
                            -prad/(m['rho']*m['T']))
                else:r[i,j]=cached[x,y,m['T'],m['rho']][4]-mixing_value(x,y)
        for label,indices in [('two_planes',[0,2]),('three_planes',[0,1,2])]:
            splines[mi,label]=PotentialSpline(axis,ys[indices],r[:,indices])
    def source_jet(t,bi):
        x,y,hx,hy=t['X'],t['Y3'],t['hx'],t['hy']
        def f(dx,dy):
            xx=x+dx*hx;yy=y+dy*hy
            return np.longdouble(sources[bi,t['material'],xx,yy]['phi'])-np.longdouble(source_mixing_value(xx,yy))
        def difference(s):
            middle=f(0,0)
            g=np.array([(f(s,0)-f(-s,0))/(2*s*hx),(f(0,s)-f(0,-s))/(2*s*hy)],dtype=float)
            xx=((f(s,0)-middle)+(f(-s,0)-middle))/(s*hx)**2
            yy=((f(0,s)-middle)+(f(0,-s)-middle))/(s*hy)**2
            xy=(f(s,s)-f(s,-s)-f(-s,s)+f(-s,-s))/(4*s*s*hx*hy)
            return g,np.array([[xx,xy],[xy,yy]],dtype=float)
        g0,h0=difference(1);g1,h1=difference(.5)
        gi,hi=mixing_derivatives(x,y)
        return gi+(4*g1-g0)/3,hi+(4*h1-h0)/3,hi+h1
    comparisons=[]
    for t in targets:
        g0,h0,_=source_jet(t,0);g1,h1,fine=source_jet(t,1)
        precision=normalized_error(h0-h1,h1);spacing=normalized_error(fine-h1,h1)
        assert precision<1e-4 and spacing<1e-4,(t,precision,spacing)
        result={**t,'material_name':material[t['material']]['name'],
                'source_gradient_over_Rgas':(g1/RGAS).tolist(),
                'source_Hessian_over_Rgas':(h1/RGAS).tolist(),
                'source_precision_error':precision,'source_spacing_error':spacing,'methods':{}}
        for name in ['two_planes','three_planes']:
            v,g,h=splines[t['material'],name].jet(t['X'],t['Y3'])
            err=normalized_error(h-h1,h1)
            result['methods'][name]={'gradient_over_Rgas':(g/RGAS).tolist(),
                                    'Hessian_over_Rgas':(h/RGAS).tolist(),
                                    'normalized_Hessian_error':err,
                                    'maximum_gradient_error_over_Rgas':float(np.max(np.abs(g-g1))/RGAS),
                                    'passes_0p5_percent_Hessian_criterion':err<.005}
        comparisons.append(result)
    report={'created_utc':datetime.now(timezone.utc).isoformat(),'scope':__doc__,
            'outcome':'completed_composition_derivative_pilot','accepted_for_evolution':False,
            'material_coordinates':material,'comparisons':comparisons,
            'new_runtime_queries':len(queries),'reused_runtime_queries':2*len(axis)*len(material)-len(queries),
            'new_direct_source_queries':len(direct_records),'source_builds':builds,
            'criterion':'Spectral norm of H_source^(-1/2) (H_candidate-H_source) H_source^(-1/2) below 0.005; source spacing and quadrature comparisons below 1e-4.',
            'remaining_work':'Temperature/density and composition cross responses, source masks and boundaries, stored composition slopes and runtime cost; chemical forces, thermal diffusion and conservative stellar coupling.',
            'inputs_sha256':{str(p):digest(p) for p in [Path(__file__),root/'scripts/composition_potential.py',
                root/'scripts/metal_eos_composition.py',parent,old_report,old_stdout,runtime_probe,family,
                root/'docs/reports/2026-09-11/evolution_latest.csv']},
            'output_sha256':{str(p):digest(p) for p in [work/'runtime_input.txt',work/'runtime_stdout.txt',work/'runtime_stderr.txt',work/'direct_source.jsonl']}}
    out=root/'docs/results/diffusion_eos_two_compositions_v1.json'
    with out.open('x') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({name:max(c['methods'][name]['normalized_Hessian_error'] for c in comparisons)
                      for name in ['two_planes','three_planes']}))


if __name__=='__main__':main()
