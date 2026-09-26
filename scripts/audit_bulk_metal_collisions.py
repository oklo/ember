#!/usr/bin/env python3
"""Check the moving GS98 group against the stationary-metal constraint.

Eliminating the new metal flux from the full response must recover the
original H/He/heat matrix. Finite differences check all six physical input
derivatives, including absent-species limits. This checks the grouping's
implementation, not the accuracy of a common settling speed for each metal.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import gzip,hashlib,json,math
from pathlib import Path
import subprocess
import numpy as np
from write_scientific_result import write_result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('stage',type=Path);ap.add_argument('output',type=Path);a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    source=Path('docs/results/thoul_diffusion_comparison_v1.json')
    table=Path('/tmp/ember-collision-eta32-lowb-192-v1.dat')
    records=json.loads(source.read_text())['records'];queries=[];base=[];stencils=[]
    for record in records:
        r=record['native']
        for z in [.02,.001,0.]:
            base.append(len(queries));queries.append([r['T'],r['rho'],r['X'],r['Y3'],z,r['screening_length']])
    # Removing hydrogen at a tenuous outer face would leave this table's
    # electron-density range. Use a supported interior state for zero-species
    # limits rather than changing any physical domain guard.
    q=queries[0]
    for x,y,z in [(0.,q[3],.02),(q[2],0.,.02),(0.,0.,.02),(0.,0.,0.)]:
        base.append(len(queries));queries.append([q[0],q[1],x,y,z,q[5]])
    for i in [0,24,48,72,78,79,80,81,82,83,84]:
        q=queries[i];helium=1-math.fsum(q[2:5])
        for k in range(6):
            if 2<=k<=4 and q[k]==0:continue
            scale=1. if k in [0,1,5] else min(q[k],helium)
            ids=[];width=1e-4
            for step in [-2,-1,1,2]:
                v=q.copy()
                if k in [0,1,5]:v[k]*=math.exp(step*width)
                else:v[k]+=step*width*scale
                ids.append(len(queries));queries.append(v)
            stencils.append(dict(index=i,coordinate=k,scale=scale,width=width,indices=ids))
    text=''.join(' '.join(format(v,'.17g') for v in q)+'\n' for q in queries)
    def run(name,args=()):
        p=subprocess.run([str(a.stage/name),str(table),*args],input=text,text=True,capture_output=True,check=True,timeout=180)
        (a.stage/(name+'.jsonl.gz')).write_bytes(gzip.compress(p.stdout.encode(),mtime=0))
        r=[json.loads(line) for line in p.stdout.splitlines()]
        if len(r)!=len(queries):raise ValueError('missing collision replies')
        return r
    with ThreadPoolExecutor(3) as pool:
        jobs=[pool.submit(run,'probe'),pool.submit(run,'old_probe',['derivatives']),pool.submit(run,'current_probe',['derivatives'])]
        replies,old,current=[job.result() for job in jobs]
    (a.stage/'queries.txt.gz').write_bytes(gzip.compress(text.encode(),mtime=0))
    failures=[];maximum=dict(constrained_response=0.,fixed_conductivity=0.,derivatives=0.,solve_backward_error=0.,old_physical_response=0.)
    def check(kind,error,tolerance,**context):
        maximum[kind]=max(maximum[kind],float(error))
        if error>tolerance:failures.append(dict(kind=kind,error=float(error),tolerance=tolerance,**context))
    def vector(r):
        return np.r_[np.array(r['mobility']).ravel(),r['enthalpy'],r['conductivity'],r['energy_scale'],r['eta'],r['ne'],r['b']]
    def scales(r):
        diag=np.diag(r['mobility']);e0=r['energy_scale']
        return np.r_[np.sqrt(abs(np.outer(diag,diag))).ravel(),np.maximum(abs(np.array(r['enthalpy'])),e0),
                     r['conductivity'],e0,max(abs(r['eta']),1),r['ne'],r['b']]
    exact_matches=0
    for i,(before,after) in enumerate(zip(old,current)):
        before.pop('seconds',None);after.pop('seconds',None)
        exact_matches+=before==after
        if 'error' not in before and 'error' not in after:
            norm=scales(before)
            error=float(max(abs(vector(before)-vector(after))/np.maximum(norm,1e-250)))
            for k in range(6):
                width=1 if k in [0,1,5] else min(queries[i][k],1-math.fsum(queries[i][2:5]))
                first,second=vector(before['partials'][k])*width,vector(after['partials'][k])*width
                error=max(error,float(max(abs(first-second)/np.maximum(norm+abs(first),1e-250))))
            check('old_physical_response',error,1e-11,index=i)
        elif before!=after:failures.append(dict(kind='old_domain_changed',index=i))
        if 'error' in replies[i]:failures.append(dict(kind='collision_domain',index=i,error=replies[i]['error']))
    checked=[]
    for i in base:
        if 'error' in replies[i]:continue
        r=replies[i];b=r['bulk'];f=r['fixed'];M=np.array(b['mobility']);F=np.array(f['mobility']);keep=[0,1,3]
        constrained=M[np.ix_(keep,keep)]
        if queries[i][4]>0:constrained=constrained-np.outer(M[keep,2],M[2,keep])/M[2,2]
        scale=np.sqrt(abs(np.outer(np.diag(F),np.diag(F))))
        error=float(np.max(abs(constrained-F)/np.maximum(scale,1e-250)))
        check('constrained_response',error,3e-11,index=i)
        check('fixed_conductivity',abs(b['conductivity']/f['conductivity']-1),3e-11,index=i)
        check('solve_backward_error',b['backward_error'],1e-12,index=i)
        active=[k for k,v in enumerate(b['active']) if v]+[3];sub=M[np.ix_(active,active)]
        norm=sub/np.sqrt(np.outer(np.diag(sub),np.diag(sub)));minimum=float(np.linalg.eigvalsh(norm)[0])
        if minimum<=0:failures.append(dict(kind='nonpositive_dissipation',index=i,minimum=minimum))
        if b['active']!=[int(queries[i][k]>0) for k in [2,3,4]]:failures.append(dict(kind='active_set',index=i))
        checked.append(dict(index=i,query=queries[i],constrained_error=error,minimum_normalized_eigenvalue=minimum))
    derivative_checks=[]
    for s in stencils:
        i,k=s['index'],s['coordinate']
        if 'error' in replies[i] or any('error' in replies[j] for j in s['indices']):continue
        r=replies[i]['bulk'];parts=[vector(replies[j]['bulk']) for j in s['indices']]
        reference=(parts[0]-8*parts[1]+8*parts[2]-parts[3])/(12*s['width'])
        exact=vector(r['partials'][k])*s['scale']
        error=float(np.max(abs(reference-exact)/np.maximum(scales(r)+abs(exact),1e-250)))
        check('derivatives',error,3e-6,index=i,coordinate=k)
        derivative_checks.append(dict(**s,error=error))
    paths=[Path(__file__),source,table,a.stage/'probe',a.stage/'old_probe',a.stage/'current_probe',
           Path('src/collision_transport.cpp'),Path('include/ember/collision_transport.hpp'),
           Path('scripts/bulk_metal_collision_probe.cpp'),a.stage/'before/collision_transport.cpp',
           a.stage/'before/collision_transport.hpp']
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    result=dict(utc=datetime.now(timezone.utc).isoformat(),outcome='bulk_metal_collision_controls_passed' if not failures else 'bulk_metal_collision_controls_failed',
                selected=False,queries=len(queries),states=len(base),old_response_exact_matches=exact_matches,
                maximum_errors=maximum,failures=failures,checks=checked,derivative_checks=derivative_checks,
                input_sha256={str(p.resolve()):sha(p) for p in paths},
                limitations=['GS98 ions share one mass velocity; their heat variables and collision charges remain separate.',
                             'Fully stripped ions; no selection for cool layers or independent element separation.',
                             'No new stellar model or atmosphere has been evolved by these controls.'])
    write_result(a.output,result)
    print(json.dumps({k:v for k,v in result.items() if k in ['outcome','queries','states','old_response_exact_matches','maximum_errors','failures']},indent=2))
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
