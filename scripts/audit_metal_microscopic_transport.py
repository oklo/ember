#!/usr/bin/env python3
"""Finite-difference and limiting checks for coupled H/He/metal drift and heat."""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import time
import numpy as np
from write_scientific_result import write_result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('probe',type=Path);p.add_argument('eos',type=Path);p.add_argument('collision',type=Path)
    p.add_argument('checkpoint',type=Path);p.add_argument('output',type=Path);p.add_argument('--scratch',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists() or a.scratch.exists():raise FileExistsError('preserve controls')
    a.scratch.mkdir(parents=True)
    root=Path(__file__).resolve().parents[1]
    paths=[a.probe,a.eos,a.collision,a.checkpoint,Path(__file__)]+[root/v for v in [
        'src/screened_metal_microscopic_transport.cpp','src/metal_microscopic_transport.cpp',
        'include/ember/metal_microscopic_transport.hpp','scripts/conditional_metal_envelope_heat.hpp']]
    sha=lambda v:hashlib.sha256(Path(v).read_bytes()).hexdigest()
    hashes={str(v.resolve()):sha(v) for v in paths}
    for v in paths:
        if v.suffix in ['.cpp','.hpp','.py']:shutil.copy2(v,a.scratch/v.name)
    shutil.copy2(a.probe,a.scratch/'probe')
    queries=[];stencils=[];checks=[]
    def add(q):queries.append(list(q));return len(queries)-1
    def point(row):
        mass,r,rho,T,L,x,y,c12,c13,n14=row
        z=.02+12*c12+13*c13+14*n14-.02*(.171836+.050335)
        return [math.log(r),math.log(rho),math.log(T),L,x,y,z]
    rows=json.loads(a.checkpoint.read_text())['model_record']['model']
    def face(i,ions=1,rad=1,mode='face'):
        return [mode,ions,1,7,rad,rows[i][0],rows[i+1][0]]+point(rows[i])+point(rows[i+1])
    # Avoid unresolved central differences for numerical differentiation;
    # include the actual burning layer and the convective envelope boundary.
    bases=[face(i,ions,rad) for i in [200,276,300,350,386] for ions,rad in [(0,0),(1,1)]]
    # A finite metal gradient independently exercises the new abundance direction.
    q=face(386);q[20]=.5*q[13];bases.append(q)
    # Inactive species and full thermal subspaces retain exactly zero flux.
    q=face(386);q[12]=q[19]=0;q[3]=5;bases.append(q)
    q=face(386);q[13]=q[20]=0;q[3]=3;bases.append(q)
    for base in bases:
        i=add(base);mask=base[3]
        for end in range(2):
            for k in range(3):
                if not mask&(1<<k):continue
                coord=11+7*end+k
                helium=1-math.fsum(base[11+7*end:14+7*end])
                scale=min(base[coord],helium)
                h=1e-4*scale;ids=[]
                for step in [-2,-1,1,2]:
                    q=base.copy();q[2]=0;q[coord]+=step*h;ids.append(add(q))
                stencils.append(dict(kind='species',base=i,end=end,coordinate=k,h=h,indices=ids))
        for kind in ['heat','total']:
            q=base.copy();q[0]=kind
            if kind=='total':q.extend([1.1e12,-2e9,-2.5e11] if mask==7 else [1.1e12,-2e9 if mask&2 else 0.,-2.5e11 if mask&4 else 0.])
            j=add(q)
            if kind=='heat':checks.append(dict(kind='same_heat',face=i,heat=j))
            for end in range(2):
                for k in range(3):
                    ids=[];h=2e-5
                    for step in [-2,-1,1,2]:
                        trial=q.copy();trial[2]=0;trial[7+7*end+k]+=step*h;ids.append(add(trial))
                    stencils.append(dict(kind='thermal',base=j,end=end,coordinate=k,h=h,indices=ids))
            if kind=='total':
                for k in range(3):
                    if not mask&(1<<k):continue
                    trial=q.copy();trial[21+k]+=1e11;checks.append(dict(kind='linear',base=j,other=add(trial),species=k,change=1e11))
        checks.append(dict(kind='wrong',index=add(['wrong']+base[1:]),accepted=False))
    # Check the cold convective prescription and its smooth temperature join.
    for temp in [1e5,2.5e6,3.5e6]:
        q=face(386,mode='envelope');q.extend([1.1e12,-2e9,-2.5e11])
        q[9]=math.log(temp);q[16]=math.log(temp*.99)
        # Supported low-density thermal source track with fixed log Q=.3.
        q[8]=math.log(10**(.3)*(temp/1e6)**1.5);q[15]=q[8]-.02
        j=add(q)
        for end in range(2):
            for k in range(3):
                h=2e-5;ids=[]
                for step in [-2,-1,1,2]:
                    trial=q.copy();trial[2]=0;trial[7+7*end+k]+=step*h;ids.append(add(trial))
                stencils.append(dict(kind='thermal',base=j,end=end,coordinate=k,h=h,indices=ids))
    q=face(386);q[9]=q[16]=math.log(1e6);checks.append(dict(kind='cold_domain',index=add(q),accepted=False))
    text=''.join(' '.join(str(v) if isinstance(v,str) else format(v,'.17g') for v in q)+'\n' for q in queries)
    start=time.monotonic()
    run=subprocess.run([str(a.scratch/'probe'),str(a.eos),str(a.collision)],input=text,text=True,capture_output=True,check=True,timeout=180)
    elapsed=time.monotonic()-start
    replies=[json.loads(v) for v in run.stdout.splitlines()]
    if len(replies)!=len(queries):raise ValueError('missing responses')
    failures=[];maximum=dict(species=0.,thermal=0.,linear=0.,same_heat=0.)
    details=[]
    for s in stencils:
        r=replies[s['base']];parts=[replies[j] for j in s['indices']]
        if 'error' in r or any('error' in v for v in parts):
            failures.append(dict(kind='domain',stencil=s,errors=[v.get('error') for v in [r]+parts]));continue
        if s['kind']=='species':
            vec=lambda v:np.array(v['rate'])
            exact=np.array(r['dright' if s['end'] else 'dleft'])[:,s['coordinate']]
            norm=np.abs(exact)+np.abs(vec(r))/queries[s['base']][11+7*s['end']+s['coordinate']]
        else:
            vec=lambda v:np.array([v['carried'],v['conductivity']])
            suffix='hi' if s['end'] else 'lo'
            exact=np.array([r['dcarried_'+suffix][s['coordinate']],r['dconductivity_'+suffix][s['coordinate']]])
            norm=np.abs(exact)+np.abs(vec(r))
        reference=sum(w*vec(v) for w,v in zip([1,-8,8,-1],parts))/(12*s['h'])
        error=float(max(abs(reference-exact)/np.maximum(norm,1.)))
        maximum[s['kind']]=max(maximum[s['kind']],error)
        details.append(dict(**s,error=error))
        if error>5e-5:failures.append(dict(kind=s['kind'],error=error,stencil=s))
    for c in checks:
        kind=c['kind']
        if 'accepted' in c:
            if ('error' not in replies[c['index']])!=c['accepted']:failures.append(c)
            continue
        left=replies[c.get('base',c.get('face'))];right=replies[c.get('other',c.get('heat'))]
        if 'error' in left or 'error' in right:failures.append(dict(check=c,errors=[left.get('error'),right.get('error')]));continue
        if kind=='same_heat':error=abs(left['carried']-right['carried'])/max(abs(left['carried']),1.)
        else:
            expected=left['enthalpy'][c['species']]*c['change']
            error=abs(right['carried']-left['carried']-expected)/max(abs(expected),abs(left['carried']),1.)
        maximum[kind]=max(maximum[kind],error)
        if error>1e-12:failures.append(dict(check=c,error=error))
    (a.scratch/'queries.txt.gz').write_bytes(gzip.compress(text.encode(),mtime=0))
    (a.scratch/'responses.jsonl.gz').write_bytes(gzip.compress(run.stdout.encode(),mtime=0))
    (a.scratch/'stderr.txt').write_text(run.stderr)
    for v,h in hashes.items():
        if sha(v)!=h:raise ValueError(f'changed input: {v}')
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed' if not failures else 'failed',
      selected_for_evolution=False,queries=len(queries),elapsed_seconds=elapsed,maximum_errors=maximum,
      failures=failures,stencils=details,input_sha256=hashes,
      artifacts_sha256={str(v):sha(v) for v in a.scratch.iterdir()},
      limitations=['Conditional fully stripped bulk metal grouping; these checks verify derivatives and heat accounting, not the ionization assumption.',
                  'Cold convective branch uses an analytic conductivity only for checking the integration interface.',
                  'Finite whole-star evolution and variable opacity response remain separate controls.'])
    write_result(a.output,report)
    print(json.dumps({k:v for k,v in report.items() if k in ['outcome','queries','elapsed_seconds','maximum_errors','failures']},indent=2))
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
