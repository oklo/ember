#!/usr/bin/env python3
"""Verify native EOS exchange enthalpies and their analytic derivatives.

The existing interpolation is evaluated in one persistent native process.
Every saved stellar state is compared with retained EOS replies; selected
thermal and composition derivatives use two independent difference widths.
No new FreeEOS source calculation is made.
"""
import argparse
import csv
from datetime import datetime,timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('probe',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--scratch',type=Path,required=True)
    p.add_argument('--reuse-raw',type=Path);args=p.parse_args()
    assert not args.output.exists() and not args.scratch.exists();args.scratch.mkdir()
    root=Path(__file__).resolve().parents[1];sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    family=Path('/private/tmp/ember-smooth-eos-refined-family-v1/freeeos300_gs98_z020.dat')
    assert sha(family)=='7d6bfa0b26264834307dad651852bc5ff693f91f3861b4934e582dee047fd7f5'
    cached_input=Path('/tmp/ember-smooth-eos-refined-physical-run-v1/retained/retained.input')
    cached_output=cached_input.with_suffix('.output')
    assert sha(cached_input)=='850c2a243e1a2fe12ddaabcf1d466e3ecb761c9e3fdada10fd1b3f5751e84e77'
    assert sha(cached_output)=='8bc24b90d0eba4a3933313178d600bac883765cd70b576cf3c3c7d9631e3b3a2'
    cached={}
    for q,line in zip(cached_input.read_text().splitlines(),cached_output.read_text().splitlines(),strict=True):
        fields=list(map(float,q.split()));reply=json.loads(line)
        if fields[4]==1 and reply['ok']:cached[tuple(fields[:4])]=reply['values']
    profile=root/'docs/reports/2026-09-11/evolution_latest_profile.csv'
    with profile.open() as f:rows=[{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]
    queries=[];meta=[]
    for i,row in enumerate(rows):
        q=[row['X'],row['Y3'],row['temperature_K'],row['density_g_cm3'],1,1]
        assert tuple(q[:4]) in cached,i
        queries.append(q);meta.append(dict(kind='saved',zone=i))
    for T,rho in [(4157.,7.2e-5),(27183.,.037),(2341000.,13.7),(9731000.,260.7)]:
        for x,y in [(0.,.001),(.05,0.),(0.,0.)]:
            queries.append([x,y,T,rho,int(x>0),int(y>0)]);meta.append(dict(kind='zero_species'))
    controls=[0,40,100,160,220,300,360,396,421,460,480,493,505,511]+list(range(512,len(queries)))
    stencils=[]
    for i in controls:
        q=queries[i]
        for k in range(4):
            if k>=2 and q[k+2]==0:continue
            scale=1. if k<2 else q[k-2]
            for width in [1e-4,5e-5]:
                indices=[]
                for step in [-2,-1,1,2]:
                    trial=q.copy()
                    if k<2:trial[k+2]*=math.exp(step*width)
                    else:trial[k-2]+=step*width*scale
                    indices.append(len(queries));queries.append(trial)
                stencils.append(dict(index=i,coordinate=k,scale=scale,width=width,replies=indices))
    text=''.join(' '.join(format(float(v),'.17g') for v in q)+'\n' for q in queries)
    cmd=[str(args.probe),str(family)]
    if args.reuse_raw:
        assert gzip.decompress((args.reuse_raw/'queries.txt.gz').read_bytes()).decode()==text
        response_text=gzip.decompress((args.reuse_raw/'responses.jsonl.gz').read_bytes()).decode()
        stderr=(args.reuse_raw/'stderr.txt').read_text();elapsed=None
    else:
        start=time.monotonic();run=subprocess.run(cmd,input=text,text=True,capture_output=True,timeout=150,check=True);elapsed=time.monotonic()-start
        response_text,stderr=run.stdout,run.stderr
    (args.scratch/'queries.txt.gz').write_bytes(gzip.compress(text.encode(),mtime=0))
    (args.scratch/'responses.jsonl.gz').write_bytes(gzip.compress(response_text.encode(),mtime=0))
    (args.scratch/'stderr.txt').write_text(stderr)
    replies=[json.loads(line) for line in response_text.splitlines()];assert len(replies)==len(queries)
    failures=[];maximum=dict(old_state=0.,old_potential=0.,enthalpy_identity=0.,material_delta=0.,coarse_derivative=0.,fine_derivative=0.,composition_join_limit=0.)
    arad=4*5.670374419e-5/2.99792458e10
    for i,m in enumerate(meta):
        r=replies[i];q=queries[i]
        if 'error' in r:failures.append(dict(index=i,reason=r['error']));continue
        if m['kind']=='saved':
            old=np.array(cached[tuple(q[:4])]);state=np.array(r['state'])
            maximum['old_state']=max(maximum['old_state'],float(np.max(abs(state-old[:10])/np.maximum(abs(old[:10]),1.))))
            current=np.r_[r['phi'],r['gradient'],np.asarray(r['hessian']).ravel(),r['gradient_T'],r['gradient_rho']]
            maximum['old_potential']=max(maximum['old_potential'],float(np.max(abs(current-old[21:32])/np.maximum(abs(old[21:32]),1.))))
            delta=(old[0]*old[5]-4*arad*q[2]**4/3)/(old[0]*old[6])
            h=q[2]*(delta*old[30:32]-old[28:30])
            maximum['material_delta']=max(maximum['material_delta'],abs(r['delta']/delta-1))
            maximum['enthalpy_identity']=max(maximum['enthalpy_identity'],float(np.max(abs(np.array(r['enthalpy'])-h)/np.maximum(abs(h),1.))))
        else:
            for k in range(2):
                if not q[4+k]:
                    assert r['enthalpy'][k] is None and all(v is None for v in r['enthalpy_partials'][k])
                    assert r['delta_partials'][2+k] is None
                    for row in r['enthalpy_partials']:assert row[2+k] is None
    checks=[];join_samples={}
    hydrogen=np.array(list(map(float,family.read_text().splitlines()[1].split()[2:])))
    for s in stencils:
        i,k=s['index'],s['coordinate'];q=queries[i];r=replies[i]
        parts=[replies[j] for j in s['replies']]
        if 'error' in r or any('error' in a for a in parts):
            failures.append(dict(index=i,reason='finite difference domain',stencil=s));continue
        active=[a for a in range(2) if q[4+a]]
        vec=lambda a:np.array([a['delta']]+[a['enthalpy'][j] for j in active])
        vm2,vm1,vp1,vp2=map(vec,parts)
        reference=(vm2-8*vm1+8*vp1-vp2)/(12*s['width'])
        exact=np.array([r['delta_partials'][k]]+[r['enthalpy_partials'][a][k] for a in active])*s['scale']
        norm=np.maximum(abs(vec(r)),1.)+abs(exact)
        error=float(np.max(abs(reference-exact)/norm));kind='fine_derivative' if s['width']==5e-5 else 'coarse_derivative'
        maximum[kind]=max(maximum[kind],error)
        at_join=k==2 and q[0] in hydrogen[1:-1]
        checks.append(dict(index=i,coordinate=k,width=s['width'],maximum_normalized_error=error,at_hydrogen_join=at_join))
        if at_join:
            join_samples.setdefault((i,k),{})[s['width']]=(reference,exact,norm)
        if error>1e-6 and not(at_join and s['width']==1e-4):
            failures.append(dict(reason='derivative tolerance',check=checks[-1]))
    join_checks=[]
    for (i,k),sample in join_samples.items():
        coarse,_,_=sample[1e-4];fine,exact,norm=sample[5e-5]
        # The H spline is C2. Its chemical gradient (and hence enthalpy) is
        # C1, so a symmetric first-derivative formula has an O(h) term at
        # these identified joins. Two widths remove that term independently
        # of the analytic derivative; no EOS value or tolerance is altered.
        extrapolated=2*fine-coarse
        error=float(np.max(abs(extrapolated-exact)/norm))
        maximum['composition_join_limit']=max(maximum['composition_join_limit'],error)
        join_checks.append(dict(index=i,coordinate=k,maximum_normalized_error=error,leading_difference_order=1))
        if error>1e-8:failures.append(dict(reason='composition join limit',check=join_checks[-1]))
    if max(maximum[k] for k in ['old_state','old_potential','enthalpy_identity','material_delta'])>1e-10:
        failures.append(dict(reason='retained value tolerance',maximum=maximum))
    sources=[args.probe,family,cached_input,cached_output,profile,Path(__file__),
             root/'src/eos_smooth_mixture.cpp',root/'include/ember/eos_smooth_mixture.hpp',root/'include/ember/detail/differential.hpp',
             root/'scripts/native_heat_probe.cpp',Path('/tmp/ember-native-heat-base-sources-v1/manifest.json')]
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed_native_heat_derivatives' if not failures else 'failed_native_heat_controls',
        accepted_for_stellar_evolution=False,maximum_errors=maximum,failures=failures,checks=checks,composition_join_checks=join_checks,
        native_queries=len(queries),retained_stellar_states=512,zero_species_states=len(meta)-512,derivative_states=len(controls),derivative_checks=len(checks),
        command=cmd,elapsed_seconds_including_load=elapsed,load_stderr=stderr,reused_runtime_from=str(args.reuse_raw) if args.reuse_raw else None,
        heat_evaluation_seconds=sum(r.get('seconds',0) for r in replies),new_FreeEOS_source_queries=0,
        retained_values_reused=True,base_results=[dict(query=queries[i],metadata=m,response=replies[i]) for i,m in enumerate(meta)],
        input_sha256={str(p):sha(p) for p in sources},artifacts_sha256={str(p):sha(p) for p in args.scratch.iterdir()},
        limitations=['Derivatives are of the existing smooth EOS interpolant; this does not improve underlying source physics.',
                     'Unselected species and their composition derivatives are undefined and returned as NaN/null, never floored.',
                     'No kinetic face coupling or selected stellar advance is established by these local EOS tests.'])
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['checks','base_results','input_sha256','artifacts_sha256','limitations']},indent=2))
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
