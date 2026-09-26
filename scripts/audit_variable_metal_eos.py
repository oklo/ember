#!/usr/bin/env python3
"""Check variable-metal material derivatives and independent FreeEOS values.

The finite-difference checks establish consistency of the interpolated
potential; held-out source compositions separately measure interpolation
error. Neither establishes full stellar coverage or selects the new EOS.
"""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time
import numpy as np
from metal_eos_composition import mixture
from write_scientific_result import write_result

R = 8.31446261815324e7
ARAD = 4*5.670374419e-5/2.99792458e10


def isotope_mixing(x, y, z):
    n3, n4 = y/3, (1-math.fsum((x, y, z)))/4
    xl = lambda n: n*math.log(n) if n > 0 else 0.
    # nuclides[].A are atomic masses, matching the selected isotope correction.
    return R*(xl(n3)+xl(n4)-xl(n3+n4)-1.5*n3*math.log(3.01602932/4.00260325))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('probe', type=Path); p.add_argument('family', type=Path)
    p.add_argument('source_probe', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--scratch', type=Path, required=True)
    p.add_argument('--center-logT', type=float, default=6.5)
    p.add_argument('--center-logQ', type=float, default=.1)
    a = p.parse_args()
    if a.output.exists() or a.scratch.exists(): raise FileExistsError('preserve completed controls')
    a.scratch.mkdir(parents=True)
    sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
    root = Path(__file__).resolve().parents[1]
    paths = [a.probe, a.family, a.source_probe, Path(__file__), root/'scripts/metal_eos_composition.py',
             root/'src/eos_variable_metal.cpp', root/'src/eos_helmholtz.cpp',
             root/'src/composition_spline.hpp', root/'include/ember/eos_variable_metal.hpp',
             root/'include/ember/composition.hpp', root/'scripts/variable_metal_eos_probe.cpp']
    hashes = {str(path.resolve()): sha(path) for path in paths}
    queries, bases, stencils, heldouts, limits = [], [], [], [], []
    def add(q): queries.append(q); return len(queries)-1
    def query(x, y, z, t=6.5125, q=.1125):
        t += a.center_logT-6.5; q += a.center_logQ-.1
        return [x, y, z, 10**t, 10**(q+1.5*(t-6)), int(x>0), int(y>0), int(z>0)]
    # Interior points, shrinking helium reservoir, and inactive-species faces.
    for x, y, z in [(.73,.008,.007),(.88,.0006,.013),(.98,.0002,.003),
                    (.994,.00005,.001),(.001,.07,.025),(0,.03,.01),(.8,0,.01),(.7,.02,0)]:
        q = query(x,y,z); i = add(q); bases.append(i)
        helium = 1-math.fsum(q[:3])
        for k in range(5):
            if k>=2 and not q[5+k-2]: continue
            scale = 1. if k<2 else min(q[k-2], helium)
            for width in [1e-4,5e-5]:
                ids = []
                for step in [-2,-1,1,2]:
                    v = q.copy()
                    if k<2: v[k+3] *= math.exp(step*width)
                    else: v[k-2] += step*width*scale
                    ids.append(add(v))
                stencils.append(dict(index=i, coordinate=k, scale=scale, width=width, indices=ids))
    # Sample independent compositions and temperatures; none are fit inputs.
    for x,y,z in [(.73,.008,.007),(.88,.0006,.013),(.98,.0002,.003),
                  (.994,.00005,.001),(.001,.07,.025),(.42,.04,.035),
                  (0,.03,.01),(.8,0,.01),(.7,.02,0)]:
        for t,qv in [(6.4875,.0875),(6.5125,.1125)]:
            heldouts.append(add(query(x,y,z,t,qv)))
    for x,y,z in [(1,0,0),(.999,0,.001),(0,0,0),(0,0,.02)]:
        q=query(x,y,z);q[5:]=[0,0,0];limits.append(dict(index=add(q),accepted=True))
    for x,y,z,active in [(.98,.02,.02,[1,1,1]),(.85,.06,.01,[1,1,1]),
                         (.7,.01,.06,[1,1,1]),(-.01,.01,.02,[0,1,1]),
                         (1,0,0,[1,0,0]),(.8,0,.01,[1,1,1]),(.8,.01,0,[1,1,1])]:
        q=query(x,y,z);q[5:]=active;limits.append(dict(index=add(q),accepted=False))
    q=query(.8,.001,.01,t=6.7); limits.append(dict(index=add(q),accepted=False))
    text=''.join(' '.join(format(float(v),'.17g') for v in q)+'\n' for q in queries)
    start=time.monotonic()
    run=subprocess.run([str(a.probe),str(a.family)],input=text,text=True,capture_output=True,check=True,timeout=60)
    elapsed=time.monotonic()-start
    (a.scratch/'queries.txt.gz').write_bytes(gzip.compress(text.encode(),mtime=0))
    (a.scratch/'responses.jsonl.gz').write_bytes(gzip.compress(run.stdout.encode(),mtime=0))
    (a.scratch/'stderr.txt').write_text(run.stderr)
    replies=[json.loads(line) for line in run.stdout.splitlines()]
    if len(replies)!=len(queries): raise ValueError('missing native replies')
    failures=[]; maximum=dict(potential_derivatives=0.,heat_derivatives=0.,enthalpy_identity=0.,
                              source_state=0.,source_potential=0.,source_gradient=0.,source_enthalpy=0.)
    def record(kind,value,limit,**context):
        maximum[kind]=max(maximum[kind],float(value))
        if value>limit: failures.append(dict(kind=kind,error=float(value),limit=limit,**context))
    for i in bases:
        q,r=queries[i],replies[i]
        if 'error' in r: failures.append(dict(index=i,error=r['error'])); continue
        for k,active in enumerate(q[5:]):
            if not active:
                if r['gradient'][k] is not None or r['enthalpy'][k] is not None:
                    failures.append(dict(index=i,kind='inactive_species_returned_value'))
                continue
            expected=q[3]*(r['delta']*r['gradient_rho'][k]-r['gradient_T'][k])
            record('enthalpy_identity',abs(r['enthalpy'][k]-expected)/max(abs(expected),1),1e-12,index=i)
    derivative_checks=[]
    for s in stencils:
        i,k=s['index'],s['coordinate'];q=queries[i];r=replies[i];parts=[replies[j] for j in s['indices']]
        if 'error' in r or any('error' in v for v in parts):
            failures.append(dict(kind='stencil_domain',stencil=s));continue
        active=[j for j in range(3) if q[5+j]]
        def vec(v):
            return np.array([v['phi']]+[v['gradient'][j] for j in active]+[v['delta']]
                            +[v['enthalpy'][j] for j in active]+[v['radiation_enthalpy'][j] for j in active])
        fm2,fm1,fp1,fp2=map(vec,parts)
        reference=(fm2-8*fm1+8*fp1-fp2)/(12*s['width'])
        if k<2:
            # Material phi_T=-E_material/T; phi_rho=P_material/(rho*T).
            pr=ARAD*q[3]**4/3
            first=-(r['state'][1]-3*pr/q[4])/q[3] if k==0 else (r['state'][0]-pr)/(q[4]*q[3])
            grad=r['gradient_T'] if k==0 else r['gradient_rho']
        else: first=r['gradient'][k-2];grad=[r['hessian'][j][k-2] for j in range(3)]
        exact=np.array([first]+[grad[j] for j in active]+[r['delta_partials'][k]]
                       +[r['enthalpy_partials'][j][k] for j in active]
                       +[r['radiation_enthalpy_partials'][j][k] for j in active])*s['scale']
        norm=np.maximum(abs(vec(r)),1)+abs(exact)
        errors=abs(reference-exact)/norm
        pe,he=float(max(errors[:1+len(active)])),float(max(errors[1+len(active):]))
        record('potential_derivatives',pe,2e-6,index=i,coordinate=k,width=s['width'])
        record('heat_derivatives',he,2e-6,index=i,coordinate=k,width=s['width'])
        derivative_checks.append(dict(**s,potential_error=pe,heat_error=he))
    for c in limits:
        accepted='error' not in replies[c['index']]
        if accepted!=c['accepted']: failures.append(dict(kind='domain_control',**c,response=replies[c['index']]))
    source_records=[]
    def source(q):
        x,y,z,T,rho=q[:5];m=mixture(x,y,metallicity=z);scale=m['source_mass_scale']
        request=' '.join(format(v,'.17g') for v in m['eps'])+'\n3 223 -2\n'+f'{math.log(rho*scale):.17g} {math.log(T):.17g}\n'
        run=subprocess.run([str(a.source_probe)],input=request,text=True,capture_output=True,check=True,timeout=60)
        row=np.array([float(v) for v in run.stdout.split()])
        if len(row)!=22 or row[0]!=0 or not np.isfinite(row).all(): raise ValueError('failed direct source')
        source_records.append(dict(query=q[:5],request=request,raw=row.tolist()))
        P,E,S=row[4],row[5]*scale,row[6]*scale-isotope_mixing(x,y,z)
        phi=E/T-S;pr=ARAD*T**4/3;Ptot=P+pr
        cr=P*row[7]/Ptot;ct=(P*row[8]+4*pr)/Ptot;delta=ct/cr
        cv=row[10]*scale/T+12*pr/(rho*T);cp=cv+Ptot/(rho*T)*ct*delta
        state=np.array([Ptot,E+3*pr/rho,S+4*pr/(rho*T),cv,cp,ct,cr,delta,Ptot/(rho*T)*delta/cp,cr*cp/cv])
        return dict(state=state,phi=phi,phi_rho=P/(rho*T),phi_T=-E/T,delta=row[15])
    source_checks=[]
    for i in heldouts+ [c['index'] for c in limits if c['accepted']]:
        q,r=queries[i],replies[i]
        if 'error' in r: failures.append(dict(kind='source_control_domain',index=i,error=r['error']));continue
        ref=source(q);state=np.array(r['state']);error=float(max(abs(state-ref['state'])/np.maximum(abs(ref['state']),1)))
        record('source_state',error,1e-3,index=i)
        record('source_potential',abs(r['phi']-ref['phi'])/max(abs(ref['phi']),1),1e-3,index=i)
        source_checks.append(dict(index=i,query=q[:5],state_error=error,native=state.tolist(),source=ref['state'].tolist()))
        if i not in heldouts:continue
        helium=1-math.fsum(q[:3])
        for k in range(3):
            if not q[5+k]: continue
            h=1e-3*min(q[k],helium);values=[]
            for step in [-2,-1,1,2]:
                trial=q.copy();trial[k]+=step*h;values.append(source(trial))
            deriv=lambda name:sum(w*v[name] for w,v in zip([1,-8,8,-1],values))/(12*h)
            gradient=deriv('phi');enthalpy=q[3]*(ref['delta']*deriv('phi_rho')-deriv('phi_T'))
            record('source_gradient',abs(r['gradient'][k]-gradient)/max(abs(gradient),R),3e-3,index=i,species=k)
            record('source_enthalpy',abs(r['enthalpy'][k]-enthalpy)/max(abs(enthalpy),R*q[3]),3e-3,index=i,species=k)
    (a.scratch/'direct_sources.json.gz').write_bytes(gzip.compress(json.dumps(source_records,allow_nan=False).encode(),mtime=0))
    for path,expected in hashes.items():
        if sha(path)!=expected: raise ValueError(f'input changed: {path}')
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),
                outcome='passed_local_variable_metal_eos_controls' if not failures else 'failed_local_variable_metal_eos_controls',
                selected_for_evolution=False,maximum_errors=maximum,failures=failures,
                native_queries=len(queries),source_queries=len(source_records),native_elapsed_seconds=elapsed,
                derivative_checks=derivative_checks,source_checks=source_checks,domain_checks=limits,
                patch_center=dict(logT=a.center_logT,logQ=a.center_logQ),
                input_sha256=hashes,artifacts_sha256={str(p):sha(p) for p in a.scratch.iterdir()},
                limitations=['Small thermal patch only; the source and derivative checks do not establish global coverage.',
                             'Fixed GS98 relative pattern; independent heavy-element changes are not represented.',
                             'Fully coupled transport, atmosphere response, and global stellar coverage remain separate work.'])
    write_result(a.output,report)
    print(json.dumps({k:v for k,v in report.items() if k in ['outcome','maximum_errors','failures','native_queries','source_queries','native_elapsed_seconds']},indent=2))
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
