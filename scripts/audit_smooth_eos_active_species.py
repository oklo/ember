#!/usr/bin/env python3
"""Check active EOS derivatives at exact absence and resolved trace limits.

Independent SciPy composition splines and analytic mixing derivatives test
three- and four-He3-plane numerical fixtures. Five retained FreeEOS H/He4
curvature references test the physical Y3=0 endpoint without new sources.
These checks do not accept thermal diffusion or stellar abundance evolution.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
import audit_smooth_eos_runtime as three
import audit_smooth_eos_four_plane_runtime as four
from composition_potential_v2 import RGAS,ISOTOPE_MASS_RATIO,mixing_value


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def query(probe,family,rows,work,label):
    i=work/(label+'.input');o=work/(label+'.output');e=work/(label+'.stderr')
    i.write_text(''.join(' '.join(map(str,r))+'\n' for r in rows))
    with i.open('rb') as stdin,o.open('xb') as stdout,e.open('xb') as stderr:
        p=subprocess.run([str(probe),str(family)],stdin=stdin,stdout=stdout,stderr=stderr,timeout=180)
    if p.returncode:raise RuntimeError(e.read_text())
    result=[json.loads(l) for l in o.read_text().splitlines()]
    if len(result)!=len(rows):raise ValueError('missing replies')
    return result


def expected(module,oracle,row):
    x,y,T,rho,mask=row;n4=(.98-x-y)/4
    out=[None]*11
    out[0]=module.physical_jet(oracle.c(x,y),T,rho)[0,0]+mixing_value(x,y)
    # Differentiate the H-only and He3-only mixing formulas explicitly. No
    # zero population appears in the active derivative's logarithm/divisor.
    mixg={};mixh={}
    if mask&1:
        mixg[0]=RGAS*(np.log(x)+1-.25*(np.log(n4)+1))
        mixh[0,0]=RGAS*(1/x+1/(16*n4))
    if mask&2:
        mixg[1]=RGAS*((np.log(y/3)+1)/3-.25*(np.log(n4)+1)-.5*np.log(ISOTOPE_MASS_RATIO))
        mixh[1,1]=RGAS*(1/(3*y)+1/(16*n4))
    if mask==3:mixh[0,1]=mixh[1,0]=RGAS/(16*n4)
    active=[k for k in range(2) if mask&(1<<k)]
    for k in active:
        g=module.physical_jet(oracle.c(x,y,int(k==0),int(k==1)),T,rho)
        out[1+k]=g[0,0]+mixg[k];out[7+k]=g[1,0];out[9+k]=g[0,1]
        for l in active:
            h=module.physical_jet(oracle.c(x,y,int(k==0)+int(l==0),int(k==1)+int(l==1)),T,rho)
            out[3+2*k+l]=h[0,0]+mixh[k,l]
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('probe','family','work','report'):p.add_argument(n,type=Path)
    a=p.parse_args();a.work.mkdir();checks=[];guards=[];limits=[];support=[]
    for label,module in [('three',three),('four',four)]:
        directory=a.work/label;directory.mkdir();family=module.fixture(directory)
        masked=module.fixture(directory,True)
        oracle=module.Oracle();rows=[]
        # Exact absence, tiny active abundances, source nodes and both smooth
        # composition directions. Empty selection is valid even at pure He4.
        compositions=[(x,0,1) for x in (1e-16,1e-10,.0006,.07,.2,.7,.75)]
        compositions += [(0,y,2) for y in (1e-16,1e-10,.003,.03,.09,.12)]
        compositions += [(0,0,0),(.12,.06,0),(.12,.06,1),(.12,.06,2),(.12,.06,3)]
        for x,y,mask in compositions:
            for lt,lq in [(5.73,-1.4),(6.13,-.4),(6.27,1.3)]:
                T=10**lt;rows.append([x,y,T,10**lq*(T/1e6)**1.5,mask])
        for k in (0,1):
            for trace in (0.,1e-4,1e-8,1e-12,1e-16):
                rows.append([.12 if k==0 else trace,trace if k==0 else .06,1e6,1.,1<<k])
        actual=query(a.probe,family,rows,directory,'active')
        for row,r in zip(rows,actual):
            want=expected(module,oracle,row);error=0.;nulls=False
            if r['ok']:
                got=r['values'][:11];nulls=all((x is None)==(y is None) for x,y in zip(want,got))
                if nulls:error=max(abs(x-y)/max(abs(x),1.) for x,y in zip(want,got) if x is not None)
            checks.append(dict(family=label,query=row,relative_error=float(error),
                               inactive_fields_unavailable=nulls,passed=bool(r['ok'] and nulls and error<2e-6)))
        # The retained derivative block has a regular, finite limit as the
        # other species tends to zero. Its singular derivative is not used.
        for k in (0,1):
            group=actual[len(rows)-10+5*k:len(rows)-5+5*k]
            indices=[1+k,3+3*k,7+k,9+k]
            error=max(abs(group[-1]['values'][i]-group[0]['values'][i])/max(abs(group[0]['values'][i]),1.) for i in indices)
            limits.append(dict(family=label,active_direction=k,trace=1e-16,
                               relative_endpoint_difference=float(error),passed=bool(error<1e-10)))
        bad=[[0,0,1e6,1,1],[0,0,1e6,1,2],[.12,0,1e6,1,3],
             [0,.06,1e6,1,3],[-1e-10,.06,1e6,1,2],[.12,-1e-10,1e6,1,1],
             [.751,0,1e6,1,1],[0,.121,1e6,1,2],[.12,0,-1,1,1],
             [.12,0,1e6,-1,1],[.12,0,1e6,101,1],[.12,0,1e5,1,1]]
        for row,r in zip(bad,query(a.probe,family,bad,directory,'guards')):
            guards.append(dict(family=label,query=row,passed=not r['ok'],response=r))
        mask_rows=[[.0006,0,1e6,.03,1],[.0006,0,1e6,1,1],[.36,0,1e6,1,1],
                   [0,.06,1e6,.03,2],[0,.06,1e6,1,2],[0,0,1e6,1,0]]
        result=query(a.probe,masked,mask_rows,directory,'masks')
        for i,(row,r) in enumerate(zip(mask_rows,result)):
            ok=r['ok']==(i in (0,2,3))
            if r['ok']:ok=ok and (r['values'][-1]>.99 if i==2 else r['values'][-1]<.101)
            support.append(dict(family=label,query=row,passed=bool(ok),response=r))
        if label=='four':
            hm=module.fixture(directory,helium_masked=True)
            hr=[[.12,0,1e6,.03,1],[.12,0,1e6,1,1],[0,0,1e6,1,0]]
            for i,(row,r) in enumerate(zip(hr,query(a.probe,hm,hr,directory,'helium_masks'))):
                ok=r['ok']==(i==0)
                if r['ok']:ok=ok and r['values'][-1]<.101
                support.append(dict(family='four_missing_He3_plane',query=row,passed=bool(ok),response=r))
    reference=Path('docs/results/diffusion_eos_composition_curvature_v2.json')
    source=json.loads(reference.read_text());controls=source['controls']
    rows=[[v['X'],0.,v['T'],v['rho'],1] for v in controls]
    replies=query(a.probe,a.family,rows,a.work,'physical_endpoint');physical=[]
    for c,row,r in zip(controls,rows,replies):
        want=c['source_curvature_over_Rgas']*RGAS
        error=abs(r['values'][3]/want-1) if r['ok'] else None
        physical.append(dict(label=c['label'],query=row,source_curvature=want,
                             relative_error=error,passed=r['ok'] and error<.005))
    ok=all(c['passed'] for c in checks+guards+limits+support+physical)
    paths=[Path(__file__),Path(three.__file__),Path(four.__file__),
           Path('scripts/composition_potential_v2.py'),reference,
           Path('docs/results/diffusion_eos_hydrogen_spline_v1.json'),
           Path('scripts/smooth_eos_active_probe.cpp'),Path('src/eos_smooth_mixture.cpp'),
           Path('include/ember/eos_smooth_mixture.hpp'),a.probe,a.family,
           a.family.parent/'sources/freeeos300_gs98_manifest.json']
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),scope=__doc__,outcome='passed' if ok else 'failed',
                accepted_for_evolution=False,new_FreeEOS_queries=0,criteria=dict(analytic_relative=2e-6,
                trace_endpoint_relative=1e-10,physical_H_curvature_relative=.005),comparisons=checks,
                guards=guards,trace_limits=limits,source_support=support,physical_endpoint=physical,
                inputs_sha256={str(q):sha(q) for q in paths},
                artifacts_sha256={str(q):sha(q) for q in a.work.rglob('*') if q.is_file()},
                limitations=['Only selected H1/He3 directions replace He4 at fixed metals.',
                             'Physical endpoint checks cover H derivatives at Y3=0; exact-X=0 tests here are analytic.',
                             'Thermal diffusion, ionization-dependent transport and conservative stellar coupling remain required.'])
    a.report.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(outcome=report['outcome'],analytic=len(checks),guards=len(guards),
                          source_support=len(support),trace_limits=len(limits),physical=len(physical))))
    if not ok:raise SystemExit(1)


if __name__=='__main__':main()
