#!/usr/bin/env python3
"""Compare analytic collision derivatives with independent scalar differences.

Replays the retained scalar physical controls, differentiates every valid
state, and checks selected saved and changed states at two difference widths.
Spline knots and the ion-extension join are included explicitly. Raw replies
are retained losslessly; no EOS or pair-scattering integrals are repeated.
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
from electron_ion_born import fermi_half, KB, ME, HBAR, E2
from ion_collision_integrals import MU
from audit_diffusion_material_regime import diagnostics


def vector(r):
    return np.r_[np.asarray(r['mobility']).ravel(),r['enthalpy'],
                 [r[k] for k in ['conductivity','energy_scale','eta','ne','b']]]


def scales(r):
    diagonal=np.sqrt(np.diag(r['mobility']))
    matrix=np.outer(diagonal,diagonal).ravel()
    matrix[matrix==0]=1
    return np.r_[matrix,np.maximum(abs(np.array(r['enthalpy'])),1.),
                 [max(abs(r[k]),1.) for k in ['conductivity','energy_scale','eta','ne','b']]]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('probe',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--scratch',type=Path,required=True)
    args=p.parse_args();assert not args.output.exists() and not args.scratch.exists()
    args.scratch.mkdir();root=Path(__file__).resolve().parents[1]
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    prior_path=root/'docs/results/collision_transport_cpp_v2.json'
    prior=json.loads(prior_path.read_text());assert prior['outcome']=='passed_physical_evaluator_comparisons'
    snapshot_path=Path('/tmp/ember-collision-transport-source-v1/manifest.json')
    snapshots=json.loads(snapshot_path.read_text())
    for name,h in prior['input_sha256'].items():
        if sha(name)!=h:
            assert snapshots[name]['sha256']==h and sha(snapshots[name]['path'])==h,name
    for name,h in prior['raw_artifacts'].items():assert sha(name)==h
    table=Path('/tmp/ember-collision-transport-table-v1.dat')
    assert sha(table)==prior['input_sha256'][str(table)]
    raw=Path('/tmp/ember-collision-transport-audit-v2/queries.txt').read_text()
    queries=[list(map(float,s.split())) for s in raw.splitlines()]
    meta=[dict(kind=r['kind'],zone=r.get('zone')) for r in prior['records']]
    # Independent construction of states at every interior pair-spline knot.
    # ne follows the Fermi number integral, rho the baryonic mixture count.
    material=json.loads((root/'docs/results/diffusion_material_regime_3890gyr_v1.json').read_text())
    composition=[.05,.001,.02]
    mr=material['records'][300]
    ye=mr['electron_density_cm3']*MU/mr['density_baryonic_g_cm3']
    h=mr['X']+2*mr['Y3']/3+mr['Y4']/2
    z=1-mr['X']-mr['Y3']-mr['Y4'];metal_e=(ye-h)/z
    def state(eta,B,T=1e7):
        X,Y3,Z=composition
        ne=(2*ME*KB*T)**1.5/(2*math.pi**2*HBAR**3)*fermi_half(eta)
        rho=ne*MU/(X+2*Y3/3+(1-X-Y3-Z)/2+Z*metal_e)
        length=math.sqrt(B*HBAR**2/(8*ME*KB*T))
        return [T,rho,X,Y3,Z,length]
    for eta in [-1.,0.,1.,2.,3.,4.,5.5]:
        queries.append(state(eta,16.));meta.append(dict(kind='eta_spline_knot',eta=eta))
    for B in [7.,16.,40.,100.]:
        queries.append(state(2.,B));meta.append(dict(kind='screening_spline_knot',B=B))
    # Pick Ni-Ni to reach the original ion-table endpoint at a plausible T.
    export=json.loads(table.with_suffix('.json').read_text())
    extension=json.loads(Path('/tmp/ember-ion-table-extension-run-v1/extended_collision_integrals.json').read_text())
    base=json.loads(Path(extension['base_table']).read_text())
    strength=10**base['log10_strength'][-1]
    B=16.;interaction=28*28*E2
    T=8*ME*interaction**2/(KB*B*HBAR**2*strength**2)
    queries.append(state(2.,B,T));meta.append(dict(kind='ion_extension_join',pair_charge_product=784.,strength=strength))
    jobs=[]
    def run(name,request,derivatives=False):
        text=''.join(' '.join(format(x,'.17g') for x in q)+'\n' for q in request)
        cmd=[str(args.probe),str(table)]+(['derivatives'] if derivatives else [])
        start=time.monotonic()
        proc=subprocess.run(cmd,input=text,capture_output=True,text=True,timeout=90,check=True)
        elapsed=time.monotonic()-start
        for suffix,content in [('queries.txt.gz',text),('responses.jsonl.gz',proc.stdout)]:
            (args.scratch/(name+'.'+suffix)).write_bytes(gzip.compress(content.encode(),mtime=0))
        replies=[json.loads(s) for s in proc.stdout.splitlines()]
        assert len(replies)==len(request)
        jobs.append(dict(name=name,command=cmd,queries=len(request),elapsed_seconds=elapsed,
                         evaluation_seconds=sum(v.get('seconds',0) for v in replies)))
        return replies
    scalar=run('scalar',queries)
    dual=run('derivatives',queries,True)
    failures=[];maxima=dict(scalar_replay=0.,derivative_value=0.,fine_derivative=0.,coarse_derivative=0.,screening=0.)
    for i,(s,d) in enumerate(zip(scalar,dual,strict=True)):
        if 'error' in s or 'error' in d:
            if i>=len(prior['records']) or 'error' not in s or 'error' not in d or prior['records'][i]['kind']!='invalid':
                failures.append(dict(index=i,kind='unexpected rejection',scalar=s,derivative=d))
            continue
        if i<len(prior['records']):
            old=prior['records'][i]['cpp']
            error=float(np.max(abs(vector(s)-vector(old))/scales(old)))
            maxima['scalar_replay']=max(maxima['scalar_replay'],error)
        error=float(np.max(abs(vector(s)-vector(d))/scales(s)))
        maxima['derivative_value']=max(maxima['derivative_value'],error)
        expected=[True,True,queries[i][2]>0,queries[i][3]>0,queries[i][4]>0,True]
        assert d['defined']==expected
        for j,part in enumerate(d['partials']):
            assert np.isfinite(vector(part)).all()
            assert np.array_equal(np.asarray(part['mobility']),np.asarray(part['mobility']).T)
            if not expected[j]:assert np.count_nonzero(vector(part))==0
    indices=[i for i,m in enumerate(meta) if m['kind'] not in ['saved','invalid']]
    indices += [2*z+k for z in [0,40,100,160,220,300,360,396,410,421] for k in [0,1]]
    perturbations=[];stencils=[]
    for i in indices:
        if 'error' in dual[i]:continue
        q=queries[i];d=dual[i]
        for k,defined in enumerate(d['defined']):
            if not defined:continue
            # Fraction derivatives are measured in the corresponding relative
            # change, resolving a tiny trace species without crossing zero.
            coordinate_scale=q[k] if 2<=k<=4 else 1.
            for width in [2e-3,1e-3]:
                stencil=[]
                for step in [-2,-1,1,2]:
                    trial=q.copy()
                    if k in [0,1,5]:
                        target=6 if k==5 and len(q)==8 else k
                        trial[target]*=math.exp(step*width)
                    else:trial[k]+=step*width*coordinate_scale
                    stencil.append(len(perturbations));perturbations.append(trial)
                stencils.append(dict(index=i,coordinate=k,scale=coordinate_scale,width=width,replies=stencil))
    finite=run('finite_differences',perturbations)
    checks=[]
    for stencil in stencils:
        i,k=stencil['index'],stencil['coordinate'];d=dual[i];q=queries[i]
        rows=[finite[j] for j in stencil['replies']]
        if any('error' in r for r in rows):
            failures.append(dict(kind='difference-domain',stencil=stencil,rows=rows));continue
        vm2,vm1,vp1,vp2=[vector(r) for r in rows]
        reference=(vm2-8*vm1+8*vp1-vp2)/(12*stencil['width'])
        exact=vector(d['partials'][k]).copy()
        if len(q)==8:
            if k==5:exact*=0
            exact+=vector(d['partials'][5])*d['screening_partials'][k]/d['screening_length']
        exact*=stencil['scale']
        normalization=scales(d)+abs(exact)
        errors=abs(reference-exact)/normalization
        error=float(np.max(errors));key='fine_derivative' if stencil['width']==1e-3 else 'coarse_derivative'
        maxima[key]=max(maxima[key],error)
        check=dict(index=i,kind=meta[i]['kind'],coordinate=k,width=stencil['width'],maximum_normalized_error=error,
                   worst_component=int(np.argmax(errors)))
        if len(q)==8:
            lm2,lm1,lp1,lp2=[r['screening_length'] for r in rows]
            reference_length=(lm2-8*lm1+8*lp1-lp2)/(12*stencil['width'])
            exact_length=d['screening_partials'][k]*stencil['scale']
            se=abs(reference_length-exact_length)/(d['screening_length']+abs(exact_length))
            check['screening_error']=se;maxima['screening']=max(maxima['screening'],se)
        checks.append(check)
        if error>1e-7:failures.append(dict(kind='derivative-tolerance',check=check))
    if maxima['scalar_replay']>2e-10 or maxima['derivative_value']>2e-10 or maxima['screening']>1e-8:
        failures.append(dict(kind='value-or-screening-tolerance',maxima=maxima))
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),
        outcome='passed_analytic_derivative_controls' if not failures else 'failed_derivative_controls',
        accepted_for_stellar_evolution=False,maximum_errors=maxima,failures=failures,checks=checks,
        state_count=len(queries),original_physical_states_replayed=len(prior['records']),
        differentiated_valid_states=sum('error' not in d for d in dual),finite_difference_states=len(indices),
        finite_difference_queries=len(perturbations),derivative_checks=len(checks),jobs=jobs,
        selected_states=[dict(index=i,metadata=meta[i],query=queries[i],derivative=dual[i]) for i in indices],
        new_EOS_queries=0,new_pair_collision_integrals=0,
        normalization='Component difference divided by the mobility diagonal scale or scalar value scale plus derivative magnitude. Fraction partials use relative composition increments.',
        input_sha256={str(p):sha(p) for p in [args.probe,table,prior_path,snapshot_path,Path(__file__),
            root/'src/collision_transport.cpp',root/'include/ember/collision_transport.hpp',root/'src/differential.hpp',
            root/'scripts/collision_transport_probe.cpp',root/'include/ember/gs98_mixture.hpp',
            root/'scripts/electron_ion_born.py',table.with_suffix('.json')]},
        artifacts_sha256={str(p):sha(p) for p in args.scratch.iterdir()},
        limitations=['Derivatives describe the conditional hot kinetic model, not an ionization or correlation prescription.',
                     'At an absent species the derivative for introducing that species is explicitly undefined; no floor is used.',
                     'The physical stellar interface still requires native EOS force/enthalpy derivatives and partial-ionization coverage.'])
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['checks','selected_states','input_sha256','artifacts_sha256','limitations']},indent=2))
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
