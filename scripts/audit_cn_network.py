#!/usr/bin/env python3
"""Independent native CN rates, source equations and stiff catalyst-update checks."""
import argparse
import hashlib
import importlib.util
import json
import math
import subprocess
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np


def independent_implicit_solve(old, rate, dt):
    with localcontext() as context:
        context.prec = 320
        x,y,z = [Decimal(v) for v in rate]
        h = Decimal(dt)
        a = [[1+h*x,Decimal(0),-h*z,Decimal(old[0])],
             [-h*x,1+h*y,Decimal(0),Decimal(old[1])],
             [Decimal(0),-h*y,1+h*z,Decimal(old[2])]]
        for k in range(3):
            pivot=max(range(k,3),key=lambda j:abs(a[j][k]))
            a[k],a[pivot]=a[pivot],a[k]
            for i in range(k+1,3):
                ratio=a[i][k]/a[k][k]
                for j in range(k,4):
                    a[i][j]-=ratio*a[k][j]
        answer=[Decimal(0)]*3
        for k in range(2,-1,-1):
            answer[k]=(a[k][3]-sum(a[k][j]*answer[j] for j in range(k+1,3)))/a[k][k]
        return [float(v) for v in answer]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--native', type=Path, required=True)
    ap.add_argument('--scratch', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    a.scratch.mkdir()
    if a.output.exists():
        raise ValueError('preserve existing audit')
    worker = root/'docs/research/fable/results/fable-restart-approximation-v1/cno_rate_v2.py'
    spec = importlib.util.spec_from_file_location('reference', worker)
    ref = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ref)
    identities = {}
    for p in [a.native, worker, Path(__file__), root/'scripts/cn_network_probe.cpp',
              root/'src/nuclear_cn.cpp', root/'src/nuclear_pp.cpp', root/'include/ember/nuclear_cn.hpp',
              root/'include/ember/nuclear.hpp', root/'include/ember/composition.hpp', root/'include/ember/gs98_mixture.hpp']:
        identities[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
    initial = np.array([.02*.171836/12, 0., .02*.050335/14])
    y = initial+initial[0]*np.array([-.6, .15, .45])
    queries, cases = [], []
    for rates in (2, 3):
        for T, rho, X, X3 in ((5e6,300.,.5,.003),(8e6,300.,.2,.001),(1.59e7,8500.,.00055,1e-9)):
            for screen in (0, 1):
                q = [T,rho,X,X3,*y,rates,screen]
                queries.append('source '+' '.join(format(v,'.17g') for v in q))
                cases.append(q)
    implicit = []
    for old in ([1.,0.,0.],[0.,1.,0.],[0.,0.,1.],list(initial)):
        for rate in ([1.,4.,.001],[1.,0.,0.],[0.,0.,1.],[0.,0.,0.],[1e-20,1e20,1e5]):
            for dt in (0.,1e-10,1.,1e5,1e30,1e100):
                implicit.append((old,rate,dt))
                queries.append('implicit '+' '.join(format(v,'.17g') for v in [dt,*rate,*old]))
    raw='\n'.join(queries)+'\n'
    (a.scratch/'queries.txt').write_text(raw)
    run=subprocess.run([str(a.native)],input=raw,capture_output=True,text=True,timeout=120,check=True)
    (a.scratch/'responses.jsonl').write_text(run.stdout)
    (a.scratch/'stderr.txt').write_text(run.stderr)
    replies=[json.loads(line) for line in run.stdout.splitlines()]
    assert len(replies)==len(queries) and not any('error' in r for r in replies)
    comparisons=[]
    names=('C12(p,g)','C13(p,g)','N14(p,g)')
    # Explicit physical changes per reaction, independently written in mass fractions.
    stoich=np.array([[-1,0,0,-12,13,0,0,0],[-1,0,0,0,-13,14,0,0],[-2,0,4,12,0,-14,0,0]],float).T
    atomic=np.array([1.00782503,3.01602932,4.00260325,12.,13.00335484,14.0030740,15.9949146,20.])
    mass=np.array([1,3,4,12,13,14,16,20],float)
    q=-np.dot(atomic/mass-1,stoich)*ref.C_L**2
    nu=np.array([.706,0,.996])*ref.MEV*ref.NA
    for case, actual in zip(cases,replies):
        T,rho,X,X3,*_=case; rates,screen=case[-2:]
        table=ref.SF2 if rates==2 else ref.SF3
        bare=[];slopes=[];screens=[]
        for name in names:
            r,s=ref.bare_rate(T,table[name],npts=24001)
            bare.append(r);slopes.append(s)
            screens.append(ref.screening(T,rho,X,X3,table[name],model='svh' if screen else 'weak',basis='baryon_mass')['logf'])
        flow=rho*X*np.array(bare)*np.exp(screens)*y
        expected={'bare':bare,'slope':slopes,'screen':screens,'reaction_rate':flow.tolist(),
                  'dXdt':(stoich@flow).tolist(),'eps':float(np.dot(q-nu,flow)),
                  'eps_neutrino':float(np.dot(nu,flow))}
        errors={}
        for name, value in expected.items():
            scale=np.maximum(np.abs(value),1e-100)
            errors[name]=float(np.max(np.abs(np.array(actual[name])-value)/scale))
        comparisons.append({'query':case,'expected':expected,'native':actual,'relative_errors':errors})
    maximum=max(e for row in comparisons for e in row['relative_errors'].values())
    assert maximum<2e-8, maximum
    # Independent high-precision solve of (I-dt A)Y_new = Y_old. The native
    # positive-adjugate formula is not reproduced by this reference.
    implicit_rows=[]
    for (old,rate,dt),actual in zip(implicit,replies[len(cases):]):
        expected=independent_implicit_solve(old,rate,dt)
        error=max(abs(u-v) for u,v in zip(actual['molality'],expected))/sum(old)
        conservation=abs(sum(actual['molality'])/sum(old)-1)
        assert error<2e-14 and conservation<2e-14 and min(actual['molality'])>=0
        implicit_rows.append({'old':old,'frequency':rate,'dt':dt,'normalized_error':error,'conservation_error':conservation})
    for p in a.scratch.iterdir():
        if p.is_file():identities[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    report={'created_utc':datetime.now(timezone.utc).isoformat(),'outcome':'passed',
            'source_comparisons':comparisons,'maximum_source_relative_error':maximum,
            'implicit_comparisons':implicit_rows,'maximum_implicit_normalized_error':max(r['normalized_error'] for r in implicit_rows),
            'input_sha256':identities,'limitations':[
                'CN-only reaction network; oxygen branches and catalyst diffusion are absent.',
                'Screening uses the supplied bulk composition, independently of the catalyst molalities.',
                'The implicit kernel holds end-of-step capture frequencies fixed; coupling to fuel, mixing and structure remains required.',
                'Short beta decays and N15 proton capture are eliminated; their timescale separation must remain valid.']}
    a.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['outcome','maximum_source_relative_error','maximum_implicit_normalized_error']},indent=2))


if __name__=='__main__':
    main()
