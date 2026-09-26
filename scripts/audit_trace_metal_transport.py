#!/usr/bin/env python3
"""Analytic finite-volume checks and the trace approximation at real bases."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import numpy as np
from trace_metal_transport import ELEMENTS,coefficients,bernoulli,step
from prepare_nongrey_sources import digest
from write_scientific_result import write_result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('output',type=Path);args=ap.parse_args()
    if args.output.exists():raise FileExistsError('preserve reports')
    source=Path('docs/results/thoul_diffusion_assessment_v1.json')
    reference=json.loads(source.read_text())
    inputs={str(p.resolve()):digest(p) for p in [Path(__file__),Path('scripts/trace_metal_transport.py'),
            Path('scripts/diffusion_burgers.py'),Path('scripts/ion_collision_integrals.py'),source]}
    checks=[]
    def check(name,error,tolerance):
        assert np.isfinite(error) and error<=tolerance,(name,error)
        checks.append(dict(name=name,error=float(error),tolerance=tolerance))
    for dt in [1e-4,1.,1e4]:
        w=np.array([.3,.7]);x=np.array([.9,.1]);f=.2;b=.4
        total=w@x;eq0=total/(w[0]+w[1]*f/b)
        expected0=eq0+(x[0]-eq0)/(1+dt*(f/w[0]+b/w[1]))
        expected=np.array([expected0,(total-w[0]*expected0)/w[1]])
        answer,error=step(x,w,[f],[b],[(0,1),(1,2)],dt)
        check('two-cell backward-Euler analytic exchange',max(abs(answer-expected)),1e-12)
        check('two-cell conserved inventory',error,1e-12)
    w=np.linspace(1.,2.,21);x=np.exp(np.linspace(-4.,4.,21));pe=np.full(20,.4)
    for dt in [1.,1e6]:
        answer,error=step(x,w,bernoulli(-pe),bernoulli(pe),[(i,i+1) for i in range(21)],dt)
        check('zero-flux exponential equilibrium',max(abs(answer/x-1)),1e-9)
    # New convective reservoir combines unequal masses. Retreating the
    # boundary subsequently releases cells with exactly their mixed value.
    x=np.linspace(.01,1.,21)
    mixed,_=step(x,w,np.zeros(20),np.zeros(20),[(i,i+1) for i in range(8)]+[(8,21)],1.)
    check('mass-weighted convection',max(abs(mixed[8:]-(w[8:]@x[8:])/sum(w[8:]))),1e-14)
    released,_=step(mixed,w,np.zeros(20),np.zeros(20),[(i,i+1) for i in range(13)]+[(13,21)],1.)
    check('retreating mixed boundary preserves inventory and abundance',max(abs(released-mixed)),1e-14)
    for pe in [-1000.,-1e-8,0.,1e-8,1000.]:
        f=float(bernoulli(np.array(-pe)));b=float(bernoulli(np.array(pe)))
        answer,error=step(np.array([0.,1.]),np.array([.3,.7]),[f],[b],[(0,1),(1,2)],1e4)
        check('large-drift positivity',max(0.,-min(answer)),0.)
    controls=[]
    for r in reference['records']:
        if not r['envelope_base']:continue
        p=Path(f"/tmp/ember-cn-main-sequence-v{r['version']}/checkpoint.json")
        inputs[str(p)]=digest(p);m=np.array(json.loads(p.read_text())['model_record']['model'])
        face=r['face']
        for name,A,Z in ELEMENTS:
            index=reference['species'].index(name)
            for heat in ['classical','suppressed']:
                small=np.array(coefficients(m[face],m[face+1],A,Z,heat,1e-10))
                larger=np.array(coefficients(m[face],m[face+1],A,Z,heat,1e-8))
                check('trace population limit',max(abs(larger/small-1)),1e-5)
                full=r['actual_gravity_controls'][heat]['velocity_cm_s'][index]
                controls.append(dict(version=r['version'],element=name,heat_flow=heat,
                    trace_velocity_cm_s=float(small[0]),diffusivity_cm2_s=float(small[1]),
                    trace_over_full_mixture=float(small[0]/full)))
    assert all(digest(p)==value for p,value in inputs.items())
    result=dict(utc=datetime.now(timezone.utc).isoformat(),outcome='trace_transport_checks_passed',
                selected=False,checks=checks,mixture_controls=controls,input_sha256=inputs,
                limitations=['Trace/full differences include neglect of other metals and H/He renormalization.',
                             'Both comparisons share published classical collision coefficients and ionization assumptions.'])
    write_result(args.output,result)
    print(json.dumps(dict(outcome=result['outcome'],checks=len(checks),
        trace_full_ratio_range=[min(x['trace_over_full_mixture'] for x in controls),
                               max(x['trace_over_full_mixture'] for x in controls)]),indent=2))


if __name__=='__main__':main()
