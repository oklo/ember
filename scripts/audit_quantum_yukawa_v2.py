#!/usr/bin/env python3
"""Check direct quantum Yukawa scattering before using its cross sections.

Numerov step, matching-radius and partial-wave controls; independent radial
Runge-Kutta solutions; weak Born and free-particle limits; direct angular
integration of the partial-wave amplitude. No stellar transport acceptance.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss,legval
from scipy.integrate import solve_ivp
from scipy.special import spherical_jn,spherical_yn

from electron_ion_born import coulomb_bracket
from quantum_yukawa import scatter_many
from fetch_tops_composition import digest


def independent_phase(a,c,ell,rmax):
    r0=1e-6/(1+c)
    h1=-c/(2*(ell+1));h2=-(a*a+c*(h1-1))/(2*(2*ell+3))
    regular=1+h1*r0+h2*r0*r0
    derivative=(ell+1)/r0+(h1+2*h2*r0)/regular
    def rhs(r,y):return [y[1],(ell*(ell+1)/r**2-a*a-c*math.exp(-r)/r)*y[0]]
    result=solve_ivp(rhs,(r0,rmax),[1.,derivative],method='DOP853',
                     rtol=1e-11,atol=1e-11,max_step=.15/max(a,1))
    if not result.success:raise ValueError(result.message)
    u,up=result.y[:,-1];x=a*rmax
    j=x*spherical_jn(ell,x);jp=a*(spherical_jn(ell,x)+x*spherical_jn(ell,x,derivative=True))
    y=x*spherical_yn(ell,x);yp=a*(spherical_yn(ell,x)+x*spherical_yn(ell,x,derivative=True))
    return math.atan2(up*j-u*jp,up*y-u*yp)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('preserve completed scattering controls')
    paths=[Path(__file__),Path('scripts/quantum_yukawa.py'),Path('scripts/quantum_yukawa_probe_v2.cpp'),
           Path('scripts/electron_ion_born.py'),a.probe]
    inputs={str(p.resolve()):digest(p) for p in paths}
    base=[(.2,.4),(1.,.5),(2.,2.),(5.,5.),(1.,20.),(5.,80.),(10.,10.)]
    queries=[]
    for k,c in base:
        lm=math.ceil(12*k+10)
        queries.extend([(k,c,.001,32.,lm),(k,c,.0005,32.,lm),
                        (k,c,.0005,24.,lm),(k,c,.0005,32.,lm+12)])
    weak=[(.2,1e-4),(1.,1e-4),(5.,1e-4)]
    for k,c in weak:queries.append((k,c,.0005,32.,math.ceil(12*k+10)))
    for k in (.2,1.,5.,10.):queries.append((k,0.,.0005,32.,math.ceil(12*k+10)))
    results=scatter_many(a.probe,queries)
    checks=[]
    def check(name,value,reference,tolerance,relative=True,**context):
        error=abs(value/reference-1) if relative else abs(value-reference)
        checks.append(dict(name=name,value=float(value),reference=float(reference),error=float(error),
                           relative=relative,tolerance=tolerance,passed=bool(np.isfinite(error) and error<tolerance),**context))
    nodes,weights=leggauss(256)
    for i,(k,c) in enumerate(base):
        coarse,fine,short,more=results[4*i:4*i+4]
        sigma=fine['cross_section_over_lambda2']
        for name,r in [('radial_step',coarse),('matching_radius',short),('partial_wave_tail',more)]:
            check(name,r['cross_section_over_lambda2'],sigma,5e-5,wavenumber=k,potential_strength=c)
        delta=np.array(fine['phase_shifts']);ell=np.arange(len(delta))
        amplitude=legval(nodes,(2*ell+1)*np.exp(1j*delta)*np.sin(delta))/k
        angular=2*math.pi*np.sum(weights*(1-nodes)*np.abs(amplitude)**2)
        check('angular_transport_sum',sigma,angular,1e-7,wavenumber=k,potential_strength=c)
        for l in (0,1,3):
            reference=independent_phase(k,c,l,32.)
            difference=math.asin(math.sin(delta[l]-reference))
            check('independent_radial_ODE',difference,0.,2e-6,relative=False,
                  wavenumber=k,potential_strength=c,ell=l)
        print(json.dumps({'converged_case_checked':[k,c]}),flush=True)
    start=4*len(base)
    for i,(k,c) in enumerate(weak):
        sigma=results[start+i]['cross_section_over_lambda2']
        born=math.pi*c*c/(2*k**4)*float(coulomb_bracket(4*k*k))
        check('weak_Born_limit',sigma,born,1e-3,wavenumber=k,potential_strength=c)
    for r in results[start+len(weak):]:
        check('free_particle',r['cross_section_over_lambda2'],0.,1e-10,relative=False,wavenumber=r['wavenumber'])
    for path,h in inputs.items():
        if digest(Path(path))!=h:raise ValueError('quantum scattering source changed')
    passed=all(c['passed'] for c in checks)
    report=dict(scope=__doc__,status='pass' if passed else 'fail',accepted_for_stellar_evolution=False,
                checks=checks,numerical_results=results,input_sha256=inputs,
                limitations=['A prescribed attractive Yukawa potential, stationary ions and nonrelativistic electrons.',
                             'No ion structure factor, dynamic screening, finite-temperature drag integral or stellar velocities.',
                             'Numerical convergence does not establish the physical accuracy of the selected interaction.'])
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=report['status'],checks=len(checks),failures=[c for c in checks if not c['passed']])),flush=True)
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
