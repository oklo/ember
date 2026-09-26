#!/usr/bin/env python3
"""Compare scattering quadrature with direct Cartesian two-body trajectories.

Reduced units: relative incoming speed 1, repulsive Yukawa acceleration
exp(-r)*(1+r)*r_vector/(2*w^2*r^3). Infinite-impact initial conditions are
approximated at a recorded finite separation, with two distances checked.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.optimize import brentq

from fetch_tops_composition import digest
from yukawa_scattering_control import gauss,cross_sections
from ion_collision_integrals import cross_section_fit


def angle_quadrature(w,b):
    r0=brentq(lambda r:r*r-b*b-r*math.exp(-r)/(w*w),b,b+40,
              xtol=1e-14,rtol=1e-14)
    # Independent direct form after u=r0/r=1-s^2, without Coulomb subtraction.
    gamma=b/r0
    def fn(s):
        if s==0:return 2/math.sqrt(2*gamma*gamma+(1-gamma*gamma)*(1+r0))
        u=1-s*s
        if u<=0:return 2*s
        d=(1-u*u)+(1-gamma*gamma)*(u*u-u*math.exp(-r0*(1/u-1)))
        return 2*s/math.sqrt(d)
    v,e=quad(fn,0,1,epsabs=1e-10,epsrel=1e-10)
    return math.pi-2*gamma*v


def orbit(w,b,distance):
    def rhs(t,y):
        r=math.hypot(y[0],y[1]);f=math.exp(-r)*(1+r)/(2*w*w*r**3)
        return [y[2],y[3],f*y[0],f*y[1]]
    def outward(t,y):return math.hypot(y[0],y[1])-math.hypot(distance,b)
    outward.direction=1;outward.terminal=True
    result=solve_ivp(rhs,(0,4*distance+100),[-distance,b,1.,0.],method='DOP853',
                     rtol=2e-11,atol=2e-12,events=outward,max_step=1.)
    if not result.success or len(result.t_events[0])!=1:raise ValueError('orbit did not exit')
    y=result.y[:,-1]
    return math.atan2(y[3],y[2]),len(result.t)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);a=p.parse_args()
    if a.output.exists():raise FileExistsError('preserve independent orbit controls')
    records=[]
    for w,b in ((.1,.1),(.1,1.),(.1,3.),(1.,.1),(1.,.5),(1.,1.),(1.,3.),(10.,.01),(10.,.1),(10.,1.)):
        angle=angle_quadrature(w,b)
        x,count=orbit(w,b,32.);y,count2=orbit(w,b,40.)
        records.append(dict(w=w,impact=b,quadrature_angle=angle,orbit_angle=y,
                            relative_error=abs(y/angle-1),distance_relative_change=abs(x/y-1),
                            steps=count2))
    sections=[]
    for w in (.1,1.,10.):
        # Integrate the impact parameter itself, independent of the turning-
        # radius Jacobian used in the main cross-section implementation.
        values=[]
        for n in (1,2):
            def fn(logb):
                b=math.exp(logb);theta=angle_quadrature(w,b)
                transfer=2*math.sin(theta/2)**2 if n==1 else math.sin(theta)**2
                return b*b*transfer
            value,error=quad(fn,-18,math.log(25.),epsabs=1e-10,epsrel=1e-7,limit=200)
            values.append(value)
        turning=cross_sections(w,angle_order=256,radius_order=512)
        sections.append(dict(w=w,impact_integral=values,turning_integral=turning.tolist(),
                             relative_difference=(np.array(values)/turning-1).tolist(),
                             source_fit=[cross_section_fit(w,n) for n in (1,2)]))
    passed=max(r['relative_error'] for r in records)<1e-7 and max(abs(v) for r in sections for v in r['relative_difference'])<1e-6
    paths=[Path(__file__),Path('scripts/yukawa_scattering_control.py'),Path('scripts/ion_collision_integrals.py')]
    report=dict(scope=__doc__,status='pass' if passed else 'fail',orbits=records,cross_sections=sections,
                input_sha256={str(p.resolve()):digest(p) for p in paths})
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
