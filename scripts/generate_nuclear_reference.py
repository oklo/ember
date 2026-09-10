#!/usr/bin/env python3
"""Independent SFII integrals and finite-degeneracy screening reference.

Requires scipy for this optional regeneration only; CTest reads the output.
Uses energy-space adaptive integration and bracketed electron inversion,
independent of ember's log-energy Gauss quadrature and Newton solver.
Sources and scope: docs/NUCLEAR.md. Numerical outputs, not source tables.
"""
import argparse
import math
from pathlib import Path
from scipy.integrate import quad
from scipy.optimize import brentq

KB, NA, H, ME, MU = 1.380649e-16, 6.02214076e23, 6.62607015e-27, 9.1093837015e-28, 1.66053906660e-24
C, E2, MEV = 2.99792458e10, (4.803204673e-10)**2, 1.602176634e-6
MASSES = [1.00782503*MU-ME, 3.01602932*MU-2*ME, 4.00260325*MU-2*ME]
REACTIONS = [(1,1,0,0,4.01e-25,4.01e-25*11.2,0), (2,2,1,1,5.21,-4.9,22), (2,2,1,2,.00056,-.00036,.000151)]

def rate(temp, reaction):
    z1,z2,i,j,s0,s1,s2 = REACTIONS[reaction]
    mu = MASSES[i]*MASSES[j]/(MASSES[i]+MASSES[j])
    kt = KB*temp
    eg = 2*mu*(2*math.pi**2*z1*z2*E2/H)**2
    peak = (eg/(4*kt))**(1/3)  # E0/kT
    # Scale away S0 and the Gamow exponential for accurate relative integrals.
    def f(x):
        if x == 0:
            return 0.
        energy = kt*x/MEV
        return (1+(s1/s0)*energy+.5*(s2/s0)*energy**2)*math.exp(3*peak-x-math.sqrt(eg/(kt*x)))
    top = max(100.,5*peak)
    opts = dict(epsabs=1e-12,epsrel=3e-13,limit=300,points=[peak])
    integral = quad(f,0,top,**opts)[0]
    moment = quad(lambda x:x*f(x),0,top,**opts)[0]
    value = NA*math.sqrt(8/(math.pi*mu*kt))*s0*MEV*1e-24*math.exp(-3*peak)*integral
    return value,-1.5+moment/integral

def screen(temp,rho,X,X3):
    ye=X+2*X3/3+(1-X-X3)/2
    ions=X+4*X3/3+(1-X-X3)
    ne=rho*NA*ye
    beta=KB*temp/(ME*C*C)
    norm=8*math.pi*(2*ME*KB*temp)**1.5/H**3
    def moments(eta):
        top=math.sqrt(max(eta,0)+60)
        def f(t, derivative=False):
            u=t*t-eta
            z=math.exp(-abs(u))
            occupation=z/(1+z) if u>0 else 1/(1+z)
            kernel=t*t*math.sqrt(1+.5*beta*t*t)*(1+beta*t*t)
            return kernel*(z/(1+z)**2 if derivative else occupation)
        points=[math.sqrt(eta)] if eta>0 else None
        return tuple(quad(lambda t:f(t,d),0,top,epsabs=1e-13,epsrel=3e-12,points=points,limit=300)[0] for d in [False,True])
    eta=brentq(lambda eta:math.log(norm*moments(eta)[0]/ne),-100,1e4,xtol=1e-12)
    n,derivative=moments(eta)
    theta=derivative/n
    ge=E2/(KB*temp)*(4*math.pi*ne/3)**(1/3)
    weak=E2/(KB*temp)*math.sqrt(4*math.pi*E2*rho*NA*(ions+theta*ye)/(KB*temp))
    strong=.9*ge*(2**(5/3)-2)
    exponent=weak*strong/math.hypot(weak,strong)
    return eta,theta,ge,weak,exponent

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    rows=['# Independent scipy adaptive quadrature; see scripts/generate_nuclear_reference.py',
          '# R T reaction NA_sigma_v dlnrate_dlnT; S T rho X X3 eta theta Gamma_e H_DH H_SVH (Z=0 baryonic)']
    for temp in [1e5,1e6,4.5e6,8e6,1.55e7,2e7]:
        for r in range(3):
            rows.append('R '+' '.join(format(v,'.17g') for v in [temp,r,*rate(temp,r)]))
    for temp,rho in [(1e7,1e-4),(4.5e6,360),(1e6,100),(1.55e7,150),(1e6,1000)]:
        for X,X3 in [(.7,0),(.35,.005)]:
            rows.append('S '+' '.join(format(v,'.17g') for v in [temp,rho,X,X3,*screen(temp,rho,X,X3)]))
    args.output.write_text('\n'.join(rows)+'\n')

if __name__=='__main__':
    main()
