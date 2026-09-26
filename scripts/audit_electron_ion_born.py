#!/usr/bin/env python3
"""Independent controls for a displaced-Fermi, screened Born drag moment.

Checks concern this explicitly limited kinetic model, not its suitability
for stellar diffusion. No heat flow, ionic correlations or abundance change.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.special import exp1, expit

from electron_ion_born import (KB, ME, HBAR, E2, cross_section,
                              coulomb_bracket, fermi_half, chemical_potential,
                              drift_integral, pair_resistance)
from fetch_tops_composition import digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('preserve collision comparisons')
    inputs={str(p.resolve()):digest(p) for p in
            (Path(__file__),Path('scripts/electron_ion_born.py'))}
    checks=[]
    def check(name,got,reference,tolerance,**context):
        error=abs(got/reference-1)
        checks.append(dict(name=name,value=float(got),reference=float(reference),
                           relative_error=float(error),tolerance=tolerance,
                           passed=bool(np.isfinite(error) and error<tolerance),**context))

    # Integrate the differential Born cross section over scattering angle.
    length=3e-9;charge=2.
    for b in np.logspace(-10,4,18):
        momentum=math.sqrt(b)*HBAR/(2*length)
        prefactor=4*ME**2*(charge*E2)**2*length**4/HBAR**4
        angular=8*math.pi*prefactor*quad(lambda u:u/(1+b*u)**2,0,1,
                    epsabs=1e-100,epsrel=2e-12,limit=250)[0]
        check('angular_cross_section',cross_section(momentum,charge,length),angular,1e-9,b=float(b))
    check('zero_momentum_cross_section',cross_section(0.,charge,length),
          16*math.pi*ME**2*(charge*E2)**2*length**4/HBAR**4,1e-13)

    # Maxwell-Boltzmann and Sommerfeld number-integral limits.
    check('classical_number_limit',fermi_half(-25),math.sqrt(math.pi)/2*math.exp(-25),1e-10)
    eta=1e4
    sommerfeld=2/3*eta**1.5+math.pi**2/(12*math.sqrt(eta))+7*math.pi**4/(960*eta**2.5)
    check('degenerate_number_limit',fermi_half(eta),sommerfeld,1e-12)
    # Independent scalar energy integral, including its lower endpoint.
    for eta in (-12.,-2.,0.,4.,18.,35.):
        direct=quad(lambda x:math.sqrt(x)*expit(eta-x),0,max(eta,0)+60,
                    epsabs=1e-100,epsrel=2e-12,limit=250)[0]
        check('number_energy_integral',fermi_half(eta),direct,1e-10,eta=eta)

    # Exact Maxwellian energy integral for the same Born potential.
    for b in (.05,.3,1.,10.,100.):
        analytic=(1+1/b)*math.exp(1/b)*exp1(1/b)-1
        check('classical_drag_limit',drift_integral(-25.,b)/math.exp(-25),analytic,2e-9,b_thermal=b)

    # Integrate a genuinely displaced distribution in momentum and angle.
    # This does not use the analytic f(1-f) linearization of the drag module.
    mu,weights=leggauss(64);mu=(mu+1)/2;weights=weights/2
    T=1e7;pth=math.sqrt(2*ME*KB*T);ni=1e25
    for eta,b in ((-8.,.1),(-8.,10.),(0.,.1),(0.,10.),(6.,.1),(6.,10.)):
        length=math.sqrt(b*HBAR**2/(8*ME*KB*T))
        ne=pth**3/(2*math.pi**2*HBAR**3)*fermi_half(eta)
        predicted=pair_resistance(electron_density=ne,ion_density=ni,temperature=T,
                    charge=2,screening_length=length,eta=eta)['resistance_g_cm3_s']
        def displaced(delta):
            def radial(t):
                energy=t*t+delta*delta-eta
                angular=np.sum(weights*mu*(expit(-energy+2*t*delta*mu)-expit(-energy-2*t*delta*mu)))
                return t**4*float(cross_section(pth*t,2,length))*angular
            value=quad(radial,0,math.sqrt(max(eta,0)+50),epsabs=1e-100,
                       epsrel=2e-10,limit=200)[0]
            return ni*pth**4/(2*math.pi**2*HBAR**3*delta)*value
        coarse,fine=displaced(1e-3),displaced(5e-4)
        check('displaced_distribution',predicted,(4*fine-coarse)/3,2e-8,eta=eta,b_thermal=b)
        check('drift_amplitude_convergence',fine,coarse,2e-6,eta=eta,b_thermal=b)

    ne=1e27;ni=ne/2;length=1e-9;pf=HBAR*(3*math.pi**2*ne)**(1/3)
    # Independent degenerate collision frequency: PBHY99 equations7--8,
    # with S(q)=1, nonrelativistic vF and Yukawa u(q).
    kf=pf/HBAR;ks=1/length
    coulomb=quad(lambda q:q**3/(q*q+ks*ks)**2,0,2*kf,
                 epsabs=1e-100,epsrel=2e-12)[0]
    nu=4*math.pi*4*E2**2*ni/(pf*pf*(pf/ME))*coulomb
    cold_reference=ME*ne*nu
    cold=[]
    for theta in (1e-2,1e-3,1e-4):
        T=theta*pf**2/(2*ME*KB)
        state=pair_resistance(electron_density=ne,ion_density=ni,temperature=T,
                             charge=2,screening_length=length)
        cold.append(dict(theta=theta,resistance=state['resistance_g_cm3_s'],
                         relative_to_zero_temperature=state['resistance_g_cm3_s']/cold_reference-1))
    check('degenerate_PBHY_limit',cold[-1]['resistance'],cold_reference,2e-7)
    check('degenerate_error_scales_with_T_squared',
          abs(cold[1]['relative_to_zero_temperature']/cold[2]['relative_to_zero_temperature']),100.,1e-3)

    # SI calculation of the same moment, with independent unit conversion.
    T=1e7;eta=chemical_potential(ne,T);moment=drift_integral(eta,8*ME*KB*T*length**2/HBAR**2)
    SI=(ni*1e6)*2*(ME*1e-3)**2*(2*E2*1e-9)**2/(3*math.pi*(HBAR*1e-7)**3)*moment
    cgs=pair_resistance(electron_density=ne,ion_density=ni,temperature=T,
                        charge=2,screening_length=length,eta=eta)['resistance_g_cm3_s']
    check('SI_cgs',SI,cgs*1e3,1e-12)
    check('longer_thermal_tail',drift_integral(eta,1.,tail=60.),drift_integral(eta,1.),1e-10)

    invalid=[]
    for label,call in [('negative_density',lambda:chemical_potential(-1,1e7)),
                       ('zero_temperature',lambda:chemical_potential(1e25,0)),
                       ('negative_momentum',lambda:cross_section(-1,1,1e-9)),
                       ('neutral_ion',lambda:cross_section(1e-18,0,1e-9)),
                       ('zero_screening',lambda:cross_section(1e-18,1,0)),
                       ('inconsistent_eta',lambda:pair_resistance(electron_density=ne,ion_density=ni,
                            temperature=T,charge=2,screening_length=length,eta=-25))]:
        try:call();raised=False
        except ValueError:raised=True
        invalid.append(dict(name=label,rejected=raised))
    for path,h in inputs.items():
        if digest(Path(path))!=h:raise ValueError('collision source changed during checks')
    passed=all(c['passed'] for c in checks) and all(c['rejected'] for c in invalid)
    result=dict(scope=__doc__,status='pass' if passed else 'fail',accepted_for_stellar_evolution=False,
                checks=checks,invalid_input_checks=invalid,degenerate_convergence=cold,
                maximum_relative_error=max(c['relative_error'] for c in checks),input_sha256=inputs,
                limitations=['First Born, nonrelativistic, stationary point ions, specified Yukawa screening.',
                             'No ion structure factor, electron heat perturbation or full collision-operator solution.',
                             'No stellar diffusion velocities or abundance evolution.'])
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=result['status'],checks=len(checks),invalid_inputs=len(invalid),
                         failures=[c for c in checks if not c['passed']])) ,flush=True)
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
