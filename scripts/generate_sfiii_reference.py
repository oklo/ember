#!/usr/bin/env python3
"""Independent Solar Fusion III rate integrals for the reduced pp network.

Acharya et al. (2025), https://arxiv.org/abs/2405.06470v3:
equations 8--9, section V.C, and equation 14. Integrate in linear energy
with adaptive quadrature; Ember uses fixed panels in logarithmic energy.
The 34 fit is never evaluated beyond its stated 1.6 MeV upper limit.
"""
import argparse
import math
from pathlib import Path
from scipy.integrate import quad

from generate_nuclear_reference import KB, NA, H, E2, MEV, MASSES


def rate(temperature, reaction):
    if not 1e5 <= temperature <= 2e7:
        raise ValueError('reference temperature outside the selected range')
    if reaction == 0:
        z1, z2, i, j, s0 = 1, 1, 0, 0, 4.09e-25
        def shape(energy):
            return 1 + 11.0*energy + 121.0*energy**2
    elif reaction == 1:
        z1, z2, i, j, s0 = 2, 2, 1, 1, 5.21
        def shape(energy):
            return (5.21 - 4.90*energy + 11.21*energy**2)/s0
    elif reaction == 2:
        z1, z2, i, j, s0 = 2, 2, 1, 2, .0005610
        def shape(energy):
            if not 0 <= energy <= 1.6:
                raise ValueError('34 fit outside its stated energy range')
            return math.exp(-.5374*energy)*(1-.4829*energy**2
                                         +.6310*energy**3-.1527*energy**4)
    else:
        raise ValueError('unknown reaction')
    mu = MASSES[i]*MASSES[j]/(MASSES[i]+MASSES[j])
    kt = KB*temperature
    eg = 2*mu*(2*math.pi**2*z1*z2*E2/H)**2
    peak = (eg/(4*kt))**(1/3)
    def integrand(x):
        if x == 0:
            return 0.
        return shape(kt*x/MEV)*math.exp(3*peak-x-math.sqrt(eg/(kt*x)))
    top = max(100., 5*peak)
    if reaction == 2 and 2*top*kt/MEV > 1.6:
        raise ValueError('tail check would leave the 34 fit domain')
    options = dict(epsabs=1e-12, epsrel=3e-13, limit=300, points=[peak])
    integral = quad(integrand, 0, top, **options)[0]
    moment = quad(lambda x: x*integrand(x), 0, top, **options)[0]
    tail = quad(integrand, top, 2*top, epsabs=1e-15, epsrel=1e-10)[0]
    if not integral > 0 or abs(tail/integral) > 1e-12:
        raise ValueError('tail criterion failed')
    value = NA*math.sqrt(8/(math.pi*mu*kt))*s0*MEV*1e-24*math.exp(-3*peak)*integral
    return value, -1.5+moment/integral


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    rows = ['# Solar Fusion III; independent adaptive energy integration.',
            '# Source: https://arxiv.org/abs/2405.06470v3, equations 8, 9, 14 and section V.C.',
            '# T_K reaction NA_sigma_v dlnrate_dlnT']
    for temperature in [1e5, 1e6, 4.5e6, 8e6, 1.55e7, 2e7]:
        for reaction in range(3):
            rows.append(' '.join(format(v, '.17g') for v in
                                 [temperature, reaction, *rate(temperature, reaction)]))
    args.output.write_text('\n'.join(rows)+'\n')


if __name__ == '__main__':
    main()
