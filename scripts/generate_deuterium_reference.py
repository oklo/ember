#!/usr/bin/env python3
"""Independent energy-space integration of Solar Fusion III's S12 table.

The native calculation integrates log energy using fixed Gauss--Legendre
quadrature. This reference uses adaptive integration in E/kT and also records
the sensitivity to monotone cubic rather than linear S-factor interpolation.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
ENERGY = np.array([0, .010, .020, .040, .080, .091, .100, .120])
S = np.array([2.028, 2.644, 3.276, 4.579, 7.31, 8.11, 8.77, 10.24])
KB, AMU, ME, H, NA, MEV = 1.380649e-16, 1.66053906660e-24, 9.1093837015e-28, 6.62607015e-27, 6.02214076e23, 1.602176634e-6
mp, md = 1.00782503*AMU-ME, 2.01410177812*AMU-ME
MU = mp*md/(mp+md)
EG = 2*MU*(np.pi*(4.803204673e-10)**2/(H/(2*np.pi)))**2


def integrate(temperature, cubic=False):
    kt = KB*temperature
    pchip = PchipInterpolator(ENERGY, S, extrapolate=False)
    def integrand(x, moment=0):
        if x == 0:
            return 0.
        energy = x*kt/MEV
        sf = float(pchip(energy)) if cubic else np.interp(energy, ENERGY, S)
        return sf*np.exp(-x-np.sqrt(EG/(x*kt)))*x**moment
    peak = (EG/(4*kt))**(1/3)
    edges = sorted(set([0., peak/2, peak, 2*peak] + list(ENERGY[1:]*MEV/kt)))
    values = [sum(quad(lambda x: integrand(x, moment), a, b, epsabs=1e-250,
                       epsrel=2e-12, limit=200)[0] for a, b in zip(edges, edges[1:]))
              for moment in (0, 1)]
    rate = NA*np.sqrt(8/(np.pi*MU*kt))*(1e-7*MEV*1e-24)*values[0]
    return rate, -1.5+values[1]/values[0]


def main():
    rows = []
    for t in [1e4, 3e4, 1e5, 2e5, 3e5, 5e5, 7e5, 1e6, 1.5e6, 2e6, 3e6, 5e6, 1e7, 2e7]:
        rate, slope = integrate(t)
        cubic, _ = integrate(t, True)
        rows.append(dict(T_K=t, molar_rate=rate, dlnrate_dlnT=slope,
                         cubic_over_linear_minus_one=cubic/rate-1))
    data = dict(source='https://arxiv.org/html/2405.06470v3#S4.T4',
                source_table='Solar Fusion III, Table IV',
                method='scipy QUADPACK, E/kT, splits at Gamow peak and S-table joins',
                energy_MeV=ENERGY.tolist(), S_in_1e_minus_7_MeV_barn=S.tolist(), rows=rows)
    out = ROOT/'tests/data'
    (out/'deuterium_reference.json').write_text(json.dumps(data, indent=2)+'\n')
    (out/'deuterium_reference.txt').write_text(''.join(
        f"{r['T_K']:.17g} {r['molar_rate']:.17g} {r['dlnrate_dlnT']:.17g}\n" for r in rows))
    print(json.dumps(dict(cases=len(rows), maximum_interpolation_difference=max(
        abs(r['cubic_over_linear_minus_one']) for r in rows),
        maximum_PMS_interpolation_difference=max(abs(r['cubic_over_linear_minus_one'])
        for r in rows if r['T_K']<=3e6)), indent=2))


if __name__ == '__main__':
    main()
