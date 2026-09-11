#!/usr/bin/env python3
"""Estimate central C/N proton-capture times without changing stellar evolution.

Solar Fusion II Table XII supplies the low-energy S-factor expansions. These
instantaneous times do not integrate the changing stellar state, establish a
global CNO luminosity bound, or evolve the metal mixture. SciPy is required.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

from scipy.integrate import quad
from numpy.polynomial.legendre import leggauss

KB, NA, H = 1.380649e-16, 6.02214076e23, 6.62607015e-27
MU, ME = 1.66053906660e-24, 9.1093837015e-28
E2, KEV, YR = (4.803204673e-10)**2, 1.602176634e-9, 3.1557e7
# S0 [keV barn], S1 [barn], S2 [barn/keV]. Atomic masses match
# Ember's existing nuclide constants; electron binding energies are omitted.
REACTIONS = {
    'C12_p_gamma': (6, 12., 1.34, .0026, .000083),
    'C13_p_gamma': (6, 13.00335484, 7.6, -.00783, .000729),
    'N14_p_gamma': (7, 14.003074, 1.66, -.0033, .000044),
}
GAUSS_X, GAUSS_W = leggauss(64)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rate(T, parameters):
    if not 1e6 <= T <= 1e7:
        raise ValueError('central diagnostic restricted to 1--10 million K')
    z, mass, s0, s1, s2 = parameters
    proton, target = 1.00782503 * MU - ME, mass * MU - z * ME
    mu = proton * target / (proton + target)
    kt = KB * T
    eg = 2 * mu * (2 * math.pi**2 * z * E2 / H)**2
    peak = (eg / (4 * kt))**(1/3)
    top = max(100., 5 * peak)
    ceiling = 130 * KEV / kt
    if top >= ceiling:
        raise ValueError('quadrature reaches the 130 keV diagnostic limit')

    def kernel(x):
        if x == 0:
            return 0.
        energy = kt * x / KEV
        return (1 + s1/s0 * energy + .5*s2/s0 * energy**2) * math.exp(
            3*peak - x - math.sqrt(eg/(kt*x)))

    integral, error = quad(kernel, 0, top, points=[peak],
                           epsabs=1e-12, epsrel=3e-13, limit=300)
    tail = quad(kernel, top, ceiling, epsabs=1e-14, epsrel=1e-10)[0]
    logarithmic_integral = 0.
    for lo, hi in [(math.log(peak)-8, math.log(peak)),
                   (math.log(peak), math.log(top))]:
        for node, weight in zip(GAUSS_X, GAUSS_W):
            x = math.exp(.5 * ((hi-lo)*node + hi+lo))
            logarithmic_integral += .5*(hi-lo)*weight*x*kernel(x)
    agreement = abs(logarithmic_integral/integral-1)
    if error/integral > 1e-10 or tail/integral > 1e-10 or agreement > 1e-10:
        raise ValueError('quadrature precision or upper integration range insufficient')
    value = NA * math.sqrt(8/(math.pi*mu*kt)) * s0 * KEV * 1e-24 * math.exp(-3*peak) * integral
    return value, eg, {'Gamow_peak_keV': peak*kt/KEV,
                       'upper_energy_keV': top*kt/KEV,
                       'quadrature_relative_error_estimate': error/integral,
                       'log_energy_gauss_relative_difference': agreement,
                       'polynomial_tail_to_130keV_fraction': tail/integral}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('probe', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('tracks', nargs='+', type=Path)
    a = p.parse_args()
    records, inputs = [], {str(a.probe): sha(a.probe), __file__: sha(__file__)}
    for path in a.tracks:
        inputs[str(path)] = sha(path)
        track = json.loads(path.read_text())
        row = dict(zip(track['columns'], track['history'][-1], strict=True))
        T, rho, X, X3 = (row[k] for k in ['central_T_K', 'central_rho', 'central_X', 'central_Y3'])
        if not X > 0 or track['mass_Msun'] != .1:
            raise ValueError('requires hydrogen-bearing 0.1 solar-mass central states')
        query = ' '.join(format(v, '.17g') for v in [T, rho, X, X3]) + '\n'
        result = subprocess.run([str(a.probe), '--gs98'], input=query, text=True,
                                capture_output=True, check=True)
        values = list(map(float, result.stdout.split()))
        if len(values) != 20 or values[:2] != [T, rho] or not all(map(math.isfinite, values)):
            raise ValueError('unexpected nuclear probe result')
        weak_pp, gamma = values[8], values[16]
        reactions = {}
        for name, parameters in REACTIONS.items():
            z = parameters[0]
            bare, eg, quadrature = rate(T, parameters)
            weak = z * weak_pp
            strong = .9 * gamma * ((z+1)**(5/3) - z**(5/3) - 1)
            exponent = weak * strong / math.hypot(weak, strong)
            zeta = gamma * 2*z/(1+z**(1/3)) / (eg/(4*KB*T))**(1/3)
            if zeta > .2:
                raise ValueError('classical-ion screening domain exceeded')
            bare_time = 1/(rho * X * bare) / YR
            reactions[name] = {'bare_molar_rate_cm3_mol_s': bare,
                               'bare_capture_time_yr': bare_time,
                               'SVH_log_screening_factor': exponent,
                               'SVH_capture_time_yr': bare_time * math.exp(-exponent),
                               'screening_zeta': zeta, 'quadrature': quadrature}
        records.append({'track': str(path), 'age_yr': row['age_yr'], 'T_K': T,
                        'rho_g_cm3': rho, 'XH': X, 'X3': X3,
                        'pp_network_heat_erg_g_s': values[13], 'reactions': reactions})
    report = {'scope': __doc__, 'source': 'https://arxiv.org/pdf/1004.2318',
              'source_locations': 'Table XII; N14 energy expansion and range in equation 50.',
              'rate_parameters': REACTIONS,
              'screening_scope': 'Same finite-degeneracy susceptibility and Salpeter--Van Horn interpolation as Ember, applied to proton/carbon and proton/nitrogen charges. This is a screening approximation, not an uncertainty bound.',
              'range_note': 'All integrated energies are below 130 keV. The omitted polynomial tail to 130 keV is checked; this is not a bound on unmodeled higher-energy cross sections.',
              'interpretation': 'A capture time comparable to the time spent at a stellar state calls for a changing-composition calculation. Dividing by the total stellar age does not determine depletion. No CNO energy or abundance change is installed.',
              'records': records, 'input_sha256': inputs}
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    for row in records:
        print(f"Age {row['age_yr']/1e12:.4g} trillion yr:",
              {k: format(v['SVH_capture_time_yr']/1e9, '.4g') + ' billion yr'
               for k, v in row['reactions'].items()})


if __name__ == '__main__':
    main()
