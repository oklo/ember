#!/usr/bin/env python3
"""Test the optical-conductivity weight implied by returned TOPS absorption.

Interpret the index-one absorption as sigma_1=c*rho*kappa_abs/(4*pi), then
compare its measured-frequency integral to pi*n_e*e^2/(2*m_e). This tests
whether those spectra alone supply the free-electron conductivity weight.
It does not recover omitted frequencies, isolate free-free from bound
absorption, infer a unique collision rate, or adopt an opacity correction.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.special import exprel

from audit_tops_spectral_means import source, digest

C = 2.99792458e10
HBAR = 1.054571817e-27
KEV = 1.602176634e-9
KB = 1.380649e-16


def integrate(energy, opacity):
    loge = np.log(energy)
    widths = np.diff(loge)
    slope = np.diff(np.log(opacity))/widths
    exact = float(np.sum(opacity[:-1]*energy[:-1]*widths*exprel((slope+1)*widths)))
    nodes, weights = leggauss(16)
    offset = widths[:, None]*(nodes+1)/2
    numerical = float(np.sum(opacity[:-1, None]*energy[:-1, None]
                     * np.exp((slope+1)[:, None]*offset)*widths[:, None]/2*weights))
    if not exact > 0 or abs(numerical/exact-1) > 1e-10:
        raise ValueError('spectral area quadrature disagrees')
    return exact, abs(numerical/exact-1), slope


def controls():
    errors = []
    for power in [-2, -1, 0, 2]:
        e = np.geomspace(.01, 100, 97)
        value, _, _ = integrate(e, e**power)
        expected = math.log(e[-1]/e[0]) if power == -1 else (e[-1]**(power+1)-e[0]**(power+1))/(power+1)
        errors.append(abs(value/expected-1))
    # Unit plasma frequency/collision rate/rho/c: Drude absorption=1/(1+w^2).
    drude = 2/math.pi*quad(lambda x: 1/(1+x*x), 0, np.inf, epsabs=1e-12)[0]
    if max(errors) > 1e-12 or abs(drude-1) > 1e-12:
        raise ValueError('analytic spectral-weight controls failed')
    return {'power_law_max_relative_error': max(errors),
            'drude_sum_rule_relative_error': abs(drude-1)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['work', 'classical_check', 'output']:
        p.add_argument(name, type=Path)
    a = p.parse_args()
    check = json.loads(a.classical_check.read_text())
    inputs = {str(v.resolve()): digest(v) for v in [Path(__file__), a.classical_check,
              Path(__file__).with_name('audit_tops_spectral_means.py')]}
    for path, checksum in check['input_sha256'].items():
        if digest(Path(path)) != checksum:raise ValueError('source input changed')
        inputs[path] = checksum
    spectra, _ = source(a.work/'cutoff-on')
    rows = []
    for row in check['records']:
        key = next(k for k in spectra if k[1] == row['density'] and
                   abs(k[0]*KEV/KB/row['temperature_K']-1) < 1e-12)
        data = spectra[key]
        energy, opacity = data[:, 0], data[:, 2]
        area, error, slope = integrate(energy, opacity)
        wp = row['classical_cutoff_u']*key[0]*KEV/HBAR
        factor = 2*C*row['density']/math.pi*KEV/HBAR/(wp*wp)
        tails = {}
        if slope[0] > -1:
            tails['low_frequency_power_law_extrapolation'] = factor*energy[0]*opacity[0]/(slope[0]+1)
        if slope[-1] < -1:
            tails['high_frequency_power_law_extrapolation'] = -factor*energy[-1]*opacity[-1]/(slope[-1]+1)
        rows.append({'temperature_keV': key[0], 'density': key[1],
                     'free_electron_density': row['estimated_electron_density'],
                     'energy_range_keV': [float(energy[0]), float(energy[-1])],
                     'integrated_absorption_keV_cm2_per_g': area,
                     'measured_weight_over_classical_free_electron_sum': factor*area,
                     'low_frequency_logarithmic_slope': float(slope[0]),
                     'high_frequency_logarithmic_slope': float(slope[-1]),
                     'formal_tail_weight_ratios': tails,
                     'tail_note': 'Extrapolations of endpoint slopes only, not physical tail models or bounds.',
                     'quadrature_relative_error': error})
    if any(digest(Path(path)) != checksum for path, checksum in inputs.items()):
        raise ValueError('source changed during calculation')
    report = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
              'references': ['https://doi.org/10.1016/j.hedp.2017.02.008',
                             'https://doi.org/10.1103/PhysRevE.91.053102'],
              'input_sha256': inputs, 'analytic_checks': controls(), 'records': rows,
              'interpretation': 'Do not use the returned absorption alone to reconstruct a complete causal dielectric function. Its finite-frequency conductivity weight must be checked with a supported treatment of missing response and existing multiple-collision corrections.'}
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps([r['measured_weight_over_classical_free_electron_sum'] for r in rows]))


if __name__ == '__main__':
    main()
