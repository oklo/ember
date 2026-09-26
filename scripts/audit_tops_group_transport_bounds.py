#!/usr/bin/env python3
"""Bound finite-band refractive transport using group means, without a scattering estimate.

For 0<=n<=1 and a nonnegative scattering fraction, n^3 <=
n^3/(1-(1-n)*f_scattering) <= n^2. Group endpoint indices therefore bound
the inverse-opacity integral even when its frequency dependence within a
group is unknown. These are bounds within the specified collisionless model
and supplied photon interval, not bounds on all missing dielectric physics.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad

from audit_tops_electron_dispersion import electron_moments, index_squared, KEV, KB, HBAR, ME, RW
from audit_tops_spectral_means import digest
from reduce_tops_group_factors import add_inputs, verify
from tops_groups import source_groups


def weight(u):
    return u**4 * np.exp(-u) / (-np.expm1(-u))**2


def group_weights(boundaries, order):
    nodes, weights = leggauss(order)
    lo, hi = boundaries[:-1], boundaries[1:]
    u = (lo[:, None] + hi[:, None]) / 2 + (hi - lo)[:, None] * nodes / 2
    return np.sum((hi - lo)[:, None] * weights * weight(u) / 2, axis=1)


def inverse_bounds(boundaries, rosseland, cutoff, vstar2, order=32):
    if boundaries[0] > cutoff:
        raise ValueError('photon interval starts above cutoff')
    indices2 = np.zeros_like(boundaries)
    valid = boundaries > cutoff
    indices2[valid] = index_squared(boundaries[valid], cutoff, vstar2)
    if np.any(np.diff(indices2) < -1e-14):
        raise ValueError('index is not monotone over this photon interval')
    inverse = group_weights(boundaries, order) / rosseland
    lower = float(np.sum(inverse * indices2[:-1]**1.5) / RW)
    upper = float(np.sum(inverse * indices2[1:]) / RW)
    if not 0 < lower <= upper or not math.isfinite(upper):
        raise ValueError('invalid finite-band transport bounds')
    return lower, upper


def independent_controls():
    boundaries = np.geomspace(.01, 200., 65)
    weights = group_weights(boundaries, 64)
    records = []
    for model in ('constant', 'varying'):
        opacity = lambda u: 2. if model == 'constant' else 2. + .1 * u**2
        integrals = [quad(lambda u: float(weight(u)) / opacity(u), lo, hi,
                          epsabs=1e-280, epsrel=1e-10)[0]
                     for lo, hi in zip(boundaries[:-1], boundaries[1:], strict=True)]
        groups = weights / np.array(integrals)
        for cutoff in (.1, 1., 5., 50.):
            lower, upper = inverse_bounds(boundaries, groups, cutoff, 0.)
            for fraction in (0., .5, 1.):
                def integrand(u):
                    n = math.sqrt(1 - (cutoff / u)**2)
                    return float(weight(u)) * n**3 / (opacity(u) * (1 - (1 - n) * fraction))
                actual = quad(integrand, cutoff, boundaries[-1], epsabs=1e-280,
                              epsrel=1e-9, limit=300)[0] / RW
                if not lower <= actual <= upper:
                    raise ValueError('independent analytic-index spectrum lies outside bounds')
                records.append({'opacity_model': model, 'cutoff_u': cutoff,
                                'scattering_fraction': fraction,
                                'inverse_mean_lower': lower, 'inverse_mean': actual,
                                'inverse_mean_upper': upper})
    return records


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('reduction', 'conduction_reference', 'conduction_probe', 'mass_probe', 'scratch', 'output'):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    if a.scratch.exists() or a.output.exists():
        raise FileExistsError('use new diagnostic outputs')
    old = json.loads(a.reduction.read_text())
    conduction_ref = json.loads(a.conduction_reference.read_text())
    inputs = dict(old['input_sha256'])
    add_inputs(inputs, conduction_ref['input_sha256'])
    paths = [a.reduction, a.conduction_reference, a.mass_probe, Path(__file__),
             Path('/tmp/ember-opacity-mass-basis-v1.cpp'),
             Path('scripts/reduce_tops_group_factors.py')]
    for path in paths:
        add_inputs(inputs, {str(path.resolve()): digest(path)})
    if inputs.get(str(a.conduction_probe.resolve())) != digest(a.conduction_probe):
        raise ValueError('conduction executable is not identified')
    verify(inputs)
    controls = independent_controls()
    selected = [r for r in old['records'] if r['density_atomic_g_cm3'] >= 10000]
    if not selected:
        raise ValueError('no dense source states')
    a.scratch.mkdir(parents=True)
    hscale, hescale, zscale = map(float, subprocess.check_output([str(a.mass_probe)], text=True).split())
    x, z = old['X'], old['Z']
    scale = 1 / (x / hscale + (1 - x - z) / hescale + z / zscale)
    xb, zb = scale * x / hscale, scale * z / zscale
    if max(abs(xb * hscale / scale - x), abs(zb * zscale / scale - z)) > 1e-14:
        raise ValueError('mass-basis round trip failed')
    lines = [f"1 1 {r['density_atomic_g_cm3'] / scale:.17g} {r['temperature_keV'] * KEV / KB:.17g} 1 {xb:.17g} 0 {1-xb-zb:.17g}\n"
             for r in selected]
    request = a.scratch / 'conduction.txt'
    request.write_text(''.join(lines))
    response = subprocess.run([str(a.conduction_probe), 'data/conduction/condtab21wd_metals.dat',
                               str(request), format(zb, '.17g')], text=True, capture_output=True, timeout=60)
    (a.scratch / 'conduction.stdout').write_text(response.stdout)
    (a.scratch / 'conduction.stderr').write_text(response.stderr)
    response.check_returncode()
    lines = response.stdout.splitlines()
    if len(lines) != len(selected) + 1 or not lines[-1].startswith('checksum '):
        raise ValueError('incomplete conduction response')
    k_conduction = [float(line.split()[0]) / scale for line in lines[:-1]]
    if any(not math.isfinite(k) or k <= 0 for k in k_conduction):
        raise ValueError('invalid conduction mean')
    sources = {root: source_groups(root) for root in sorted({r['source'] for r in selected})}
    records = []
    for row, kc in zip(selected, k_conduction, strict=True):
        t, rho = row['temperature_keV'], row['density_atomic_g_cm3']
        boundaries, groups, grey = sources[row['source']]
        ne = grey[t, rho]['electron_density_cm3']
        e = electron_moments(t * KEV / KB, ne)
        cutoff = HBAR * math.sqrt(4 * math.pi * 4.80320471257e-10**2 * ne / ME) / (t * KEV)
        cutoff *= math.sqrt(e['plasma_frequency_squared_ratio'])
        pair = [inverse_bounds(boundaries / t, groups[t, rho][:, 1], cutoff, e['vstar_squared'], n)
                for n in (32, 64)]
        lower, upper = pair[1]
        refinement = max(abs(x / y - 1) for x, y in zip(pair[0], pair[1], strict=True))
        if refinement > 1e-7:
            raise ValueError('bound weight integration did not converge')
        estimate = 1 / row['rosseland_atomic_cm2_g']
        if not lower <= estimate <= upper:
            raise ValueError('group transport estimate lies outside its conditional bounds')
        normalization = row['uncut_source_rosseland'] / row['uncut_recombined_rosseland']
        lower_normalized, upper_normalized = lower / normalization, upper / normalization
        records.append({'temperature_keV': t, 'density_atomic_g_cm3': rho,
                        'source': row['source'], 'cutoff_u': cutoff,
                        'photon_max_u': float(boundaries[-1] / t),
                        'inverse_opacity_lower': lower, 'inverse_opacity_upper': upper,
                        'inverse_opacity_estimate': estimate,
                        'bound_quadrature_relative_change': refinement,
                        'native_uncut_normalization': normalization,
                        'conductive_opacity_atomic_cm2_g': kc,
                        'normalized_finite_band_radiative_fraction_upper': upper_normalized / (1 / kc + upper_normalized),
                        'normalized_finite_band_radiative_fraction_lower': lower_normalized / (1 / kc + lower_normalized),
                        'normalized_combined_opacity_relative_span': (1 / kc + upper_normalized) / (1 / kc + lower_normalized) - 1})
    verify(inputs)
    result = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
              'source_X_atomic': x, 'source_Z_atomic': z,
              'mass_scale_atomic_over_baryonic': scale,
              'independent_synthetic_controls': controls,
              'dense_source_states': len(records),
              'maximum_bound_quadrature_relative_change': max(r['bound_quadrature_relative_change'] for r in records),
              'finite_band_radiative_fraction_below_threshold_counts': {
                  str(threshold): sum(r['normalized_finite_band_radiative_fraction_upper'] < threshold for r in records)
                  for threshold in (1e-2, 1e-3, 1e-4)},
              'limitations': ['Bounds apply within the collisionless free-electron prescription and supplied frequency interval.',
                              'No independent upper bound on the omitted high-frequency inverse-opacity integral is included.',
                              'Native uncut normalization is a construction, not a separate accuracy proof.',
                              'Conductive transport has the physical assumptions of the selected conduction table.',
                              'These bounds alone do not authorize omitting radiation or accepting a complete opacity family.'],
              'input_sha256': inputs, 'records': records,
              'scratch_sha256': {p.name: digest(p) for p in a.scratch.iterdir() if p.is_file()}}
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('scope', 'records', 'input_sha256',
                     'scratch_sha256', 'independent_synthetic_controls', 'limitations')}), flush=True)


if __name__ == '__main__':
    main()
