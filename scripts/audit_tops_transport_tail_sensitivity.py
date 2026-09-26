#!/usr/bin/env python3
"""Quantify a conditional high-energy tail bound for saved transport groups.

This does not derive a dense-plasma opacity beyond the supplied photon grid.
It states a lower-opacity envelope explicitly, tests its normalization over
twelve orders of magnitude, and reports the resulting transport bounds. The
envelope is a sensitivity assumption, not a new measured source value.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import logsumexp

from audit_tops_spectral_means import digest
from audit_tops_electron_dispersion import RW
from reduce_tops_group_factors import add_inputs, verify
from tops_groups import source_groups


def log_tail_upper(u0, k0, power, scale):
    """Bound integral w(u)/k(u)/RW for k>=scale*k0*(u/u0)**(-power).

    Above u0, (1-exp(-u))**-2 <= (1-exp(-u0))**-2. The remaining
    exponential times integer power integrates to an incomplete gamma
    function, evaluated here as its finite sum in logarithmic form.
    """
    if not (u0 > 0 and k0 > 0 and scale > 0 and
            isinstance(power, int) and 0 <= power <= 8):
        raise ValueError('invalid tail-envelope parameters')
    if not all(math.isfinite(v) for v in (u0, k0, scale)):
        raise ValueError('nonfinite tail-envelope parameters')
    n = 4 + power
    log_gamma = (math.lgamma(n + 1) - u0 +
                 logsumexp([j * math.log(u0) - math.lgamma(j + 1)
                            for j in range(n + 1)]))
    return float(log_gamma - power * math.log(u0) -
                 2 * math.log(-math.expm1(-u0)) - math.log(k0) -
                 math.log(scale) - math.log(RW))


def controls():
    """Independent semi-infinite quadrature, with exp(-u0) scaled out."""
    records = []
    for u0 in (1., 5., 10., 100., 1000., 8000.):
        for power in (0, 1, 4, 8):
            k0, scale = 2.3, 1e-6
            f = lambda y: (1 + y / u0)**(4 + power) * math.exp(-y) / (-math.expm1(-u0-y))**2
            integral, error = quad(f, 0, np.inf, epsabs=1e-10,
                                   epsrel=2e-12, limit=200)
            actual_log = (-u0 + 4 * math.log(u0) + math.log(integral) -
                          math.log(k0 * scale * RW))
            upper_log = log_tail_upper(u0, k0, power, scale)
            excess = math.expm1(upper_log - actual_log)
            if excess < -2e-11 or error / integral > 1e-8:
                raise ValueError('tail bound failed independent quadrature')
            # Changing an opacity normalization must change its inverse exactly.
            scaling_error = abs((log_tail_upper(u0, k0, power, 1.) -
                                 upper_log) - math.log(scale))
            if scaling_error > 3e-12:
                raise ValueError('incorrect tail normalization dependence')
            records.append({'photon_max_u': u0, 'power': power,
                            'log_inverse_tail_integral': actual_log,
                            'log_inverse_tail_upper': upper_log,
                            'relative_bound_excess': excess,
                            'quadrature_relative_error_estimate': error / integral,
                            'log_normalization_error': scaling_error})
    for args in ((0., 1., 4, 1.), (1., -1., 4, 1.), (1., 1., 9, 1.),
                 (1., 1., 4, 0.), (math.inf, 1., 4, 1.)):
        try:
            log_tail_upper(*args)
        except ValueError:
            continue
        raise ValueError('invalid envelope was not rejected')
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bounds', type=Path)
    parser.add_argument('scratch', type=Path, help='completed bounds scratch directory')
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('use a new diagnostic result')
    old = json.loads(args.bounds.read_text())
    inputs = dict(old['input_sha256'])
    for name, expected in old['scratch_sha256'].items():
        add_inputs(inputs, {str((args.scratch / name).resolve()): expected})
    for path in (args.bounds, Path(__file__), Path('scripts/tops_groups.py'),
                 Path('scripts/reduce_tops_group_factors.py')):
        add_inputs(inputs, {str(path.resolve()): digest(path)})
    verify(inputs)
    checks = controls()
    # Parse each saved response once. No network requests or repeated conduction.
    sources = {root: source_groups(root) for root in
               sorted({r['source'] for r in old['records']})}
    rows = []
    for row in old['records']:
        t, rho = row['temperature_keV'], row['density_atomic_g_cm3']
        boundaries, groups, _ = sources[row['source']]
        u0 = float(boundaries[-1] / t)
        if u0 != row['photon_max_u'] or row['cutoff_u'] >= u0:
            raise ValueError('inconsistent photon support')
        # K0 is a group mean used ONLY as the scale of the stated envelope.
        # A group mean is not a measured endpoint or a lower bound beyond it.
        k0 = float(groups[t, rho][-1, 1])
        kc = row['conductive_opacity_atomic_cm2_g']
        normalization = row['native_uncut_normalization']
        band_upper = row['inverse_opacity_upper'] / normalization
        band_lower = row['inverse_opacity_lower'] / normalization
        cases = []
        for scale in (1., 1e-6, 1e-12):
            log_tail = log_tail_upper(u0, k0, 4, scale) - math.log(normalization)
            log_impact = log_tail - math.log(1 / kc + band_lower)
            log_omit = float(np.logaddexp(math.log(band_upper), log_tail) + math.log(kc))
            omission = math.exp(log_omit)
            cases.append({'opacity_envelope_scale': scale,
                          'log10_inverse_tail_upper': log_tail / math.log(10),
                          'log10_combined_opacity_relative_tail_upper': log_impact / math.log(10),
                          'radiation_omission_relative_opacity_upper': omission,
                          'radiative_fraction_upper': omission / (1 + omission)})
        rows.append({'temperature_keV': t, 'density_atomic_g_cm3': rho,
                     'photon_max_keV': float(boundaries[-1]), 'photon_max_u': u0,
                     'last_group_rosseland_atomic_cm2_g': k0,
                     'native_uncut_normalization': normalization,
                     'cases': cases})
    verify(inputs)
    result = {
        'scope': __doc__, 'accepted_for_stellar_opacity': False,
        'source_X_atomic': old['source_X_atomic'],
        'source_Z_atomic': old['source_Z_atomic'],
        'dense_source_states': len(rows),
        'tail_opacity_envelope': {
            'inequality': 'kappa(u) >= scale * last_group_Rosseland * (u/u_max)^(-4)',
            'domain': 'u >= u_max; positive total transport extinction',
            'normalizations': [1., 1e-6, 1e-12],
            'status': 'Explicit sensitivity assumption, not inferred from the last group mean.',
            'physical_context': 'High-energy scattering, absorption and pair production provide extinction; no dense-plasma lower bound is imported from neutral-atom tables.',
            'context_source': 'https://physics.nist.gov/PhysRefData/Xcom/Text/intro.html'},
        'independent_controls': checks,
        'negative_controls_passed': 5,
        'minimum_photon_max_u': min(r['photon_max_u'] for r in rows),
        'maximum_log10_combined_opacity_relative_tail_upper': max(
            c['log10_combined_opacity_relative_tail_upper'] for r in rows for c in r['cases']),
        'radiation_omission_relative_opacity_below_threshold_counts': {
            str(threshold): sum(r['cases'][-1]['radiation_omission_relative_opacity_upper'] < threshold for r in rows)
            for threshold in (.01, .005, .001, .0001)},
        'limitations': [
            'No omitted high-energy opacity is measured or reconstructed here.',
            'Tail bounds are conditional on the stated opacity envelope; scale changes are sensitivity tests.',
            'Finite-band bounds retain the specified collisionless refractive model and conductive table.',
            'The imposed native-uncut normalization is retained explicitly.',
            'Source nodes at X=0,Z=0.02 do not establish coverage between nodes or for other mixtures.',
            'Relative opacity change on omitting radiation equals kappa_conduction/kappa_radiation, not the radiative flux fraction.',
            'Logarithmic tail values remain finite even where their exponential would underflow.'],
        'input_sha256': inputs, 'records': rows}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in
                     ('scope', 'tail_opacity_envelope', 'independent_controls',
                      'limitations', 'input_sha256', 'records')}), flush=True)


if __name__ == '__main__':
    main()
