#!/usr/bin/env python3
"""Compare compact group transport with verified full-spectrum hot opacities."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from audit_tops_spectral_means import digest
from audit_tops_electron_dispersion import electron_moments, index_squared, KEV, KB, HBAR, ME, RW
from tops_groups import source_groups, transport


def constant_control():
    t, ne, absorption, scattering = 1.25, 1e26, 10., 2.
    e = electron_moments(t*KEV/KB, ne)
    up = HBAR*math.sqrt(4*math.pi*4.80320471257e-10**2*ne/ME)/(t*KEV)
    up *= math.sqrt(e['plasma_frequency_squared_ratio'])
    def f(u):
        n = math.sqrt(float(index_squared(np.array([u]), up, e['vstar_squared'])[0]))
        return u**4*math.exp(-u)/(-math.expm1(-u))**2*n**3/(absorption+n*scattering)
    reference = RW/quad(f, up, 200/t, epsabs=1e-11, epsrel=1e-11, points=[1., 4., 10., 30.])[0]
    boundaries = np.geomspace(.01, 200., 129)
    groups = np.column_stack([boundaries[:-1], np.full(128, absorption+scattering), np.full(128, absorption)])
    value = transport(boundaries, groups, t, ne)['rosseland_component_estimate']
    error = abs(value/reference-1)
    if error > 1e-8:
        raise ValueError('constant-component quadrature control failed')
    return {'independent_quadrature_relative_error': error, 'reference': reference}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'reference', 'output'):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    if a.output.exists(): raise FileExistsError('use a new group report')
    spec = json.loads(a.plan.read_text()); old = json.loads(a.reference.read_text())
    inputs = dict(old['input_sha256'])
    paths = [a.plan, a.reference, Path(__file__), Path(__file__).with_name('tops_groups.py'),
             Path(__file__).with_name('fetch_tops_group_plan.py')]
    inputs.update({str(path.resolve()): digest(path) for path in paths})
    for name, h in inputs.items():
        if digest(Path(name)) != h: raise ValueError('group reference dependency changed: '+name)
    reference = {(r['temperature_keV'], r['density_atomic_g_cm3']): r for r in old['records']}
    rows, seen, source_bytes = [], set(), 0
    control = constant_control()
    for job in spec['requests']:
        root = Path(job['work']); boundaries, data, grey = source_groups(root)
        inputs.update({str(path.resolve()): digest(path) for path in
                       [root/n for n in ('receipt.json', 'request.json', 'source.txt', 'recipe.json')]})
        source_bytes += (root/'source.txt').stat().st_size
        for key, groups in data.items():
            if key in seen: raise ValueError('duplicate group state')
            seen.add(key)
            pair = [transport(boundaries, groups, key[0], grey[key]['electron_density_cm3'], order) for order in (16, 32)]
            error = abs(pair[0]['rosseland_component_estimate']/pair[1]['rosseland_component_estimate']-1)
            if error > 1e-6: raise ValueError('group quadrature did not converge')
            row = {'temperature_keV': key[0], 'density_atomic_g_cm3': key[1],
                   'source': str(root), 'electron_density_cm3': grey[key]['electron_density_cm3'],
                   'quadrature_relative_change': error, **pair[1],
                   'full_spectrum_reference_available': key in reference}
            if key in reference:
                ref = reference[key]
                row.update(reference_rosseland=ref['rosseland_atomic_cm2_g'],
                    component_relative_error=pair[1]['rosseland_component_estimate']/ref['rosseland_atomic_cm2_g']-1,
                    common_index_relative_error=pair[1]['rosseland_common_index']/ref['rosseland_atomic_cm2_g']-1,
                    electron_density_relative_difference=grey[key]['electron_density_cm3']/ref['free_electron_density_cm3']-1)
            rows.append(row)
    if not set(reference) <= seen:
        raise ValueError('group plan does not cover every reference state')
    comparisons = [r for r in rows if r['full_spectrum_reference_available']]
    maximum = max(abs(r['component_relative_error']) for r in comparisons)
    if any(digest(Path(name)) != h for name, h in inputs.items()):
        raise ValueError('group input changed during the comparison')
    result = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
        'constant_component_control': control, 'source_states': len(rows), 'full_spectrum_comparisons': len(comparisons),
        'groups_per_state': spec['photon_boundaries']-1, 'source_text_bytes': source_bytes,
        'maximum_component_relative_error': maximum, 'relative_comparison_criterion': .005,
        'group_representation_check_passed': maximum <= .005,
        'maximum_quadrature_relative_change': max(r['quadrature_relative_change'] for r in rows),
        'records': rows, 'input_sha256': inputs,
        'assumptions': ['Group Rosseland means are harmonic averages of total index-one opacity with the source thermal weight.',
            'Within a narrow group, max(0, 1-Planck_absorption/Rosseland_total) estimates the scattering fraction. It is not an independently supplied scattering mean.',
            'The correction retains the source group Rosseland inverse mean when the index is one, and evaluates the index within each group.',
            'The finite-temperature free-electron dispersion and full thermal normalization are unchanged from the full-spectrum calculation.',
            'Agreement at these hot X=0, Z=0.02 states does not validate other temperatures, compositions or densities.']}
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('records', 'input_sha256', 'assumptions')}), flush=True)


if __name__ == '__main__': main()
