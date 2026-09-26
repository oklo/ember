#!/usr/bin/env python3
"""Assess OP source density support without evaluating or selecting an opacity.

At each native (T, ne), charge neutrality gives rho = m_u ne / sum(w_k q_k),
where w_k is atoms per reference mass in atomic mass units. Stellar weights
use the tracked integer isotope masses. Missing P, Cl, K, Ti contribute a
reported neutral-to-fully-ionized charge interval; this is a conditional
coverage screen, not a bound on their omitted opacity or negative-ion physics.
"""
import csv
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np

MU = 1.66053906660e-24
KB = 1.380649e-16
KEV = 1.602176634e-9
ROOT = Path('/Users/greglaughlin/Projects/ember')
SMALL = Path('/tmp/ember-op-native-grid-inspection-v1')
SUMMARY = SMALL / 'OP4STARS_1.3/mono'
PROFILE = ROOT / 'docs/reports/2026-09-11/evolution_latest_profile.csv'
CONTROL = ROOT / 'docs/results/tops_cool_full_spectra_v4.json'
METALS = ROOT / 'include/ember/gs98_mixture.hpp'
OUTPUT = ROOT / 'docs/results/op_native_stellar_domain_v2.json'


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load_summaries():
    manifest = json.loads((SMALL / 'extraction_manifest.json').read_text())
    expected = {r['path']: r['sha256'] for r in manifest['files']}
    data = {}; inputs = {}
    for p in sorted(SUMMARY.glob('m*.smry')):
        inputs[str(p)] = digest(p)
        if inputs[str(p)] != expected[str(p)]:
            raise ValueError('native summary changed')
        lines = iter(p.read_text().splitlines())
        header = next(lines).split(); element = int(header[0])
        start, stop, step = map(int, next(lines).split()); grid = {}
        for it in range(start, stop + 1, step):
            label, j0, j1, dj = map(int, next(lines).split())
            if label != it:
                raise ValueError('temperature label mismatch')
            rows = np.array([list(map(float, next(lines).split()))
                             for _ in range(j0, j1 + 1, dj)])
            if not np.array_equal(rows[:, 0], np.arange(j0, j1 + 1, dj)):
                raise ValueError('density labels mismatch')
            if not np.isfinite(rows).all() or np.any(rows[:, 2] < 0) or np.any(rows[:, 3] <= 0):
                raise ValueError('unusable summary means')
            grid[it] = rows
        if list(lines):
            raise ValueError('trailing summary data')
        data[element] = grid
    if len(data) != 17:
        raise ValueError('incomplete element summaries')
    its = np.array(sorted(data[1]))
    for grids in data.values():
        if not np.array_equal(sorted(grids), its):
            raise ValueError('different temperature grids')
        for it in its:
            if not np.array_equal(grids[it][:, 0], data[1][it][:, 0]):
                raise ValueError('different electron density grids')
    return data, its, inputs


def reference_weights(row, metals):
    x, y3, y4 = (float(row[k]) for k in ('X', 'Y3', 'Y4'))
    z = 1 - x - y3 - y4
    if min(x, y3, y4, z) < 0:
        raise ValueError('invalid saved composition')
    weights = {1: x, 2: y3 / 3 + y4 / 4}
    weights.update({int(e[0]): z * e[3] / e[1] for e in metals})
    return weights, z


def support(data, its, weights, temperature, density, count):
    logt = math.log10(temperature); target = logt / .025
    upper = int(np.searchsorted(its, target, side='right'))
    if target < its[0] or target > its[-1]:
        return dict(supported=False, reason='temperature_outside_native_grid')
    # The stencil is contained in the source grid; requested T is never clipped.
    lower = min(max(upper - 1 - (count - 2) // 2, 0), len(its) - count)
    selected = list(map(int, its[lower:lower + count]))
    missing = {z: w for z, w in weights.items() if z not in data}
    missing_max = sum(z * w for z, w in missing.items())
    planes = []
    for it in selected:
        indices = data[1][it][:, 0]
        ne = 10 ** (.25 * indices)
        q = sum(w * data[z][it][:, 1] for z, w in weights.items() if z in data)
        # The complete n_e interval must map to positive, ordered densities for
        # every charge assumption here. Do not bridge an invalid interval.
        valid = q > 0
        where = np.flatnonzero(valid)
        if len(where) < 2 or not np.array_equal(where, np.arange(where[0], where[-1] + 1)):
            return dict(supported=False, reason='nonpositive_or_disconnected_charge', temperature_indices=selected)
        r_neutral = MU * ne[where] / q[where]
        r_ionized = MU * ne[where] / (q[where] + missing_max)
        if np.any(np.diff(r_neutral) <= 0) or np.any(np.diff(r_ionized) <= 0):
            return dict(supported=False, reason='nonmonotonic_density_mapping', temperature_indices=selected)
        planes.append(dict(temperature_index=it, log_temperature=.025 * it,
                           electron_index_min=int(indices[where[0]]), electron_index_max=int(indices[where[-1]]),
                           density_min_interval=[float(r_ionized[0]), float(r_neutral[0])],
                           density_max_interval=[float(r_ionized[-1]), float(r_neutral[-1])],
                           missing_charge_fraction_at_dense_edge=float(missing_max / (q[where[-1]] + missing_max))))
    lo = max(p['density_min_interval'][1] for p in planes)
    hi = min(p['density_max_interval'][0] for p in planes)
    # A sufficient rectangle shared by every native T plane and both stated
    # missing-element charge limits. This does not reproduce cubic OPserver.
    okay = lo <= density <= hi
    return dict(supported=okay, reason='inside_conditional_rectangle' if okay else
                ('density_above_rectangle' if density > hi else 'density_below_rectangle'),
                temperature_indices=selected, shared_density_min=lo, shared_density_max=hi,
                requested_density_over_upper_limit=density / hi,
                missing_charge_per_reference_mass=missing_max, planes=planes)


def main():
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    data, its, inputs = load_summaries()
    for p in (Path(__file__).resolve(), PROFILE, CONTROL, METALS, SMALL / 'extraction_manifest.json'):
        inputs[str(p)] = digest(p)
    metals = [tuple(map(float, row.split(','))) for row in
              re.findall(r'^\s*\{([^{}]+)\}, //', METALS.read_text(), re.M)]
    if len(metals) != 19 or abs(sum(e[3] for e in metals) - 1) > 2e-15:
        raise ValueError('GS98 composition source changed')
    # Independent fully ionized limits check both isotope accounting and the
    # reference-mass convention, using a native high-temperature dilute row.
    analytic = []
    for x, y3, y4 in ((1, 0, 0), (0, 0, 1), (0, 1, 0), (.1, .03, .87)):
        weights, _ = reference_weights(dict(X=x, Y3=y3, Y4=y4), metals)
        weights = {z: w for z, w in weights.items() if w}
        q = sum(w * data[z][320][0, 1] for z, w in weights.items())
        expected = x + 2 * y3 / 3 + y4 / 2
        assert abs(q / expected - 1) < 1e-14
        analytic.append(dict(X=x, Y3=y3, Y4=y4, electron_weight=q, expected=expected))
    with PROFILE.open() as f:
        profile = list(csv.DictReader(f))
    records = []
    for zone, row in enumerate(profile):
        weights, z = reference_weights(row, metals)
        t = float(row['temperature_K']); rho = float(row['density_g_cm3'])
        records.append(dict(zone=zone, temperature_K=t, density_baryonic_g_cm3=rho,
                            X=float(row['X']), Y3=float(row['Y3']), Y4=float(row['Y4']), Z=z,
                            mass_g=float(row['mass_g']), missing_metal_mass_fraction=z * sum(e[3] for e in metals if int(e[0]) not in data),
                            two_temperature=support(data, its, weights, t, rho, 2),
                            four_temperature=support(data, its, weights, t, rho, 4)))
    # TOPS controls retain their printed atom-number fractions and source mass
    # calibration; do not substitute the integer-isotope stellar convention.
    refs = json.loads(CONTROL.read_text()); sources = {}; comparisons = []
    for ref in refs['records']:
        root = Path(ref['full_spectrum_source'])
        if str(root) not in sources:
            receipt = json.loads((root / 'receipt.json').read_text())
            if digest(root / 'source.txt') != receipt['sha256']:
                raise ValueError('TOPS spectral source changed')
            text = (root / 'source.txt').read_text()
            body = text.split('No. Fraction Mass Fraction  At. No.  Chem. Sym.  Mat ID.\n')[1].split('Temperature grid')[0]
            elements = [r.split() for r in body.splitlines() if len(r.split()) == 5]
            helium = next(r for r in elements if r[3] == 'He')
            mean_mass = 4.002602 * float(helium[0]) / float(helium[1])
            sources[str(root)] = {int(r[2]): float(r[0]) / mean_mass for r in elements}
            for p in (root / 'source.txt', root / 'receipt.json'):
                inputs[str(p)] = digest(p)
        t = ref['temperature_keV'] * KEV / KB; rho = ref['density_atomic_g_cm3']
        comparisons.append(dict(X=ref['X'], Z=ref['Z'], temperature_keV=ref['temperature_keV'],
                                density_atomic_g_cm3=rho, source=str(root),
                                two_temperature=support(data, its, sources[str(root)], t, rho, 2),
                                four_temperature=support(data, its, sources[str(root)], t, rho, 4)))
    counts = {key: dict(supported=sum(r[key]['supported'] for r in records),
                        unsupported_zones=[r['zone'] for r in records if not r[key]['supported']],
                        supported_cool_layers=sum(r[key]['supported'] and r['temperature_K'] < .05 * KEV / KB for r in records),
                        cool_layers=sum(r['temperature_K'] < .05 * KEV / KB for r in records))
              for key in ('two_temperature', 'four_temperature')}
    for p, h in inputs.items():
        if digest(p) != h:
            raise ValueError('input changed during domain assessment')
    result = dict(scope=__doc__, outcome='completed_conditional_domain_assessment', accepted_for_stellar_opacity=False,
                  input_sha256=inputs, analytic_charge_controls=analytic, stellar_summary=counts,
                  zero_printed_planck_counts={str(z): sum(int(np.count_nonzero(r[:, 2] == 0)) for r in grid.values()) for z, grid in data.items()},
                  stellar_records=records, tops_reference_cases=comparisons,
                  tops_supported_counts={k: sum(r[k]['supported'] for r in comparisons) for k in counts},
                  limitations=['The original parser rejected zero printed Planck means in five elements. The domain calculation uses electron counts only; v2 records those zeros, retains positive-Rosseland and finite-value checks, and makes no Planck-opacity acceptance claim.',
                               'Density support is a sufficient rectangle on native temperature planes; failure of this rectangle alone need not exclude a different interpolation.',
                               'No native cross section is interpolated, no selected opacity is changed, and no missing elemental spectrum is approximated.',
                               'Missing-element charge interval assumes neutral through fully ionized populations; negative ions are not bounded by that interval.',
                               'Printed summary electron counts have four significant figures; these are source-domain diagnostics, not precision EOS comparisons.',
                               'No molecular/grain completeness, dense-plasma correction or spectral tail accuracy is established.'])
    OUTPUT.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(dict(stellar=counts, tops_supported_counts=result['tops_supported_counts'], tops_cases=len(comparisons))), flush=True)


if __name__ == '__main__':
    main()
