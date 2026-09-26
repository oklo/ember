#!/usr/bin/env python3
"""Compare candidate refractive ratios with independent composition sources.

Ember evaluates the bracketing single-composition tables in temperature and
density. Log-linear averaging in X and Z follows its composition convention.
This tests joint X/Z/density interpolation at native temperatures; it does not
test a resampled final opacity family or between-native-temperature accuracy.
"""
import argparse
import bisect
import json
import math
from pathlib import Path

from audit_tops_group_factor_tables import query
from fetch_tops_composition import digest
from reduce_tops_group_factors import add_inputs, verify


def bracket(axis, value):
    if not axis[0] <= value <= axis[-1]:
        raise ValueError('composition outside candidate family')
    i = min(max(0, bisect.bisect_right(axis, value) - 1), len(axis) - 2)
    return axis[i], axis[i + 1], (value - axis[i]) / (axis[i + 1] - axis[i])


def supported(table, t, rho):
    tt = table['temperatures_keV']
    if t not in tt:
        return False
    i = tt.index(t)
    count = min(table['density_prefix_sizes'][max(0, i - 2):min(len(tt), i + 3)])
    rr = table['densities_atomic_g_cm3']
    return rr[0] <= rho <= rr[count - 1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate', type=Path, nargs='+', required=True)
    p.add_argument('--reference', type=Path, nargs='+', default=[])
    p.add_argument('--saved-reference', type=Path, action='append', default=[],
                   help='Reuse factor_records from a saved composition comparison')
    for name in ('probe', 'probe_reference', 'scratch', 'output'):
        p.add_argument('--' + name.replace('_', '-'), type=Path, required=True)
    a = p.parse_args()
    if a.scratch.exists() or a.output.exists():
        raise FileExistsError('preserve previous comparison outputs')
    candidates, refs, inputs = {}, [], {}
    for role, paths in (('candidate', a.candidate), ('reference', a.reference)):
        for path in paths:
            report = json.loads(path.read_text())
            add_inputs(inputs, report['input_sha256'])
            add_inputs(inputs, report['output_sha256'])
            add_inputs(inputs, {str(path.resolve()): digest(path)})
            if role == 'candidate':
                key = report['X'], report['Z']
                if key in candidates:
                    raise ValueError('duplicate candidate composition')
                candidates[key] = report
            else:
                refs.append(report)
    xs = sorted({x for x, z in candidates})
    zs = sorted({z for x, z in candidates})
    if len(xs) < 2 or len(zs) < 2 or set(candidates) != {(x, z) for x in xs for z in zs}:
        raise ValueError('incomplete candidate composition rectangle')
    for path in a.saved_reference:
        report = json.loads(path.read_text())
        add_inputs(inputs, report['input_sha256'])
        add_inputs(inputs, {str(path.resolve()): digest(path)})
        grouped = {}
        for row in report['factor_records']:
            x, z = row['X'], row['Z']
            if x in xs or z in zs:
                continue
            grouped.setdefault((x, z), []).append({
                'temperature_keV': row['temperature_keV'],
                'density_atomic_g_cm3': row['density_atomic_g_cm3'],
                'refractive_ratio': row['factor'], 'source': row['source']})
        refs.extend({'X': x, 'Z': z, 'records': rows} for (x, z), rows in grouped.items())
    if not refs:
        raise ValueError('no independent reference compositions supplied')
    if len({(r['X'], r['Z']) for r in refs}) != len(refs):
        raise ValueError('duplicate reference composition')
    build = json.loads(a.probe_reference.read_text())
    add_inputs(inputs, build['input_sha256'])
    for path in (a.probe_reference, Path(__file__),
                 Path('scripts/audit_tops_group_factor_tables.py')):
        add_inputs(inputs, {str(path.resolve()): digest(path)})
    if inputs.get(str(a.probe.resolve())) != digest(a.probe):
        raise ValueError('unidentified runtime probe')
    verify(inputs)
    selected, omitted, requests = [], [], {key: set() for key in candidates}
    for ref in refs:
        x, z = ref['X'], ref['Z']
        if x in xs or z in zs:
            raise ValueError('reference must lie between both X and Z nodes')
        xl, xh, fx = bracket(xs, x)
        zl, zh, fz = bracket(zs, z)
        corners = [(xl, zl), (xh, zl), (xl, zh), (xh, zh)]
        for row in ref['records']:
            t, rho = row['temperature_keV'], row['density_atomic_g_cm3']
            point = t, rho
            record = {'X': x, 'Z': z, 'temperature_keV': t,
                      'density_atomic_g_cm3': rho, 'direct_ratio': row['refractive_ratio'],
                      'reference_source': row['source']}
            if not all(supported(candidates[k], t, rho) for k in corners):
                omitted.append({**record, 'reason': 'outside conservative support intersection'})
                continue
            selected.append((record, corners, fx, fz))
            for key in corners:
                requests[key].add(point)
    if not selected:
        raise ValueError('no supported comparisons')
    a.scratch.mkdir(parents=True)
    values = {}
    for (x, z), points in requests.items():
        if not points:
            continue
        table = candidates[x, z]
        path = next(Path(p) for p in table['output_sha256'] if Path(p).name == 'factor.dat')
        points = sorted(points)
        rows = query(a.probe, path, points, x, a.scratch, f'x{x:g}-z{z:g}')
        values[x, z] = {point: row[0] for point, row in zip(points, rows, strict=True)}
    comparisons = []
    for record, corners, fx, fz in selected:
        point = record['temperature_keV'], record['density_atomic_g_cm3']
        ll, hl, lh, hh = [math.log(values[key][point]) for key in corners]
        actual = math.exp((1 - fz) * ((1 - fx) * ll + fx * hl)
                          + fz * ((1 - fx) * lh + fx * hh))
        off_grid = any(point[1] not in candidates[k]['densities_atomic_g_cm3'] for k in corners)
        comparisons.append({**record, 'interpolated_ratio': actual,
                            'relative_error': actual / record['direct_ratio'] - 1,
                            'off_candidate_density_grid': off_grid})
    maximum = max(abs(r['relative_error']) for r in comparisons)
    per_x = {}
    for x in sorted({r['X'] for r in comparisons}):
        rows = [r for r in comparisons if r['X'] == x]
        dense = [r for r in rows if r['density_atomic_g_cm3'] > 10000]
        per_x[str(x)] = {'comparisons': len(rows),
                         'maximum_relative_error': max(abs(r['relative_error']) for r in rows),
                         'dense_comparisons': len(dense),
                         'dense_maximum_relative_error': max((abs(r['relative_error']) for r in dense), default=None)}
    verify(inputs)
    result = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
              'candidate_X': xs, 'candidate_Z': zs, 'source_comparisons': len(comparisons),
              'reference_reports': [str(p) for p in a.reference + a.saved_reference],
              'reference_states_omitted': len(omitted), 'comparison_by_X': per_x,
              'comparisons_off_density_grid': sum(r['off_candidate_density_grid'] for r in comparisons),
              'relative_criterion': .005, 'maximum_relative_error': maximum,
              'sampled_comparison_passed': maximum <= .005,
              'limitations': ['Reference compositions are not candidate nodes.',
                              'This measures joint composition and density interpolation at native temperatures.',
                              'The final grey-times-ratio opacity family requires its own runtime checks.',
                              'Source-normalized output does not repair failed strict group/native recovery diagnostics.'],
              'input_sha256': inputs, 'comparisons': comparisons, 'omitted': omitted,
              'scratch_sha256': {p.name: digest(p) for p in a.scratch.iterdir() if p.is_file()}}
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in
                     ('scope', 'limitations', 'input_sha256', 'comparisons', 'omitted', 'scratch_sha256')}), flush=True)
    if not result['sampled_comparison_passed']:
        raise ValueError('joint composition/density comparison exceeds retained criterion')


if __name__ == '__main__':
    main()
