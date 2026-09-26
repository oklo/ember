#!/usr/bin/env python3
"""Compare interpolated refractive ratios with saved direct middle-Z sources.

The runtime evaluates each bracketing table; log-linear metallicity averaging
then follows MixtureOpacity. The reference composition was not a node of these
brackets. These checks include the chosen density compression, but do not
measure between-native-temperature errors or accept a stellar opacity family.
"""
import argparse
import json
import math
from pathlib import Path

from audit_tops_group_factor_tables import query
from fetch_tops_composition import digest
from reduce_tops_group_factors import add_inputs, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('lower', 'upper', 'reference', 'probe', 'probe_reference', 'scratch', 'output'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    if args.scratch.exists() or args.output.exists():
        raise FileExistsError('use new comparison outputs')
    lo, hi, ref, build = [json.loads(p.read_text()) for p in
                           (args.lower, args.upper, args.reference, args.probe_reference)]
    if not (lo['X'] == hi['X'] == ref['X'] and lo['Z'] < ref['Z'] < hi['Z']):
        raise ValueError('invalid composition bracket')
    if lo['temperatures_keV'] != hi['temperatures_keV']:
        raise ValueError('different bracketing temperature grids')
    inputs = {}
    for report in (lo, hi, ref, build):
        add_inputs(inputs, report['input_sha256'])
        add_inputs(inputs, report.get('output_sha256', {}))
    for path in (args.lower, args.upper, args.reference, args.probe_reference,
                 Path(__file__), Path('scripts/audit_tops_group_factor_tables.py')):
        add_inputs(inputs, {str(path.resolve()): digest(path)})
    if inputs.get(str(args.probe.resolve())) != digest(args.probe):
        raise ValueError('runtime executable is not identified')
    verify(inputs)
    args.scratch.mkdir(parents=True)
    tt = lo['temperatures_keV']
    selected, omitted = [], []
    for row in ref['records']:
        t, rho = row['temperature_keV'], row['density_atomic_g_cm3']
        if t not in tt:
            omitted.append({'temperature_keV': t, 'density_atomic_g_cm3': rho,
                            'reason': 'outside bracketing temperature grid'})
            continue
        it = tt.index(t)
        limits = []
        for table in (lo, hi):
            count = min(table['density_prefix_sizes'][max(0, it-2):min(len(tt), it+3)])
            limits.append((table['densities_atomic_g_cm3'][0], table['densities_atomic_g_cm3'][count-1]))
        if not max(v[0] for v in limits) <= rho <= min(v[1] for v in limits):
            omitted.append({'temperature_keV': t, 'density_atomic_g_cm3': rho,
                            'reason': 'outside conservative value/derivative support intersection'})
        else:
            selected.append(row)
    if not selected:
        raise ValueError('no supported reference states')
    points = [(r['temperature_keV'], r['density_atomic_g_cm3']) for r in selected]
    values = []
    for label, table in (('lower', lo), ('upper', hi)):
        path = next(Path(p) for p in table['output_sha256'] if Path(p).name == 'factor.dat')
        values.append(query(args.probe, path, points, ref['X'], args.scratch, label))
    f = (ref['Z'] - lo['Z']) / (hi['Z'] - lo['Z'])
    rows = []
    for row, low, high in zip(selected, *values, strict=True):
        direct = row['refractive_ratio']
        actual = math.exp((1-f)*math.log(low[0]) + f*math.log(high[0]))
        rows.append({'temperature_keV': row['temperature_keV'],
                     'density_atomic_g_cm3': row['density_atomic_g_cm3'],
                     'direct_ratio': direct, 'interpolated_ratio': actual,
                     'relative_error': actual / direct - 1,
                     'on_both_bracket_density_grids': all(row['density_atomic_g_cm3'] in
                         table['densities_atomic_g_cm3'] for table in (lo, hi))})
    maximum = max(abs(r['relative_error']) for r in rows)
    hot = [r for r in rows if r['temperature_keV'] >= .5 and r['density_atomic_g_cm3'] >= 1000]
    verify(inputs)
    result = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
              'X': ref['X'], 'reference_Z': ref['Z'], 'bracket_Z': [lo['Z'], hi['Z']],
              'source_comparisons': len(rows), 'reference_states_omitted': len(omitted),
              'comparisons_off_bracket_density_grid': sum(not r['on_both_bracket_density_grids'] for r in rows),
              'maximum_relative_error': maximum, 'relative_criterion': .005,
              'sampled_comparison_passed': maximum <= .005,
              'hot_dense_comparisons': len(hot),
              'hot_dense_maximum_relative_error': max((abs(r['relative_error']) for r in hot), default=None),
              'maximum_relative_error_by_temperature': {str(t): max(abs(r['relative_error']) for r in rows
                   if r['temperature_keV'] == t) for t in sorted({r['temperature_keV'] for r in rows})},
              'limitations': ['Reference sources are saved direct group calculations, not newly fetched spectra.',
                              'Density compression was chosen using earlier X=0.2,Z=0.02 results.',
                              'This joint Z/density comparison does not isolate either source of error.',
                              'It does not test interpolation between native temperature rows.',
                              'Ratios do not replace native uncut grey means; their errors must be considered together.'],
              'input_sha256': inputs, 'comparisons': rows, 'omitted': omitted,
              'scratch_sha256': {p.name: digest(p) for p in args.scratch.iterdir() if p.is_file()}}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in
                     ('scope', 'limitations', 'input_sha256', 'comparisons', 'omitted',
                      'scratch_sha256', 'maximum_relative_error_by_temperature')}), flush=True)
    if not result['sampled_comparison_passed']:
        raise ValueError('joint metallicity/density comparison exceeds its retained criterion')


if __name__ == '__main__':
    main()
