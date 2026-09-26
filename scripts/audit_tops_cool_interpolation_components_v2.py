#!/usr/bin/env python3
"""Separate grey, density-ratio and composition-ratio interpolation differences.

Use new ratio planes at the exact reference compositions, with reference
densities withheld, to isolate density interpolation. Compare native uncut
grey tables independently. The coarse joint ratio comparison includes both its coarser density grid and composition interpolation; it does not isolate composition alone. The product of the refined grey and ratio interpolants is a
diagnostic, not a resampled stellar opacity family or a physical acceptance.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from audit_tops_group_factor_tables import query
from fetch_tops_composition import digest
from reduce_tops_group_factors import add_inputs, verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scratch', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    a = p.parse_args()
    if a.scratch.exists() or a.report.exists():
        raise FileExistsError('preserve completed interpolation comparisons')
    inputs = {}

    def pin(path):
        path = Path(path)
        add_inputs(inputs, {str(path.resolve()): digest(path)})
        return path

    def report(path):
        r = json.loads(pin(path).read_text())
        add_inputs(inputs, r['input_sha256'])
        add_inputs(inputs, r.get('output_sha256', {}))
        return r

    joint = report('docs/results/tops_cool_ratio_hydrogen_comparison_v1.json')
    probe = Path('/tmp/ember-tops-spectral-runtime-probe-v1')
    if inputs.get(str(probe.resolve())) != digest(probe):
        raise ValueError('unidentified table probe')
    for name in ('audit_tops_group_factor_tables.py', 'reduce_tops_group_factors.py',
                 'audit_tops_cool_interpolation_components_v2.py'):
        pin(Path(__file__).with_name(name))
    selected = {}
    for label in ('010', '030', '070'):
        grey = report(f'docs/results/tops_cool_grey_refined_x{label}-z020_v1.json')
        factor = report(f'docs/results/tops_cool_density_refinement_x{label}-z020_reduction_v1.json')
        reference = report(f'docs/results/tops_cool_ratio_x{label}-z020-control_reduction_v1.json')
        if (grey['X'], grey['Z']) != (factor['X'], factor['Z']) or (factor['X'], factor['Z']) != (reference['X'], reference['Z']):
            raise ValueError('comparison compositions differ')
        selected[label] = grey, factor, reference
    verify(inputs)
    a.scratch.mkdir()
    previous = {(r['X'], r['Z'], r['temperature_keV'], r['density_atomic_g_cm3']): r
                for r in joint['comparisons']}
    records = []
    node_checks = []
    for label, (grey, factor, reference) in selected.items():
        x, z = factor['X'], factor['Z']
        refs = [r for r in reference['records'] if r['temperature_keV'] < .025]
        if len(refs) != 140:
            raise ValueError('expected all 140 new independent cold reference states')
        points = [(r['temperature_keV'], r['density_atomic_g_cm3']) for r in refs]
        grey_table = Path(next(iter(grey['output_sha256'])))
        factor_table = next(Path(p) for p in factor['output_sha256'] if Path(p).name == 'factor.dat')
        g = query(probe, grey_table, points, x, a.scratch, f'x{label}-grey-independent')
        f = query(probe, factor_table, points, x, a.scratch, f'x{label}-factor-independent')
        native = [r for r in factor['records'] if r['temperature_keV'] < .025]
        npoints = [(r['temperature_keV'], r['density_atomic_g_cm3']) for r in native]
        nvalues = query(probe, factor_table, npoints, x, a.scratch, f'x{label}-factor-nodes')
        error = max(abs(v[0]/r['refractive_ratio']-1) for v, r in zip(nvalues, native, strict=True))
        if error > 1e-10:
            raise ValueError('native factor values do not reproduce')
        node_checks.append({'X': x, 'queries': len(native), 'maximum_relative_error': error})
        for ref, gv, fv in zip(refs, g, f, strict=True):
            t, rho = ref['temperature_keV'], ref['density_atomic_g_cm3']
            old = previous[x, z, t, rho]
            records.append({'X': x, 'Z': z, 'temperature_keV': t,
                            'density_atomic_g_cm3': rho, 'reference_source': ref['source'],
                            'off_grey_density_grid': rho not in grey['densities_atomic_g_cm3'],
                            'off_factor_density_grid': rho not in factor['densities_atomic_g_cm3'],
                            'grey_relative_error': gv[0]/ref['uncut_source_rosseland']-1,
                            'density_only_ratio_relative_error': fv[0]/ref['refractive_ratio']-1,
                            'joint_ratio_relative_error': old['relative_error'],
                            'coarse_joint_ratio_relative_change': old['interpolated_ratio']/fv[0]-1,
                            'diagnostic_product_relative_error': gv[0]*fv[0]/ref['native_uncut_times_ratio_atomic_cm2_g']-1})
    fields = ('grey_relative_error', 'density_only_ratio_relative_error',
              'joint_ratio_relative_error', 'coarse_joint_ratio_relative_change',
              'diagnostic_product_relative_error')
    summaries = {}
    for name, rows in [('all_cold', records),
                       *[(f'rho_at_most_{limit:g}', [r for r in records if r['density_atomic_g_cm3'] <= limit])
                         for limit in (.01, .1, 1., 10.)]]:
        summaries[name] = {'comparisons': len(rows), **{
            field: {'maximum_absolute_difference': max(abs(r[field]) for r in rows),
                    'count_above_0p5_percent': sum(abs(r[field]) > .005 for r in rows)}
            for field in fields}}
    # Compare relevant stellar structures without rerunning their evolution.
    benchmark = json.loads(pin('docs/results/hermite_value_evolution_benchmark_v1.json').read_text())
    early_path = pin('/tmp/ember-hermite-value-benchmark-v1/run-1.json')
    if digest(early_path) != benchmark['runs'][1]['output_sha256']:
        raise ValueError('early-model output differs from completed benchmark')
    early = json.loads(early_path.read_text())
    if not early['converged'] or early['history'][-1][0] != 1e10:
        raise ValueError('unexpected early-model age or convergence')
    profile_path = pin('docs/reports/2026-09-11/evolution_latest_profile.csv')
    runtime = report('docs/results/tops_refractive_hot_family_refined_runtime_v1.json')
    if runtime['input_sha256'].get(str(profile_path.resolve())) != digest(profile_path):
        raise ValueError('unidentified evolved profile')
    profiles = [('early_10_Gyr', [dict(zip(early['profile_columns'], row, strict=True)) for row in early['profile']]),
                ('evolved_3890_Gyr', [{k: float(v) for k, v in row.items()} for row in csv.DictReader(profile_path.open())])]
    profile_rows = []
    for label, rows in profiles:
        rows = sorted(rows, key=lambda r: r['temperature_K'])
        tt = np.log([r['temperature_K'] for r in rows])
        rr = np.log([r['density_g_cm3'] for r in rows])
        for T in (3e4, 5e4, 1e5, 3e5, 5e5):
            if not tt[0] <= math.log(T) <= tt[-1]:
                raise ValueError('profile interpolation would extrapolate')
            profile_rows.append({'model': label, 'temperature_K': T,
                                 'density_baryonic_g_cm3': float(np.exp(np.interp(math.log(T), tt, rr)))})
    verify(inputs)
    result = {'scope': __doc__, 'analysis_complete': True, 'accepted_for_stellar_opacity': False,
              'reference_comparisons': len(records), 'node_checks': node_checks,
              'summaries': summaries, 'comparisons': records, 'profile_samples': profile_rows,
              'limitations': ['The early profile is at 10 Gyr, not the exact initial structure.',
                              'Profile densities use baryonic mass; opacity references use atomic source mass.',
                              'Density-only comparisons share composition with the reference; separate X/Z controls remain necessary.',
                              'The diagnostic product is not the final resampled opacity interpolation.',
                              'These comparisons test source-table interpolation, not absolute opacity physics or full-spectrum accuracy.'],
              'input_sha256': inputs,
              'output_sha256': {str(p.resolve()): digest(p) for p in a.scratch.iterdir() if p.is_file()}}
    a.report.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'summaries': summaries, 'profile_samples': profile_rows}), flush=True)


if __name__ == '__main__':
    main()
