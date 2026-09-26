#!/usr/bin/env python3
"""Check candidate group tables in Ember and compare saved independent spectra.

Supports nonzero hydrogen and temperature-dependent density limits. Tests
conversion, derivatives and source exclusions without accepting stellar inputs.
Optional density compression compares omitted source nodes, not fresh sources.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess

from audit_tops_electron_dispersion import KEV, KB
from audit_tops_spectral_means import digest
from reduce_tops_group_factors import add_inputs, verify, table_text


def query(probe, table, points, x, work, label, reject=False):
    request = ''.join(f'{x:.17g} {t * KEV / KB:.17g} {r:.17g}\n' for t, r in points)
    (work / (label + '.request.txt')).write_text(request)
    result = subprocess.run([str(probe), str(table), str(table)], input=request,
                            text=True, capture_output=True, timeout=60)
    (work / (label + '.stdout')).write_text(result.stdout)
    (work / (label + '.stderr')).write_text(result.stderr)
    if reject:
        if result.returncode == 0 or 'outside table' not in result.stderr:
            raise ValueError('unsupported point was not rejected by the domain guard')
        return None
    result.check_returncode()
    rows = [list(map(float, line.split())) for line in result.stdout.splitlines()]
    if len(rows) != len(points) or any(len(r) != 3 or not all(map(math.isfinite, r)) or r[0] <= 0 for r in rows):
        raise ValueError('invalid probe output')
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('reduction', 'probe', 'probe_reference', 'scratch', 'output'):
        p.add_argument(name, type=Path)
    p.add_argument('--spectral-reference', type=Path)
    p.add_argument('--compress-base-density', action='store_true')
    a = p.parse_args()
    if a.scratch.exists() or a.output.exists():
        raise FileExistsError('use new check outputs')
    old = json.loads(a.reduction.read_text())
    build = json.loads(a.probe_reference.read_text())
    inputs = dict(old['input_sha256'])
    add_inputs(inputs, old['output_sha256'])
    add_inputs(inputs, build['input_sha256'])
    for path in (a.reduction, a.probe_reference, Path(__file__),
                 Path(__file__).with_name('reduce_tops_group_factors.py')):
        add_inputs(inputs, {str(path.resolve()): digest(path)})
    if inputs.get(str(a.probe.resolve())) != digest(a.probe):
        raise ValueError('runtime executable identity is missing')
    spectral = None
    if a.spectral_reference:
        if old['X'] != 0 or old['Z'] != .02:
            raise ValueError('saved contraction spectra are X=0, Z=.02')
        spectral = json.loads(a.spectral_reference.read_text())
        add_inputs(inputs, spectral['input_sha256'])
        add_inputs(inputs, {str(a.spectral_reference.resolve()): digest(a.spectral_reference)})
    verify(inputs)
    a.scratch.mkdir(parents=True)
    tt, rr, prefixes = old['temperatures_keV'], old['densities_atomic_g_cm3'], old['density_prefix_sizes']
    records = {(r['temperature_keV'], r['density_atomic_g_cm3']): r for r in old['records']}
    tables = {Path(name).name: Path(name) for name in old['output_sha256']}
    # Exact temperature knots can round to either side. Restrict node checks
    # to the intersection of both possible four-row interpolation supports.
    points = [(t, r) for it, t in enumerate(tt)
              for r in rr[:min(prefixes[max(0, it - 2):min(len(tt), it + 3)])]]
    node_checks = []
    derivative_checks = []
    for filename, field in (('factor.dat', 'refractive_ratio'),
                             ('normalized_transport.dat', 'native_uncut_times_ratio_atomic_cm2_g'),
                             ('group_transport.dat', 'rosseland_atomic_cm2_g')):
        values = query(a.probe, tables[filename], points, old['X'], a.scratch, filename + '-nodes')
        node_checks.append({'table': filename, 'queries': len(points),
                            'maximum_relative_error': max(abs(v[0] / records[k][field] - 1)
                                for k, v in zip(points, values, strict=True))})
        centers = []
        for it in (0, 5, 10, 15, 20, 25, 30, 34):
            t = math.sqrt(tt[it] * tt[it + 1])
            count = min(prefixes[max(0, it - 1):min(len(tt), it + 3)])
            for ir in sorted({5, 25, 45, 60, count - 3}):
                if ir + 1 < count:
                    centers.append((t, math.sqrt(rr[ir] * rr[ir + 1])))
        queries = []
        h = 2e-5
        for t, r in centers:
            queries.extend([(t, r), (t * math.exp(-h), r), (t * math.exp(h), r),
                            (t, r * math.exp(-h)), (t, r * math.exp(h))])
        values = query(a.probe, tables[filename], queries, old['X'], a.scratch, filename + '-derivatives')
        for i, (t, r) in enumerate(centers):
            v = values[5 * i:5 * i + 5]
            for axis, j, k, d in (('temperature', 1, 2, 1), ('density', 3, 4, 2)):
                finite = (math.log(v[k][0]) - math.log(v[j][0])) / (2 * h)
                error = abs(finite - v[0][d]) / max(1., abs(finite), abs(v[0][d]))
                derivative_checks.append({'table': filename, 'temperature_keV': t,
                    'density_atomic_g_cm3': r, 'axis': axis,
                    'runtime_derivative': v[0][d], 'finite_difference': finite, 'scaled_error': error})
    exclusions = [(tt[0], rr[prefixes[0] - 1] * 1.01), (tt[-1], rr[-1] * 1.01)]
    for i, point in enumerate(exclusions):
        query(a.probe, tables['factor.dat'], [point], old['X'], a.scratch, f'excluded-{i}', reject=True)
    comparisons = []
    if spectral:
        reference = spectral['records']
        points_s = [(r['temperature_keV'], r['density_atomic_g_cm3']) for r in reference]
        values = query(a.probe, tables['normalized_transport.dat'], points_s, old['X'], a.scratch, 'spectral')
        for r, v in zip(reference, values, strict=True):
            key = r['temperature_keV'], r['density_atomic_g_cm3']
            comparisons.append({'temperature_keV': key[0], 'density_atomic_g_cm3': key[1],
                'on_group_table_coordinate': key in records,
                'reference_opacity': r['rosseland_atomic_cm2_g'], 'runtime_opacity': v[0],
                'relative_error': v[0] / r['rosseland_atomic_cm2_g'] - 1})
    compression = None
    if a.compress_base_density:
        if len(rr) != 71 or any(n != 71 for n in prefixes):
            raise ValueError('base compression expects the full 71-density rectangle')
        keep = [0, 20, 30] + list(range(40, 71))
        retained = {(t, rr[i]): records[t, rr[i]] for t in tt for i in keep}
        contents, _ = table_text(retained, tt, [rr[i] for i in keep], old['X'], old['Z'],
                                'refractive_ratio', 'DIMENSIONLESS refractive ratio; omitted-node comparison')
        path = a.scratch / 'compressed-factor.dat'
        path.write_text(contents)
        omitted = [(t, r) for t in tt for r in rr if (t, r) not in retained]
        values = query(a.probe, path, omitted, old['X'], a.scratch, 'omitted-density')
        compression = {'retained_density_indices': keep, 'source_density_nodes': len(rr),
                       'retained_density_nodes': len(keep), 'omitted_source_comparisons': len(omitted),
                       'fresh_independent_sources': False,
                       'maximum_relative_error': max(abs(v[0] / records[k]['refractive_ratio'] - 1)
                           for k, v in zip(omitted, values, strict=True))}
    verify(inputs)
    node_max = max(c['maximum_relative_error'] for c in node_checks)
    derivative_max = max(c['scaled_error'] for c in derivative_checks)
    spectral_max = max((abs(c['relative_error']) for c in comparisons), default=None)
    result = {'scope': __doc__, 'accepted_for_stellar_opacity': False, 'X': old['X'], 'Z': old['Z'],
              'source_states': old['states'], 'source_coordinates_tested_per_table': len(points),
              'source_coordinates_omitted_for_derivative_support': old['states'] - len(points),
              'source_node_checks': node_checks, 'derivative_checks': derivative_checks,
              'maximum_scaled_derivative_error': derivative_max, 'excluded_queries_rejected': len(exclusions),
              'runtime_numerical_checks_passed': node_max < 1e-12 and derivative_max < 1e-5,
              'source_uncut_recovery_check_passed': old['uncut_recovery_check_passed'],
              'spectral_comparisons': comparisons, 'maximum_spectral_relative_error': spectral_max,
              'spectral_relative_criterion': .005,
              'sampled_spectral_comparisons_passed': None if spectral_max is None else spectral_max <= .005,
              'density_compression': compression, 'input_sha256': inputs,
              'scratch_sha256': {p.name: digest(p) for p in a.scratch.iterdir() if p.is_file()}}
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in (
        'input_sha256', 'scratch_sha256', 'derivative_checks', 'spectral_comparisons')}), flush=True)
    if not result['runtime_numerical_checks_passed'] or result['sampled_spectral_comparisons_passed'] is False:
        raise ValueError('candidate table check failed')


if __name__ == '__main__':
    main()
