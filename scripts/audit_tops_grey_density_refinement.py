#!/usr/bin/env python3
"""Test independently sourced density refinement before expanding all mixtures.

Compare uncut grey means and grey-times-refractive-ratio tables through Ember's
actual table interpolator. Reference density points remain withheld. This is
a three-composition pilot, not acceptance of a complete stellar opacity family.
"""
import argparse
import gzip
import json
import math
from pathlib import Path

from assemble_tops_refractive_family import ratio_corners
from audit_refractive_opacity_family_v2 import stable_change
from audit_tops_electron_dispersion import KEV, KB
from audit_tops_group_factor_tables import query
from fetch_tops_composition import digest
from import_tops_composition import read
from reduce_tops_group_factors import add_inputs, verify, table_text


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base-plan', type=Path, required=True)
    p.add_argument('--supplement-plan', type=Path, action='append', required=True)
    p.add_argument('--scratch', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    a = p.parse_args()
    if a.scratch.exists() or a.report.exists():
        raise FileExistsError('retain previous pilot outputs')
    inputs = {}

    def pin(path):
        path = Path(path)
        add_inputs(inputs, {str(path.resolve()): digest(path)})
        return path

    runtime_path = pin('docs/results/tops_refractive_hot_family_runtime_v2.json')
    runtime = json.loads(runtime_path.read_text())
    add_inputs(inputs, runtime['input_sha256'])
    archive = pin('/tmp/ember-refractive-family-runtime-v1/independent-groups.json.gz')
    if runtime['input_sha256'].get(str(archive.resolve())) != digest(archive):
        raise ValueError('unidentified prior conduction query archive')
    old = json.loads(gzip.decompress(archive.read_bytes()))
    if old['returncode'] != 0:
        raise ValueError('failed prior conduction query')
    old_rows = [json.loads(line) for line in old['stdout'].splitlines()]
    requests = [list(map(float, line.split())) for line in old['request'].splitlines()]
    if len(requests) != len(old_rows) or any(q != r['query'] for q, r in zip(requests, old_rows)):
        raise ValueError('prior conduction queries are misaligned')
    conduction = {tuple(r['query']): r for r in old_rows}
    probe = Path('/tmp/ember-tops-spectral-runtime-probe-v1')
    if inputs.get(str(probe.resolve())) != digest(probe):
        raise ValueError('unidentified actual table probe')
    for name in ('audit_tops_grey_density_refinement.py', 'assemble_tops_refractive_family.py',
                 'audit_refractive_opacity_family_v2.py', 'audit_tops_group_factor_tables.py',
                 'reduce_tops_group_factors.py', 'import_tops_composition.py'):
        pin(Path(__file__).with_name(name))
    ratios = {}
    for x in ('000', '020', '050', '075'):
        for z in ('010', '030'):
            path = pin(f'docs/results/tops_ratio_family_x{x}-z{z}_reduction_v1.json')
            r = json.loads(path.read_text())
            add_inputs(inputs, r['input_sha256']); add_inputs(inputs, r['output_sha256'])
            ratios[r['X'], r['Z']] = r
    tt = next(iter(ratios.values()))['temperatures_keV']
    refs = {}
    for x in ('010', '030', '070'):
        path = pin(f'docs/results/tops_ratio_hydrogen_control_x{x}-z020_reduction_v1.json')
        r = json.loads(path.read_text())
        add_inputs(inputs, r['input_sha256']); add_inputs(inputs, r['output_sha256'])
        refs[r['X'], r['Z']] = r
    plans = [json.loads(pin(path).read_text()) for path in [a.base_plan, *a.supplement_plan]]
    if any(plan['baseline_manifest_sha256'] != digest(Path(plan['baseline_manifest'])) for plan in plans):
        raise ValueError('baseline manifest changed')
    baseline = json.loads(pin(plans[0]['baseline_manifest']).read_text())
    cells, added, sources = {}, {}, []
    for plan_index, plan in enumerate(plans):
        for j in plan['requests']:
            key = j['X'], j['Z']
            if key not in refs or max(j['densities_g_cm3']) > 10000 or not any(t in tt for t in j['temperatures_keV']):
                continue
            root = Path(j['work'])
            receipt = json.loads(pin(root / 'receipt.json').read_text())
            request = json.loads(pin(root / 'request.json').read_text())
            recipe = json.loads(pin(root / 'recipe.json').read_text())
            if ((receipt['X'], receipt['Z']) != key or recipe['job'] != j or
                    digest(root / 'request.json') != receipt['request_sha256'] or
                    request['datype'] != 'gray' or request['plasnu'] != 'off' or request['lib'] != 'new'):
                raise ValueError('wrong source recipe, composition or opacity convention')
            native_t, native_r, part, excluded = read(pin(root / 'source.txt'), receipt,
                dimensions=(len(j['temperatures_keV']), len(j['densities_g_cm3'])))
            if native_t != j['temperatures_keV'] or any(abs(r/q-1) > 5e-5 for r, q in zip(native_r, j['densities_g_cm3'], strict=True)):
                raise ValueError('returned source axes differ from requested axes')
            add_inputs(inputs, recipe.get('input_sha256', {}))
            dest = cells.setdefault(key, {})
            for point, value in part.items():
                if point[0] not in tt:
                    continue
                if point in excluded:
                    raise ValueError('pilot base domain contains a source exclusion')
                if point in dest:
                    raise ValueError('duplicate source coordinate in pilot')
                dest[point] = {'grey': value}
                if plan_index:
                    added.setdefault(key, set()).add(point)
            sources.append({'source': str(root), 'states': len(part), 'supplement': bool(plan_index)})
    for key, ref in refs.items():
        if any((r['temperature_keV'], r['density_atomic_g_cm3']) in added[key] for r in ref['records']):
            raise ValueError('an independent reference point was inserted in the refined grid')
        old_plane = next(r for r in baseline['planes'] if (r['X'], r['Z']) == key)
        path = pin(Path(plans[0]['baseline_manifest']).parent / old_plane['file'])
        bt, br, _, _ = read(path, old_plane)
        actual = sorted({r for t, r in cells[key]})
        if [t for t in bt if t >= .025] != tt or not set(br).issubset(actual):
            raise ValueError('original native grid was not retained')
        if set(cells[key]) != {(t, r) for t in tt for r in actual}:
            raise ValueError('incomplete merged pilot rectangle')
    verify(inputs)
    a.scratch.mkdir(parents=True)
    all_points = sorted(set().union(*(set(c) for c in cells.values())))
    factor_values = {}
    for key, r in ratios.items():
        table = next(Path(name) for name in r['output_sha256'] if Path(name).name == 'factor.dat')
        values = query(probe, table, all_points, key[0], a.scratch, f'ratio-x{key[0]:g}-z{key[1]:g}')
        factor_values[key] = {point: row[0] for point, row in zip(all_points, values, strict=True)}
    xs, zs = sorted({x for x, z in ratios}), sorted({z for x, z in ratios})
    comparisons, omissions, node_checks, planes = [], [], [], []
    for key, records in sorted(cells.items()):
        x, z = key
        corners = ratio_corners(x, z, xs, zs)
        for point, row in records.items():
            factor = math.exp(sum(w*math.log(factor_values[c][point]) for c, w in corners))
            if factor < 1-1e-12:
                raise ValueError('invalid refractive ratio')
            row['normalized'] = row['grey'] * factor
        rr = sorted({r for t, r in records})
        ref = refs[key]
        selected = [r for r in ref['records'] if r['density_atomic_g_cm3'] <= rr[-1]]
        omissions += [{'X': x, 'Z': z, 'temperature_keV': r['temperature_keV'],
                       'density_atomic_g_cm3': r['density_atomic_g_cm3'],
                       'reason': 'dense extension outside this base-density pilot'}
                      for r in ref['records'] if r['density_atomic_g_cm3'] > rr[-1]]
        points = [(r['temperature_keV'], r['density_atomic_g_cm3']) for r in selected]
        queried = {}
        for field in ('grey', 'normalized'):
            table = a.scratch / f'x{x:g}-z{z:g}-{field}.dat'
            contents, prefixes = table_text(records, tt, rr, x, z, field, 'Density refinement pilot; not selected stellar opacity')
            table.write_text(contents)
            queried[field] = query(probe, table, points, x, a.scratch, f'x{x:g}-{field}-independent')
            nodes = sorted(records)
            values = query(probe, table, nodes, x, a.scratch, f'x{x:g}-{field}-nodes')
            error = max(abs(v[0]/records[q][field]-1) for q, v in zip(nodes, values, strict=True))
            node_checks.append({'X': x, 'field': field, 'queries': len(nodes), 'maximum_relative_error': error})
        planes.append({'X': x, 'Z': z, 'temperatures': len(tt), 'densities': len(rr), 'added_source_states': len(added[key])})
        for i, ref in enumerate(selected):
            t, rho = points[i]
            direct = ref['native_uncut_times_ratio_atomic_cm2_g']
            delta = queried['normalized'][i][0]/direct-1
            prior = conduction[x, z, 0., t*KEV/KB, rho]
            change = stable_change(delta, direct, prior['conductive_opacity']) if prior.get('conduction_covered') else None
            comparisons.append({'X': x, 'Z': z, 'temperature_keV': t, 'density_atomic_g_cm3': rho,
                'grey_reference': ref['uncut_source_rosseland'], 'grey_interpolated': queried['grey'][i][0],
                'grey_relative_error': queried['grey'][i][0]/ref['uncut_source_rosseland']-1,
                'normalized_reference': direct, 'normalized_interpolated': queried['normalized'][i][0],
                'radiative_relative_error': delta, 'combined_relative_error': change,
                'combined_absolute_error_bound': abs(delta) if change is None else abs(change),
                'off_merged_density_grid': rho not in rr, 'source': ref['source']})
    high = [r for r in comparisons if math.log10(r['temperature_keV']*KEV/KB) >= 5.7]
    maximum = max(r['combined_absolute_error_bound'] for r in high)
    result = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
        'planes': planes, 'sources': sources, 'comparisons': comparisons, 'omissions': omissions,
        'node_checks': node_checks, 'maximum_node_relative_error': max(r['maximum_relative_error'] for r in node_checks),
        'maximum_grey_relative_error': max(abs(r['grey_relative_error']) for r in comparisons),
        'maximum_radiative_relative_error': max(abs(r['radiative_relative_error']) for r in comparisons),
        'fully_selected_high_comparisons': len(high), 'high_combined_absolute_error_upper_bound': maximum,
        'relative_comparison_criterion': .005, 'high_combined_comparison_passed': maximum <= .005,
        'high_failing_comparisons': sum(r['combined_absolute_error_bound'] > .005 for r in high),
        'independent_reference_points_added_to_grid': False,
        'limitations': ['Three compositions at native temperatures only; this does not test the full X/Z family or temperature blend.',
            'Conduction is reused from identified production-class queries at exactly the same material states.',
            'References test interpolation against the same underlying source library, not the absolute accuracy of its opacity physics.'],
        'input_sha256': inputs, 'output_sha256': {str(p.resolve()): digest(p) for p in a.scratch.iterdir() if p.is_file()}}
    verify(inputs)
    a.report.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('scope', 'sources', 'comparisons', 'omissions', 'input_sha256', 'output_sha256', 'limitations')}), flush=True)
    if result['maximum_node_relative_error'] > 1e-10 or maximum > .005:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
