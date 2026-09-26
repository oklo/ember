#!/usr/bin/env python3
"""Assemble the hot native-uncut opacity family times a refractive ratio.

Retain the fine grey composition/density grids. Evaluate each needed ratio
coordinate once through Ember, then interpolate log ratio in composition.
Missing sources, exclusions, and conservative ratio support remain explicit.
The existing lower-temperature tables are copied unchanged; their connection
to the corrected hot tables requires a separate physical/runtime comparison.
No generated family is accepted automatically for stellar evolution.
"""
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import shlex
import shutil

from audit_tops_electron_dispersion import KEV, KB
from audit_tops_factor_hydrogen import bracket
from audit_tops_group_factor_tables import query
from fetch_tops_composition import digest
from import_tops_composition import read
from reduce_tops_group_factors import add_inputs, verify


def ratio_corners(x, z, xs, zs):
    xl, xh, fx = bracket(xs, x)
    zl, zh, fz = bracket(zs, z)
    return [((xx, zz), wx*wz) for xx, wx in ((xl, 1-fx), (xh, fx))
            for zz, wz in ((zl, 1-fz), (zh, fz)) if wx*wz > 0]


def ratio_limit(table, t):
    tt = table['temperatures_keV']
    i = tt.index(t)
    n = min(table['density_prefix_sizes'][max(0, i-2):min(len(tt), i+3)])
    return table['densities_atomic_g_cm3'][n-1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('grey_plan', 'lower_temperature_family', 'probe', 'probe_reference',
                 'hydrogen_comparison', 'scratch', 'output', 'report'):
        p.add_argument('--' + name.replace('_', '-'), type=Path, required=True)
    p.add_argument('--ratio', type=Path, nargs='+', required=True)
    a = p.parse_args()
    if any(path.exists() for path in (a.scratch, a.output, a.report)):
        raise FileExistsError('preserve existing candidate outputs')
    inputs = {}
    plan = json.loads(a.grey_plan.read_text())
    baseline = Path(plan['baseline_manifest'])
    if digest(baseline) != plan['baseline_manifest_sha256']:
        raise ValueError('grey baseline manifest changed')
    expected = json.loads(baseline.read_text())['planes']
    ratios = {}
    for path in a.ratio:
        r = json.loads(path.read_text())
        key = r['X'], r['Z']
        if key in ratios:
            raise ValueError('duplicate ratio composition')
        ratios[key] = r
        add_inputs(inputs, r['input_sha256'])
        add_inputs(inputs, r['output_sha256'])
        add_inputs(inputs, {str(path.resolve()): digest(path)})
    xs = sorted({x for x, z in ratios})
    zs = sorted({z for x, z in ratios})
    if len(xs) < 2 or len(zs) < 2 or set(ratios) != {(x, z) for x in xs for z in zs}:
        raise ValueError('incomplete ratio composition rectangle')
    build = json.loads(a.probe_reference.read_text())
    check = json.loads(a.hydrogen_comparison.read_text())
    if not check['high_domain_sampled_ratio_transport_check_passed']:
        raise ValueError('independent high-domain ratio comparison failed')
    for r in (build, check):
        add_inputs(inputs, r['input_sha256'])
    for path in (a.grey_plan, baseline, a.probe_reference, a.hydrogen_comparison,
                 Path(__file__), Path('scripts/audit_tops_factor_hydrogen.py'),
                 Path('scripts/audit_tops_group_factor_tables.py')):
        add_inputs(inputs, {str(path.resolve()): digest(path)})
    if inputs.get(str(a.probe.resolve())) != digest(a.probe):
        raise ValueError('unidentified runtime executable')
    verify(inputs)
    tt = next(iter(ratios.values()))['temperatures_keV']
    if any(r['temperatures_keV'] != tt for r in ratios.values()) or tt[0] != .025:
        raise ValueError('unexpected hot ratio temperature axes')
    compositions = {(v['X'], v['Z']): v for v in expected}
    if len(compositions) != len(expected):
        raise ValueError('duplicate grey composition in baseline')
    required_z = sorted({z for x, z in compositions})
    required_x = sorted({x for x, z in compositions})
    if set(compositions) != {(x, z) for x in required_x for z in required_z}:
        raise ValueError('incomplete fine grey composition rectangle')
    corners = {key: ratio_corners(*key, xs, zs) for key in compositions}
    limits = {(x, z, t): min(ratio_limit(ratios[key], t) for key, weight in corners[x, z])
              for x, z in compositions for t in tt}
    jobs = list(plan['requests'])
    for retained in plan.get('retained_requests', []):
        path = Path(retained['plan'])
        if digest(path) != retained['plan_sha256']:
            raise ValueError('retained source plan changed')
        add_inputs(inputs, {str(path.resolve()): digest(path)})
        jobs.append(json.loads(path.read_text())['requests'][retained['request_index']])
    selected = []
    outside = []
    for j in jobs:
        x, z = j['X'], j['Z']
        if (x, z) not in compositions:
            raise ValueError('unexpected grey source composition')
        dense = max(j['densities_g_cm3']) > 10000
        use = any(t in tt and (not dense or limits[x, z, t] > 10000)
                  for t in j['temperatures_keV'])
        (selected if use else outside).append(j)
    missing = [j['work'] for j in selected if not (Path(j['work']) / 'receipt.json').is_file()]
    if missing:
        raise ValueError(f'{len(missing)} required source batches are not complete; first: {missing[:4]}')
    cells, exclusions, sources = {}, {}, []
    overlaps = 0
    for j in selected:
        source = Path(j['work'])
        receipt = json.loads((source / 'receipt.json').read_text())
        request = json.loads((source / 'request.json').read_text())
        recipe = json.loads((source / 'recipe.json').read_text())
        x, z = j['X'], j['Z']
        if ((receipt['X'], receipt['Z']) != (x, z) or recipe['job'] != j or
                digest(source / 'request.json') != receipt['request_sha256'] or
                request['datype'] != 'gray' or request['plasnu'] != 'off' or request['lib'] != 'new'):
            raise ValueError('wrong grey scientific source or request')
        native_t, native_r, part, excluded = read(source / 'source.txt', receipt,
            dimensions=(len(j['temperatures_keV']), len(j['densities_g_cm3'])))
        if native_t != j['temperatures_keV'] or any(abs(r/q-1) > 5e-5 for r, q in
                                                   zip(native_r, j['densities_g_cm3'], strict=True)):
            raise ValueError('grey source coordinates differ from request')
        used = 0
        for (t, rho), value in part.items():
            if t not in tt or rho > limits[x, z, t]:
                continue
            key = x, z, t, rho
            flag = (t, rho) in excluded
            if key in cells:
                if exclusions[key] != flag or (not flag and cells[key] != value):
                    raise ValueError('inconsistent overlapping source state')
                overlaps += 1
            else:
                cells[key], exclusions[key] = value, flag
            used += 1
        for name in ('receipt.json', 'request.json', 'recipe.json', 'source.txt'):
            path = source / name
            add_inputs(inputs, {str(path.resolve()): digest(path)})
        add_inputs(inputs, recipe.get('input_sha256', {}))
        sources.append({'work': str(source), 'returned_states': len(part),
                        'source_exclusions': len(excluded), 'states_in_requested_hot_domain': used})
    # The fine native base grid and the ratio's fine dense grid give an
    # independent requirement for the assembled density axes.
    base_axes = []
    for key, old in compositions.items():
        path = baseline.parent / old['file']
        native_t, native_r, _, _ = read(path, old)
        if [t for t in native_t if t >= .025] != tt:
            raise ValueError('hot temperatures differ from the native source grid')
        base_axes.append(native_r)
        add_inputs(inputs, {str(path.resolve()): digest(path)})
    if any(rr != base_axes[0] for rr in base_axes):
        raise ValueError('fine native density axes differ')
    dense_axes = [tuple(rho for rho in r['densities_atomic_g_cm3'] if rho > 10000)
                  for r in ratios.values() if r['densities_atomic_g_cm3'][-1] > 10000]
    if not dense_axes or len(set(dense_axes)) != 1:
        raise ValueError('missing or inconsistent dense ratio axes')
    rr = sorted(set(base_axes[0]) | set(dense_axes[0]))
    requests = defaultdict(set)
    retained_rows = {}
    omitted_cells = []
    for x, z in sorted(compositions):
        for t in tt:
            prefix = []
            ended = False
            for rho in rr:
                key = x, z, t, rho
                planned = rho <= limits[x, z, t]
                if planned and key not in cells:
                    raise ValueError(f'missing required grey source state {key}')
                excluded = planned and exclusions[key]
                if not planned or excluded:
                    ended = True
                    if excluded:
                        omitted_cells.append({'X': x, 'Z': z, 'temperature_keV': t,
                                              'density_atomic_g_cm3': rho, 'reason': 'source exclusion'})
                else:
                    if ended:
                        raise ValueError('valid grey source values do not form a density prefix')
                    prefix.append(rho)
                    for corner, weight in corners[x, z]:
                        requests[corner].add((t, rho))
            if len(prefix) < 4:
                raise ValueError('insufficient supported density prefix')
            retained_rows[x, z, t] = prefix
    verify(inputs)
    a.scratch.mkdir(parents=True)
    ratios_at = {}
    for (x, z), points in sorted(requests.items()):
        points = sorted(points)
        r = ratios[x, z]
        table = next(Path(path) for path in r['output_sha256'] if Path(path).name == 'factor.dat')
        values = query(a.probe, table, points, x, a.scratch, f'x{x:g}-z{z:g}')
        if any(row[0] < 1-1e-12 for row in values):
            raise ValueError('interpolated refractive ratio is below unity')
        ratios_at[x, z] = {key: row[0] for key, row in zip(points, values, strict=True)}
    a.output.mkdir(parents=True)
    copied, relocated = [], []
    for name in ('aesopus21_gs98_mixture.dat', 'tops_gs98_mixture_low.dat'):
        manifest = a.lower_temperature_family / name
        lines = manifest.read_text().splitlines()
        fields = lines[0].split()
        if fields[:2] != ['EMBER_OPACITY_MIXTURE', '1'] or int(fields[2]) != len(lines)-1:
            raise ValueError('invalid lower-temperature family manifest')
        add_inputs(inputs, {str(manifest.resolve()): digest(manifest)})
        rewritten = [lines[0]]
        for line in lines[1:]:
            z, filename = shlex.split(line)
            path = (manifest.parent / filename).resolve(strict=True)
            if not path.is_relative_to(a.lower_temperature_family.resolve().parent):
                raise ValueError('lower-temperature source escapes its opacity directory')
            target = a.output / path.name
            if target.exists():
                if digest(target) != digest(path):
                    raise ValueError('conflicting lower-temperature table names')
            else:
                shutil.copyfile(path, target)
            add_inputs(inputs, {str(path.resolve()): digest(path)})
            if digest(target) != digest(path):
                raise ValueError('lower-temperature copy changed')
            copied.append(str(target.resolve()))
            rewritten.append(f'{z} "{path.name}"')
        (a.output / name).write_text('\n'.join(rewritten) + '\n')
        relocated.append({'source': str(manifest.resolve()),
                          'output': str((a.output / name).resolve()),
                          'change': 'Only referenced table paths are relocated; header and metallicities are retained.'})
    high_manifest = [f'EMBER_OPACITY_MIXTURE 1 {len(required_z)} logRho ATOMIC_native_uncut_times_refractive_ratio']
    planes = []
    for z in required_z:
        lines = [f'EMBER_OPACITY_TABLE 2 {len(required_x)} {len(tt)} {len(rr)} Native UNCUT grey times refractive ratio; conservative support',
                 ' '.join(format(math.log10(rho), '.17g') for rho in rr),
                 ' '.join(format(math.log10(t*KEV/KB), '.17g') for t in tt)]
        for x in required_x:
            lines.append(f'{x:.17g} {z:.17g}')
            prefixes = []
            for t in tt:
                prefix = retained_rows[x, z, t]
                values = []
                for rho in prefix:
                    logf = sum(weight*math.log(ratios_at[corner][t, rho])
                               for corner, weight in corners[x, z])
                    logk = math.log(cells[x, z, t, rho]) + logf
                    if not math.isfinite(logk):
                        raise ValueError('nonfinite corrected opacity')
                    values.append(logk / math.log(10))
                prefixes.append(len(prefix))
                lines.append(str(len(prefix)) + ' ' + ' '.join(format(v, '.17g') for v in values))
            planes.append({'X': x, 'Z': z, 'density_prefix_sizes': prefixes,
                           'states': sum(prefixes)})
        name = f'tops_gs98_mixture_z{round(z*1000):03d}_high.dat'
        (a.output / name).write_text('\n'.join(lines) + '\n')
        high_manifest.append(f'{z:.17g} "{name}"')
    (a.output / 'tops_gs98_mixture_high.dat').write_text('\n'.join(high_manifest) + '\n')
    verify(inputs)
    result = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
              'complete_native_composition_rectangle': True, 'X': required_x, 'Z': required_z,
              'temperature_keV': tt, 'density_atomic_g_cm3': rr,
              'hot_states_written': sum(r['states'] for r in planes),
              'grey_source_batches_used': len(selected),
              'grey_source_batches_outside_requested_hot_support': len(outside),
              'exact_overlapping_source_states': overlaps,
              'unique_ratio_runtime_queries': sum(len(points) for points in requests.values()),
              'ratio_evaluations_if_repeated_for_every_grey_state': sum(
                  len(corners[x, z]) * len(prefix) for (x, z, t), prefix in retained_rows.items()),
              'copied_lower_temperature_files': copied, 'planes': planes,
              'relocated_lower_temperature_manifests': relocated,
              'source_exclusions': omitted_cells, 'sources': sources,
              'limitations': ['Lower-temperature opacities are unchanged and their physical approximation must be assessed separately.',
                              'The strict global radiative-ratio and group/native recovery failures remain explicit.',
                              'The final assembled interpolation, temperature blend, isotope mapping and stellar response are not checked by assembly.',
                              'Cooler dense states outside the conservative ratio support are not filled or extrapolated.'],
              'input_sha256': inputs,
              'output_sha256': {str(path.resolve()): digest(path) for path in a.output.rglob('*') if path.is_file()},
              'scratch_sha256': {path.name: digest(path) for path in a.scratch.iterdir() if path.is_file()}}
    a.report.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key not in
                     ('scope', 'X', 'Z', 'temperature_keV', 'density_atomic_g_cm3', 'planes',
                      'sources', 'source_exclusions', 'limitations', 'input_sha256', 'output_sha256',
                      'scratch_sha256')}), flush=True)


if __name__ == '__main__':
    main()
