#!/usr/bin/env python3
"""Extend a checked mixed H/He atmosphere table from completed source columns.

The plan names a parent directory, new columns and independent check columns.
Existing source values and interpolation intervals are preserved. The output
remains a candidate until its native interpolation has also been checked.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import shlex

from assemble_nongrey_grid import physical_identity
from generate_nongrey_grid import sequence, temperatures
from import_nongrey_grid import (read_text, source_diagnostics_match,
                                source_inputs, source_state)
from prepare_nongrey_sources import digest


def read_parent(directory):
    manifest = json.loads((directory / 'assembly.json').read_text())
    table = directory / 'mixed_tau100.dat'
    if digest(table) != manifest['table_sha256']:
        raise ValueError('parent table differs from its recorded source values')
    lines = table.read_text().splitlines()
    if lines[0] != 'EMBER_COMPOSITION_ATMOSPHERE 4':
        raise ValueError('expected a table with explicit interpolation cells')
    end = lines.index('data') + 1
    header = dict((v[0], v[1:]) for v in map(shlex.split, lines[1:end-1]))
    if header['helium3'] != ['1', '0']:
        raise ValueError('this extension uses zero source He3; isotope mapping is in the caller')
    axes = manifest['axes']
    for i, name in enumerate(['hydrogen', 'log_teff', 'log_g']):
        actual = list(map(float, header[name][1:]))
        expected = [math.log10(x) if i == 1 else x for x in axes[i]]
        if actual != expected or int(header[name][0]) != len(actual):
            raise ValueError('parent source and stored axes disagree')
    count = math.prod(map(len, axes))
    rows = dict(zip(itertools.product(*axes), lines[end:end+count], strict=True))
    source_rows = {tuple(v['coordinates']): v['diagnostics']
                   for v in manifest['source_columns']}
    for key, line in rows.items():
        fields = line.split()
        if fields == ['0']:
            if key in source_rows:
                raise ValueError('parent source manifest contradicts a missing row')
            continue
        if len(fields) != 3 or fields[0] != '1' or key not in source_rows:
            raise ValueError('parent table and source columns disagree')
        actual = list(map(float, fields[1:]))
        expected = [math.log10(source_rows[key][k]) for k in ['T', 'Pgas']]
        if any(not math.isfinite(v) or abs(v - w) > 1e-12
               for v, w in zip(actual, expected)):
            raise ValueError('parent matching state differs from its source record')
    if lines[end+count] != 'cells':
        raise ValueError('missing parent cell flags')
    keys = list(cell_keys(axes))
    flags = lines[end+count+1:]
    if len(flags) != len(keys) or any(v not in ['0', '1'] for v in flags):
        raise ValueError('invalid parent cell flags')
    return manifest, lines[:end], rows, dict(zip(keys, map(int, flags)))


def cell_keys(axes):
    return itertools.product(*(zip(a[:-1], a[1:]) for a in axes))


def merged_axes(old, coordinates):
    axes = []
    for i, axis in enumerate(old):
        new = sorted(set(axis) | {k[i] for k in coordinates})
        if any(not math.isfinite(v) for v in new):
            raise ValueError('nonfinite source coordinate')
        if any(axis[0] < v < axis[-1] and v not in axis for v in new):
            raise ValueError('inserting an interior coordinate would change old interpolation')
        axes.append(new)
    return axes


def render(header, old_axes, old_rows, old_cells, additions):
    rows = dict(old_rows)
    for key, state in additions.items():
        if key in rows and rows[key] != '0':
            raise ValueError('cannot replace an existing source column')
        if not all(math.isfinite(state[k]) and state[k] > 0 for k in ['T', 'Pgas']):
            raise ValueError('invalid matching temperature or pressure')
        rows[key] = '1 ' + ' '.join(format(math.log10(state[k]), '.17g')
                                   for k in ['T', 'Pgas'])
    axes = merged_axes(old_axes, additions)
    output = list(header)
    for name, axis, previous in zip(['hydrogen', 'log_teff', 'log_g'], axes, old_axes):
        if axis == previous:
            continue
        values = [math.log10(v) if name == 'log_teff' else v for v in axis]
        i = next(i for i, line in enumerate(output) if line.startswith(name + ' '))
        output[i] = f'{name} {len(values)} ' + ' '.join(format(v, '.17g') for v in values)
    output += [rows.get(k, '0') for k in itertools.product(*axes)] + ['cells']
    cells = []
    for key in cell_keys(axes):
        corners = list(itertools.product(*key))
        present = all(rows.get(k, '0') != '0' for k in corners)
        # An intentionally disabled cell stays disabled. A cell missing only
        # source columns can become usable when those columns are supplied.
        if key in old_cells and not old_cells[key]:
            if all(old_rows.get(k, '0') != '0' for k in corners):
                present = False
        if old_cells.get(key) and not present:
            raise ValueError('an existing interpolation cell was lost')
        output.append(str(int(present)))
        cells.append(dict(hydrogen=key[0], teff=key[1], log_g=key[2], supported=present))
    return '\n'.join(output) + '\n', axes, cells


def load_column(directory, expected_physics):
    def archived(name):
        path = directory / name
        return read_text(path if path.exists() else path.with_name(name + '.gz'))

    result_path = directory / 'result.json'
    result = json.loads(result_path.read_text())
    if not result['accepted'] or not result.get('control_pass', True):
        raise ValueError(f'unaccepted atmosphere: {directory}')
    for name, checksum in result['retained_sha256'].items():
        if digest(directory / name) != checksum:
            raise ValueError(f'changed source output: {directory / name}')
    spec = json.loads((directory / 'specification.json').read_text())
    prepared = json.loads((directory / 'provenance.json').read_text())
    if physical_identity(spec, prepared) != expected_physics:
        raise ValueError(f'atmosphere physics differs: {directory}')
    if spec['helium3'] != [0.] or any(len(spec[k]) != 1
                                    for k in ['hydrogen', 'teff_K', 'log_g']):
        raise ValueError('one zero-He3 source column is required per directory')
    x, t, g = (spec[k][0] for k in ['hydrogen', 'teff_K', 'log_g'])
    opacity = directory / 'opacity.bin'
    if digest(opacity) != (directory / 'opacity.sha256').read_text().strip():
        raise ValueError('source opacity differs from the completed column')
    log = archived('run.log')
    inputs = {k: (directory / name).read_text() for k, name in
              [('atmosphere_input', 'fort.5'), ('parameters', 'tas'),
               ('element_masses', 'ember-masses.dat'), ('initial_structure', 'fort.8')]}
    source_inputs(inputs, spec, x, 0., t, g, log)
    state = source_state(log, archived('fort.9'), t, g,
                         {'temperature_K': temperatures(spec),
                          'density_g_cm3': sequence(spec['log_density'])}, spec['tau'])
    if not source_diagnostics_match(result['diagnostics'], state):
        raise ValueError('source diagnostics differ from the completion record')
    return (x, t, g), state, dict(coordinates=[x, t, g], directory=str(directory),
                                 result_sha256=digest(result_path),
                                 opacity_sha256=digest(opacity), diagnostics=state)


def interpolate(axes, nodes, key):
    intervals, weights = [], []
    for i, (a, x) in enumerate(zip(axes, key)):
        j = next((j for j in range(len(a)-1) if a[j] <= x <= a[j+1]), None)
        if j is None:
            raise ValueError('independent check lies outside the assembled axes')
        lo, hi = a[j:j+2]
        w = math.log(x/lo)/math.log(hi/lo) if i == 1 else (x-lo)/(hi-lo)
        intervals.append([lo, hi]); weights.append(w)
    value = {}
    for field in ['T', 'Pgas']:
        total = 0.
        for corner in itertools.product([0, 1], repeat=3):
            q = tuple(a[j] for a, j in zip(intervals, corner))
            if q not in nodes:
                raise ValueError('independent check needs a missing source column')
            w = math.prod(u if j else 1-u for u, j in zip(weights, corner))
            total += w * math.log(nodes[q][field])
        value[field] = math.exp(total)
    return value


def assemble(plan_path, output):
    plan_path = plan_path.resolve(strict=True)
    plan = json.loads(plan_path.read_text())
    resolve = lambda p: (plan_path.parent / p).resolve(strict=True)
    parent = resolve(plan['parent'])
    base, header, old_rows, old_cells = read_parent(parent)
    if (output / 'mixed_tau100.dat').exists() or (output / 'assembly.json').exists():
        raise FileExistsError('preserve the existing output table')
    columns = {tuple(v['coordinates']): v for v in base['source_columns']}
    physics = base.get('physical_identity')
    if physics is None:
        reference = Path(next(iter(columns.values()))['directory'])
        physics = physical_identity(json.loads((reference/'specification.json').read_text()),
                                    json.loads((reference/'provenance.json').read_text()))
    added = {}
    for path in plan['sources']:
        key, state, record = load_column(resolve(path), physics)
        if key in added:
            raise ValueError('duplicate new source column')
        added[key] = state
        columns[key] = record
    if set(added) != {tuple(k) for k in plan['required_coordinates']}:
        raise ValueError('new columns differ from the requested extension')
    text, axes, cells = render(header, base['axes'], old_rows, old_cells, added)
    nodes = {k: v['diagnostics'] for k, v in columns.items()}
    checks = []
    for path in plan['independent_checks']:
        key, state, record = load_column(resolve(path), physics)
        if key in nodes:
            raise ValueError('an independent check duplicates a table column')
        if not any(c['supported'] and all(lo <= x <= hi for x, (lo, hi) in
                   zip(key, [c['hydrogen'], c['teff'], c['log_g']])) for c in cells):
            raise ValueError('independent check lies in an unsupported cell')
        value = interpolate(axes, nodes, key)
        checks.append(dict(hydrogen=key[0], teff_K=key[1], log_g=key[2],
                           source=record['directory'],
                           relative_error={k: value[k]/state[k]-1 for k in value}))
    if not checks:
        raise ValueError('include at least one independent interpolation check')
    output.mkdir(parents=True, exist_ok=True)
    table = output / 'mixed_tau100.dat'
    table.write_text(text)
    report = dict(axes=axes, source_columns=list(columns.values()),
                  supported_cells=cells, held_out=checks, table_sha256=digest(table),
                  selected_for_evolution=False, matching_tau=base['matching_tau'],
                  limitations=base['limitations'], physical_identity=physics,
                  parent=str(parent),
                  parent_table_sha256=base['table_sha256'],
                  preserved_valid_vertices=sum(v != '0' for v in old_rows.values()),
                  new_source_columns=len(added), plan_sha256=digest(plan_path),
                  assembler_sha256=digest(__file__))
    (output / 'assembly.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('plan', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    result = assemble(args.plan, args.output)
    print(json.dumps({k: result[k] for k in
                     ['table_sha256', 'new_source_columns', 'held_out']}, indent=2))
