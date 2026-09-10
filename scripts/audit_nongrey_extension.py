#!/usr/bin/env python3
"""Check an assembled gas extension with Ember's runtime interpolator and EOS.

Reassembles pinned sources, checks retained boundary values, reports source/EOS
density differences, and compares independent canonical heldout atmospheres.
Local discrepancies are diagnostics, not global physical error bounds.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

from assemble_nongrey_grid import load_continuation, physical_identity
from import_nongrey_grid import import_grid
from prepare_nongrey_sources import digest


def runtime(probe, eos, table, points):
    inputs = ''.join(' '.join(format(v, '.17g') for v in p)+'\n' for p in points)
    result = subprocess.run([str(probe.resolve()), str(eos), str(table)],
                            input=inputs, text=True, capture_output=True, check=True)
    rows = [json.loads(line) for line in result.stdout.splitlines()]
    if len(rows) != len(points):
        raise ValueError('runtime probe returned an incomplete response')
    for row in rows:
        if type(row.get('covered')) is not bool:
            raise ValueError('invalid runtime coverage flag')
        if row['covered']:
            if any(not math.isfinite(row[k]) or row[k] <= 0 for k in ['T', 'Pgas', 'P', 'rho']):
                raise ValueError('invalid runtime state')
            for key in ['logarithmic_thermal_derivatives', 'composition_derivatives']:
                if len(row[key]) != 4 or not all(math.isfinite(v) for v in row[key]):
                    raise ValueError('invalid runtime derivatives')
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['table', 'eos_family', 'eos_probe', 'grid_probe', 'output']:
        p.add_argument(name, type=Path)
    p.add_argument('--heldouts', type=Path, nargs='+', default=[])
    a = p.parse_args()
    manifest = a.table.with_suffix('.manifest.json')
    assembled = json.loads(manifest.read_text())
    if digest(a.table) != assembled['table_sha256']:
        raise ValueError('assembled table changed')
    dependencies = {**assembled['input_files_sha256'],
                    str(a.table.resolve()): digest(a.table),
                    str(manifest.resolve()): digest(manifest),
                    str(a.eos_probe.resolve()): digest(a.eos_probe),
                    str(a.grid_probe.resolve()): digest(a.grid_probe)}
    for path in a.eos_family.parent.glob('*.dat'):
        dependencies[str(path.resolve())] = digest(path)
    for name, expected in dependencies.items():
        if digest(name) != expected:
            raise ValueError('assembled source dependency changed')
    base_manifest = Path(assembled['base_manifest'])
    base = json.loads(base_manifest.read_text())
    identity = physical_identity(base, base['provenance'])
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        old_table = root/'base.dat'
        states = import_grid(base_manifest, old_table)
        old_points = list(states)
        old_count = len(states)
        subprocess.run([sys.executable, '-B', str(Path(__file__).with_name('assemble_nongrey_grid.py')),
                        str(base_manifest), str(root/'reassembled.dat'), '--continuations',
                        *[r['work'] for r in assembled['extension']]],
                       capture_output=True, text=True, check=True)
        if digest(root/'reassembled.dat') != digest(a.table):
            raise ValueError('source reassembly does not reproduce the candidate')
        for record in assembled['extension']:
            key, state, _, _, files = load_continuation(record['work'], identity)
            dependencies.update(files)
            states[key] = state
        axes = [base[k] for k in ['hydrogen', 'helium3', 'teff_K', 'log_g']]
        # One off-centre point in every original cell, in addition to all
        # original knots and closed edges. The latter test masked-edge logic.
        interior = []
        fraction = (3-math.sqrt(5))/2
        for indices in itertools.product(*(range(len(axis)-1) for axis in axes)):
            point = []
            for k, (axis, j) in enumerate(zip(axes, indices)):
                lo, hi = axis[j:j+2]
                point.append(math.exp((1-fraction)*math.log(lo)+fraction*math.log(hi))
                             if k == 2 else (1-fraction)*lo+fraction*hi)
            interior.append(tuple(point))
        reference_points = old_points+interior
        before = runtime(a.grid_probe, a.eos_family, old_table, reference_points)
        after = runtime(a.grid_probe, a.eos_family, a.table, reference_points)
        if any(not row['covered'] for row in before+after):
            raise ValueError('candidate lost original supported boundary states')
        retained = {k: max(abs(new[k]/old[k]-1) for old, new in zip(before, after))
                    for k in ['T', 'Pgas', 'P', 'rho']}
        if max(retained.values()) > 1e-12:
            raise ValueError('candidate changed the original boundary interpolant')
        derivative_changes = []
        for point, old, new in zip(reference_points, before, after):
            changes = {k: max(abs(a-b) for a, b in zip(new[k], old[k]))
                       for k in ['logarithmic_thermal_derivatives', 'composition_derivatives']}
            if max(changes.values()) > 1e-12:
                derivative_changes.append({'coordinates': point, 'max_abs_changes': changes})
        covered = runtime(a.grid_probe, a.eos_family, a.table, list(states))
        inputs = ''.join(' '.join(format(v, '.17g') for v in [*key[:2], state['T'], state['Pgas'], state['source_density']])+'\n'
                         for key, state in states.items())
        result = subprocess.run([str(a.eos_probe.resolve()), '--family', str(a.eos_family)],
                                input=inputs, capture_output=True, text=True, check=True)
        eos = []
        for (key, state), line in zip(states.items(), result.stdout.splitlines(), strict=True):
            values = list(map(float, line.split()))
            expected = [*key[:2], state['T'], state['Pgas'], state['source_density']]
            if len(values) != 7 or values[:5] != expected or values[5] <= 0 or not all(math.isfinite(v) for v in values):
                raise ValueError('invalid EOS probe response')
            eos.append({'coordinates': key, 'source_density': values[4],
                        'ember_density': values[5], 'relative_difference': values[6]})
        comparisons = []
        for work in a.heldouts:
            key, state, _, record, files = load_continuation(work, identity)
            if key in states:
                raise ValueError('heldout coordinate is already a source node')
            dependencies.update(files)
            estimate = runtime(a.grid_probe, a.eos_family, a.table, [key])[0]
            if not estimate['covered']:
                raise ValueError('heldout is outside a complete runtime stencil')
            comparisons.append({'coordinates': key, 'work': str(work.resolve()),
                                'direct': state, 'runtime': estimate,
                                'relative_difference': {k: estimate[k]/state[k]-1 for k in ['T', 'Pgas']},
                                'density_difference_including_EOS': estimate['rho']/state['source_density']-1,
                                'validation_sha256': digest(work/'final/validated.json')})
    if any(digest(path) != expected for path, expected in dependencies.items()):
        raise ValueError('source inputs changed during the audit')
    report = {'scope': __doc__, 'script_sha256': digest(__file__),
              'retained_models': old_count, 'source_models': len(states),
              'runtime_supported_source_nodes': sum(r['covered'] for r in covered),
              'complete_cells': assembled['complete_cells'],
              'retained_runtime_queries': len(reference_points),
              'retained_max_relative_difference': retained,
              'retained_knot_derivative_changes': derivative_changes,
              'density_comparisons': eos,
              'density_relative_difference_range': [min(r['relative_difference'] for r in eos), max(r['relative_difference'] for r in eos)],
              'heldout_comparisons': comparisons,
              'heldout_status': 'local comparisons recorded' if comparisons else 'pending; no extension acceptance implied',
              'input_files_sha256': dependencies}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: report[k] for k in ['source_models', 'runtime_supported_source_nodes', 'retained_max_relative_difference', 'density_relative_difference_range', 'heldout_status']}))


if __name__ == '__main__':
    main()
