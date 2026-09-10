#!/usr/bin/env python3
"""Separate local composition and temperature/gravity interpolation differences.

Canonical controls at the heldout Teff/gravity bracket its composition. The
runtime interpolates log T and log Pgas, so differences add in logarithms.
This is a local diagnostic of the selected gas sources, not a physical error
bound, a global grid acceptance, or an atmosphere installation.
"""
import argparse
import itertools
import json
import math
from pathlib import Path

from assemble_nongrey_grid import load_continuation, physical_identity
from audit_nongrey_extension import runtime
from prepare_nongrey_sources import digest


def composition_weights(coordinates, heldout, axes):
    """Require the four corners of the runtime cell, even at zero weight."""
    bounds = []
    for axis, value in zip(axes[:2], heldout[:2]):
        intervals = [(lo, hi) for lo, hi in zip(axis, axis[1:]) if lo < value < hi]
        if len(intervals) != 1:
            raise ValueError('heldout must lie strictly inside a composition cell')
        bounds.append(intervals[0])
    expected = set(itertools.product(*bounds))
    if len(coordinates) != 4 or set(coordinates) != expected:
        raise ValueError('controls must supply all four unique runtime composition corners')
    fractions = [(value-lo)/(hi-lo) for value, (lo, hi) in zip(heldout[:2], bounds)]
    return [math.prod(fraction if value == hi else 1-fraction
                      for value, (lo, hi), fraction in zip(point, bounds, fractions))
            for point in coordinates]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['table', 'eos_family', 'grid_probe', 'heldout', 'output']:
        p.add_argument(name, type=Path)
    p.add_argument('--controls', type=Path, nargs=4, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('use a new diagnostic report')
    manifest = a.table.with_suffix('.manifest.json')
    assembled = json.loads(manifest.read_text())
    if digest(a.table) != assembled['table_sha256']:
        raise ValueError('assembled table changed')
    inputs = dict(assembled['input_files_sha256'])
    for path in [a.table, manifest, a.grid_probe, *a.eos_family.parent.glob('*.dat')]:
        inputs[str(path.resolve())] = digest(path)
    if any(digest(path) != expected for path, expected in inputs.items()):
        raise ValueError('pinned grid input changed')
    base = json.loads(Path(assembled['base_manifest']).read_text())
    identity = physical_identity(base, base['provenance'])
    key, direct, spec, _, files = load_continuation(a.heldout, identity)
    inputs.update(files)
    controls = []
    resolution = ['depths', 'atmosphere_frequencies', 'convection_derivative_step']
    nodes = set(itertools.product(*(base[k] for k in ['hydrogen', 'helium3', 'teff_K', 'log_g'])))
    nodes.update(tuple(row['coordinates']) for row in assembled['extension'])
    for work in a.controls:
        point, state, control_spec, _, files = load_continuation(work, identity)
        if point[2:] != key[2:]:
            raise ValueError('controls must share heldout temperature and gravity')
        if any(control_spec[k] != spec[k] for k in resolution):
            raise ValueError('canonical numerical resolution differs')
        inputs.update(files)
        controls.append({'work': str(work.resolve()), 'coordinates': point, 'direct': state})
    if key in nodes or any(tuple(r['coordinates']) in nodes for r in controls):
        raise ValueError('comparison must use independent source coordinates')
    weights = composition_weights([r['coordinates'][:2] for r in controls], key, assembled['axes'])
    estimates = runtime(a.grid_probe, a.eos_family, a.table,
                        [key, *[r['coordinates'] for r in controls]])
    if any(not r['covered'] for r in estimates):
        raise ValueError('comparison lies outside complete runtime stencils')
    for row, weight, estimate in zip(controls, weights, estimates[1:]):
        row.update(weight=weight, runtime=estimate,
                   relative_difference={k: estimate[k]/row['direct'][k]-1 for k in ['T', 'Pgas']})
    quantities = {}
    for name in ['T', 'Pgas']:
        composition_log = sum(w*math.log10(r['direct'][name]) for w, r in zip(weights, controls))
        runtime_controls_log = sum(w*math.log10(r[name]) for w, r in zip(weights, estimates[1:]))
        grid_log, direct_log = math.log10(estimates[0][name]), math.log10(direct[name])
        closure = runtime_controls_log-grid_log
        if abs(closure) > 1e-12:
            raise ValueError('runtime composition interpolation does not close')
        differences = {'total': grid_log-direct_log,
                       'composition_at_fixed_temperature_gravity': composition_log-direct_log,
                       'temperature_gravity_at_composition_corners': runtime_controls_log-composition_log}
        quantities[name] = {
            'direct': direct[name], 'composition_only': 10**composition_log,
            'full_grid': estimates[0][name], 'log10_difference': differences,
            'relative_difference': {k: math.expm1(math.log(10)*v) for k, v in differences.items()},
            'log10_additivity_residual': differences['total']-sum(v for k, v in differences.items() if k != 'total'),
            'runtime_composition_closure_log10': closure}
    if any(digest(path) != expected for path, expected in inputs.items()):
        raise ValueError('source input changed during audit')
    report = {'scope': __doc__, 'coordinates': key, 'heldout_work': str(a.heldout.resolve()),
              'direct': direct, 'controls': controls, 'quantities': quantities,
              'interpretation': 'Relative differences multiply; logarithmic differences add. The temperature/gravity term is averaged over these four composition corners. Depth, source opacity interpolation and physical prescription errors are not independently bounded here.',
              'script_sha256': digest(__file__), 'input_files_sha256': inputs}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'coordinates': key, 'quantities': quantities}))


if __name__ == '__main__':
    main()
