#!/usr/bin/env python3
"""Compare a pair-grid interpolator with retained independent stellar integrals."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from electron_pair_table import ElectronPairTable
from electron_mixture_collision import positive_solve, symmetric


def response(ei, ee, ratio):
    joint = symmetric(ei + ratio * ee)
    return symmetric(positive_solve(joint, np.eye(len(joint))[:, :2])[:2])


def relative_operator_error(trial, reference):
    """Largest relative quadratic-form error, including mixed currents."""
    chol = np.linalg.cholesky(reference)
    left = np.linalg.solve(chol, trial-reference)
    scaled = np.linalg.solve(chol, left.T).T
    return float(np.max(abs(np.linalg.eigvalsh(symmetric(scaled)))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    start = time.process_time()
    root = Path(__file__).resolve().parents[1]
    inputs = {}
    def read(name):
        path = root/'docs/results'/name
        inputs[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        return json.loads(path.read_text())
    direct = read('electron_electron_energy_refinement_v1.json')
    mixture = read('electron_mixture_collision_v1.json')
    heat = read('stellar_electron_heat_3890gyr_v1.json')
    read('electron_pair_table_sources_v1.json')
    for name in ['electron_pair_table.py', 'electron_mixture_collision.py',
                 'audit_electron_pair_interpolation.py']:
        path = root/'scripts'/name
        inputs[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    results = []
    for degree in (1, 3):
        table = ElectronPairTable(root/'docs/results/electron_pair_table_sources_v1.json',
                                  interpolation_degree=degree)
        grid_error = max(relative_operator_error(
            table.matrices(s['eta'], s['b_thermal'])[1:, 1:],
            np.asarray(s['collision_matrix'])[1:, 1:]) for s in table.sources)
        invalid_rejections = 0
        for eta, b in [(-2.001, 10), (7.001, 10), (0, 2.999), (0, 150.001),
                       (float('nan'), 10), (0, float('inf'))]:
            try:
                table.matrices(eta, b)
            except ValueError:
                invalid_rejections += 1
        controls = []
        for mr in mixture['records']:
            s = direct['calculations'][f"zone_{mr['zone']}"]
            ee = np.asarray(s['collision_matrix'])
            ei = np.asarray(s['electron_ion_matrix'])
            trial = table.matrices(s['eta'], s['b_thermal'])
            for ratio in [0.25, 0.5, 1., mr['electron_electron_prefactor_ratio']]:
                ref = response(ei, ee, ratio)
                found = response(ei, trial, ratio)
                heat_ref = ref[1, 1] - ref[0, 1]**2/ref[0, 0]
                heat_found = found[1, 1] - found[0, 1]**2/found[0, 0]
                controls.append(dict(zone=mr['zone'], eta=s['eta'], b_thermal=s['b_thermal'],
                    prefactor_ratio=ratio, response_relative_operator_error=relative_operator_error(found, ref),
                    heat_response_relative_difference=float(heat_found/heat_ref-1),
                    pair_matrix_relative_operator_error=relative_operator_error(trial[1:, 1:], ee[1:, 1:])))
        etas = [r['eta_nonrelativistic'] for r in heat['records'] for c in r['cases']]
        bs = [c['b_thermal'] for r in heat['records'] for c in r['cases']]
        tick = time.perf_counter()
        matrices = table.matrices(etas, bs)
        elapsed = time.perf_counter()-tick
        np.linalg.cholesky(matrices[:, 1:, 1:])
        assert np.all(matrices[:, 0] == 0) and np.all(matrices[:, :, 0] == 0)
        worst = max(c['response_relative_operator_error'] for c in controls)
        results.append(dict(degree=degree, source_reproduction_error=grid_error,
            outside_domain_rejections=invalid_rejections, controls=controls,
            maximum_independent_response_error=worst, response_tolerance=0.005,
            passed=bool(worst <= 0.005 and grid_error < 2e-10 and invalid_rejections == 6),
            stellar_queries=len(etas), stellar_query_elapsed_seconds=elapsed,
            every_stellar_heat_block_positive=True, exact_momentum_null_mode=True))
    selected = next((r['degree'] for r in results if r['degree'] == 3 and r['passed']), None)
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(),
        outcome='passed_retained_controls' if selected else 'needs_refinement',
        accepted_for_stellar_evolution=False, candidate_degree=selected, results=results,
        input_sha256=inputs, cpu_seconds=time.process_time()-start, new_pair_integrals=0,
        limitations=['Three independent source states test interpolation, not full-domain physical accuracy.',
                     'New off-grid controls are required before selecting an interpolated transport model.',
                     'The source collision model omits correlations and relativistic corrections.',
                     'No conservative stellar abundance, energy or structure evolution is performed.'])
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({k: report[k] for k in ['outcome', 'candidate_degree', 'cpu_seconds']}))
    for r in results:
        print(json.dumps({k: r[k] for k in ['degree', 'maximum_independent_response_error',
            'source_reproduction_error', 'stellar_query_elapsed_seconds', 'passed']}))


if __name__ == '__main__':
    main()
