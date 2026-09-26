#!/usr/bin/env python3
"""Check native EOS conservative updates from exact zero H or He3 abundances.

Conductances are prescribed for this interface control; no stellar diffusion
rate or evolutionary age is inferred from the test step.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np

from conservative_material_transport_v2 import implicit_step
from native_material_transport_v2 import NativeMaterialThermodynamics, RGAS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    assert args.work.is_dir() and not args.output.exists()
    root = Path(__file__).resolve().parents[1]
    old = np.array([[0., .002, np.log(1e7)], [.001, 0., np.log(1e7)],
                    [0., .002, np.log(1e7)], [.001, 0., np.log(1e7)]])
    mass = np.array([.1, .2, .3, .4])
    matrix = .001*np.array([[1., .2, .1], [.2, .6, -.15], [.1, -.15, .9]])
    conductance = np.broadcast_to(matrix, (3, 3, 3)).copy()
    thermo = NativeMaterialThermodynamics(root, args.work, np.full(4, 1000.), RGAS*1e7)
    results = []
    try:
        for mixing in (.01, .1):
            guess = old.copy()
            guess[:, :2] = (1-mixing)*old[:, :2]+mixing*(mass@old[:, :2])
            result = implicit_step(old_primitives=old, initial_guess=guess, mass=mass,
                conductance=conductance, dt=1., thermodynamics=thermo)
            results.append(dict(initial_guess_mixing=mixing, iterations=result['iterations'],
                residual=result['residual'], final_primitives=result['primitives'].tolist(),
                final_conserved=result['conserved'].tolist(),
                integrated_conservation_error=result['integrated_conservation_error'].tolist(),
                entropy_change=result['entropy_change'], entropy_bound=result['backward_euler_entropy_bound'],
                minimum_newly_populated_fraction=float(result['primitives'][:, :2][old[:, :2] == 0].min())))
    finally:
        thermo.close()
    checks = []
    def check(name, error, tolerance=2e-9):
        checks.append(dict(name=name, error=float(error), tolerance=tolerance,
                          passed=bool(np.isfinite(error) and abs(error) <= tolerance)))
    for result in results:
        check('closed native species and internal energy', max(abs(np.asarray(result['integrated_conservation_error']))))
        check('native entropy inequality', min(result['entropy_change']-result['entropy_bound'], 0.))
        check('positive incoming isotope', min(result['minimum_newly_populated_fraction'], 0.))
        assert result['minimum_newly_populated_fraction'] > 0
    check('independence of initial Newton guess', np.max(abs(np.asarray(results[0]['final_conserved'])-
                                                          np.asarray(results[1]['final_conserved']))))
    check('one thermal old-state query per distinct state', thermo.thermal_only_queries-2, 0.)
    # States 0 and 2, and states 1 and 3 are exact repeats, independent of mass.
    passed = all(c['passed'] for c in checks)
    inputs = [Path(__file__), root/'scripts/native_material_transport.py',
              root/'scripts/native_material_transport_v2.py', root/'scripts/conservative_material_transport_v2.py']
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(), outcome='passed' if passed else 'failed',
        accepted_for_stellar_evolution=False, old_primitives=old.tolist(), normalized_masses=mass.tolist(),
        prescribed_conductance=matrix.tolist(), records=results, checks=checks,
        thermal_only_native_queries=thermo.thermal_only_queries, new_native_table_queries=thermo.new_queries,
        new_FreeEOS_source_queries=0, input_sha256={**thermo.input_sha256,
            **{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}},
        limitations=['Native thermodynamics with prescribed conductances; this control does not supply a physical stellar transport rate.',
                     'Both isotopes have nonzero total inventory; globally absent species still need removal or a nuclear source.',
                     'No abundance floor is imposed; small-tail accuracy remains set by the residual tolerance.'])
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False);stream.write('\n')
    print(json.dumps(dict(outcome=report['outcome'], checks=checks,
        thermal_queries=thermo.thermal_only_queries, total_queries=thermo.new_queries), indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
