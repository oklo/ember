#!/usr/bin/env python3
"""Check changing face coefficients against exact discrete Fick diffusion."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np

from audit_conservative_material_transport import AnalyticMixture
from conservative_material_transport_v2 import implicit_step


def fick_conductance(q):
    x = q[:, 0]
    logit = np.log(x)-np.log1p(-x)
    difference = np.diff(x)
    potential_difference = np.diff(logit)
    mobility = np.empty_like(difference)
    separated = abs(difference) > 1e-7*np.minimum(x[:-1], x[1:])
    mobility[separated] = difference[separated]/potential_difference[separated]
    middle = (x[:-1]+x[1:])/2
    mobility[~separated] = middle[~separated]*(1-middle[~separated])
    result = np.zeros((len(x)-1, 2, 2))
    result[:, 0, 0] = mobility
    result[:, 1, 1] = 1.
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    assert not args.output.exists()
    checks, records = [], []
    def check(name, error, tolerance=2e-9, **meta):
        checks.append(dict(name=name, error=float(error), tolerance=tolerance,
                          passed=bool(np.isfinite(error) and abs(error) <= tolerance), **meta))
    for n in (2, 12, 32):
        mass = np.linspace(1., 2., n);mass /= sum(mass)
        eos = AnalyticMixture([1, 1], [1, 1], [0., 0.], np.ones(n))
        old = np.zeros((n, 2));old[:n//2, 0] = 1.
        mean = float(mass@old[:, 0])
        guess = old.copy();guess[:, 0] = .9*old[:, 0]+.1*mean
        laplacian = np.zeros((n, n))
        for i in range(n-1):
            laplacian[i, i] += 1;laplacian[i+1, i+1] += 1
            laplacian[i, i+1] -= 1;laplacian[i+1, i] -= 1
        for dt in (1e-4, .01, 1., 100.):
            result = implicit_step(old_primitives=old, initial_guess=guess,
                mass=mass, conductance=fick_conductance, dt=dt, thermodynamics=eos,
                max_iterations=150)
            reference = np.linalg.solve(np.diag(mass)+dt*laplacian, mass*old[:, 0])
            meta = dict(cells=n, dt=dt)
            check('exact discrete Fick diffusion', np.max(abs(result['primitives'][:, 0]-reference)), **meta)
            check('constant temperature with identical species', np.max(abs(result['primitives'][:, 1])), **meta)
            check('species and energy conservation', np.max(abs(result['integrated_conservation_error'])), **meta)
            check('entropy inequality with changing mobility', min(result['entropy_change']-result['backward_euler_entropy_bound'], 0), **meta)
            records.append(dict(**meta, iterations=result['iterations'], residual=result['residual'],
                minimum_final_fraction=float(min(result['primitives'][:, 0].min(), 1-result['primitives'][:, 0].max())),
                maximum_absolute_difference_from_Fick=float(np.max(abs(result['primitives'][:, 0]-reference)))))
    passed = all(c['passed'] for c in checks)
    root = Path(__file__).resolve().parents[1]
    inputs = [Path(__file__), root/'scripts/conservative_material_transport.py',
              root/'scripts/conservative_material_transport_v2.py', root/'scripts/audit_conservative_material_transport.py']
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(), outcome='passed' if passed else 'failed',
        accepted_for_stellar_evolution=False, checks=checks, records=records,
        input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        limitations=['Equal-mass binary Fick limit with analytic ideal thermodynamics; not a stellar transport prescription.',
                     'State-dependent iteration convergence is tested here, not guaranteed for arbitrary coefficient functions.'])
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False);stream.write('\n')
    print(json.dumps(dict(outcome=report['outcome'], checks=len(checks),
        failed=[c for c in checks if not c['passed']], records=records), indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
