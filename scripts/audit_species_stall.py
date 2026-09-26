#!/usr/bin/env python3
"""Independently analyze retained abundance-solver failure and trial states."""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from decimal import Decimal, localcontext

import numpy as np

from audit_full_star_diffusion import digest


def precise_solution(record, residual):
    """Ordinary block Gaussian elimination, distinct from native flux elimination.

    Decimal arithmetic retains a weak cell's diagonal when large face blocks
    are added. The recorded double inputs are converted exactly.
    """
    with localcontext() as context:
        context.prec = 70
        D = lambda x: Decimal.from_float(float(x))
        add = lambda a, b: [[a[i][j] + b[i][j] for j in range(2)] for i in range(2)]
        sub = lambda a, b: [[a[i][j] - b[i][j] for j in range(2)] for i in range(2)]
        mul = lambda a, b: [[sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
        mv = lambda a, b: [sum(a[i][j] * b[j] for j in range(2)) for i in range(2)]
        def inv(a):
            determinant = a[0][0] * a[1][1] - a[0][1] * a[1][0]
            assert determinant != 0
            return [[a[1][1] / determinant, -a[0][1] / determinant],
                    [-a[1][0] / determinant, a[0][0] / determinant]]
        E, A, B = ([[[D(x) for x in row] for row in block] for block in record[k]]
                   for k in ('E', 'A', 'B'))
        n = len(E)
        zero = [[D(0), D(0)], [D(0), D(0)]]
        diagonal = [add(add(E[i], B[i-1] if i else zero), A[i] if i+1 < n else zero) for i in range(n)]
        upper = [sub(zero, a) for a in B]
        lower = [sub(zero, a) for a in A]
        rhs = [[-D(v) for v in row] for row in residual]
        reduced = [[row[:] for row in block] for block in diagonal]
        right = [row[:] for row in rhs]
        for i in range(1, n):
            factor = mul(lower[i-1], inv(reduced[i-1]))
            reduced[i] = sub(reduced[i], mul(factor, upper[i-1]))
            correction = mv(factor, right[i-1])
            right[i] = [right[i][j] - correction[j] for j in range(2)]
        answer = [None] * n
        answer[-1] = mv(inv(reduced[-1]), right[-1])
        for i in range(n-2, -1, -1):
            product = mv(upper[i], answer[i+1])
            answer[i] = mv(inv(reduced[i]), [right[i][j] - product[j] for j in range(2)])
        error = D(0)
        for i in range(n):
            for j in range(2):
                terms = [diagonal[i][j][k] * answer[i][k] for k in range(2)]
                if i:
                    terms += [lower[i-1][j][k] * answer[i-1][k] for k in range(2)]
                if i+1 < n:
                    terms += [upper[i][j][k] * answer[i+1][k] for k in range(2)]
                denominator = sum(map(abs, terms)) + abs(rhs[i][j])
                if denominator:
                    error = max(error, abs(sum(terms) - rhs[i][j]) / denominator)
        return np.asarray([[float(x) for x in row] for row in answer]).ravel(), float(error)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    assert not args.output.exists()
    root = Path(__file__).resolve().parents[1]
    scratch = Path('/tmp/ember-species-stall-diagnostic-v1')
    records = [json.loads(line) for line in (scratch / 'stderr.jsonl').read_text().splitlines()
               if line.startswith('{')]
    assert len(records) == 8
    base = records[0]
    assert base['label'] == 'current_derivatives'
    mass = np.asarray(base['mass'])
    direction = np.asarray(base['Newton_direction']).ravel()
    independent, _ = precise_solution(base, base['residual'])
    direction_error = float(np.max(abs(independent - direction)) / np.max(abs(direction)))
    checks = [dict(check='independent Newton direction', error=direction_error, tolerance=1e-10,
                   passed=direction_error <= 1e-10)]
    rows = []
    for record in records:
        residual = np.asarray(record['residual'])
        frozen_direction, error = precise_solution(base, residual)
        full_direction = precise_solution(record, residual)[0] if record['label'].endswith('_derivatives') else None
        checks.append(dict(check='independent linear balance', label=record['label'],
                           error=error, tolerance=1e-10, passed=error <= 1e-10))
        rows.append(dict(label=record['label'], rate_residual_norm=record['norm'],
                         frozen_Newton_correction=float(np.max(abs(frozen_direction))),
                         updated_Newton_correction=None if full_direction is None else float(np.max(abs(full_direction))),
                         maximum_rate_residual_cell=list(map(int, np.unravel_index(
                             np.argmax(abs(residual) / mass[:, None]), residual.shape))),
                         integrated_species_balance=np.sum(residual, axis=0).tolist()))
    old, full = rows[0], rows[2]
    checks.extend([
        dict(check='full trial reduces abundance correction below original tolerance',
             passed=full['updated_Newton_correction'] < 1e-15 < old['updated_Newton_correction']),
        dict(check='full trial improves species conservation', passed=max(map(abs, full['integrated_species_balance']))
             < max(map(abs, old['integrated_species_balance']))),
        dict(check='unscaled rate merit rejects that full trial', passed=full['rate_residual_norm'] > old['rate_residual_norm']),
    ])
    source_manifest = Path('/tmp/ember-full-star-diffusion-source-v1/manifest.json')
    original = json.loads(source_manifest.read_text())[str(root / 'src/species_transport.cpp')]
    assert digest(original['retained']) == original['sha256']
    paths = [Path(__file__), root / 'src/species_transport.cpp', source_manifest,
             Path(original['retained']),
             Path('/tmp/ember-species-stall-diagnostic-run-v1/receipt.json'),
             root / 'docs/results/full_star_timestep_v1.json']
    paths += [p for p in scratch.iterdir() if p.is_file()]
    paths += [root / 'docs/results/species_timestep_stall_v1.json',
              Path('/tmp/ember-species-stall-analysis-source-v1.py')]
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  outcome='verified_abundance_merit_stall' if all(c['passed'] for c in checks) else 'failed_stall_analysis',
                  accepted_for_stellar_evolution=False, checks=checks, records=rows,
                  data_scope='One retained failure; eight unchanged-state and Newton-trial evaluations.',
                  independent_arithmetic='70-digit Decimal block Gaussian elimination; exact conversion of retained double inputs.',
                  supersedes='species_timestep_stall_v1.json: two ordinary sparse-double linear solves missed the independent backward-error target. Native inputs and trials are reused.',
                  original_species_source_sha256=original['sha256'],
                  input_sha256={str(p): digest(p) for p in paths},
                  interpretation='A smaller abundance error need not decrease a stiff tiny-cell rate residual. '
                                 'A fallback measured with the current Jacobian can retain this trial without relaxing final tolerances.',
                  limitations=['The diagnostic preserves the failed scientific result and does not repeat completed stellar steps.',
                               'The updated nonlinear algorithm must separately pass full stellar and analytic controls.'])
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: result[k] for k in ('outcome', 'checks', 'records')}, indent=2))
    if not all(c['passed'] for c in checks):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
