#!/usr/bin/env python3
"""Check a new EOS against an already accepted atmosphere table.

The accepted source and atmosphere convergence checks are reused by checksum.
Only the changed EOS support, density inversion and runtime responses are
re-evaluated; the unchanged source atmospheres are not regenerated.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import shlex

from audit_atmosphere_retention import read
from audit_nongrey_extension import runtime
from prepare_nongrey_sources import digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['acceptance', 'eos_family', 'grid_probe', 'output']:
        p.add_argument(name, type=Path)
    a = p.parse_args()
    accepted = json.loads(a.acceptance.read_text())
    if accepted.get('accepted_for_gas_trajectory') is not True:
        raise ValueError('atmosphere source has not been accepted')
    table = Path(accepted['table'])
    if digest(table) != accepted['table_sha256']:
        raise ValueError('accepted atmosphere table changed')
    inputs = {str(path.resolve()): digest(path) for path in [a.acceptance, table, a.eos_family,
                                                           a.grid_probe, Path(__file__)]}
    for name, expected in accepted['reports_sha256'].items():
        path = Path(name)
        if digest(path) != expected:
            raise ValueError('an accepted atmosphere check changed')
        inputs[str(path.resolve())] = expected
    for line in a.eos_family.read_text().splitlines()[3:]:
        path = a.eos_family.parent/shlex.split(line)[0]
        inputs[str(path.resolve())] = digest(path)
    _, axes, states = read(table)
    def point(key):
        return [key[0], key[1], 10**key[2], key[3]]
    keys = [key for key, value in states.items() if value[0]]
    if len(keys) != accepted['source_models']:
        raise ValueError('accepted source count differs')
    points = [point(key) for key in keys]
    for indices in itertools.product(*(range(len(axis)-1) for axis in axes)):
        corners = list(itertools.product(*[[axis[i], axis[i+1]] for axis, i in zip(axes, indices)]))
        if all(states[key][0] for key in corners):
            # Use an off-centre interior point, including composition and g.
            f = (3-math.sqrt(5))/2
            points.append(point([(1-f)*axis[i]+f*axis[i+1] for axis, i in zip(axes, indices)]))
    rows = runtime(a.grid_probe, a.eos_family, table, points, report_domain_errors=True)
    unsupported = [{'coordinates': q, 'response': r} for q, r in zip(points, rows, strict=True) if not r['covered']]
    knot_difference = 0.
    for key, row in zip(keys, rows):
        if row['covered']:
            # Matching values and the effective-temperature axis use log10.
            expected = states[key][1:]
            knot_difference = max(knot_difference, abs(math.log10(row['T'])-expected[0]),
                                  abs(math.log10(row['Pgas'])-expected[1]))
    if any(digest(path) != expected for path, expected in inputs.items()):
        raise ValueError('an input changed during the EOS support check')
    passed = not unsupported and knot_difference < 1e-10
    report = dict(scope=__doc__, passed=passed, input_sha256=inputs,
                  source_nodes=len(keys), interior_queries=len(points)-len(keys),
                  supported_queries=sum(row['covered'] for row in rows), unsupported=unsupported,
                  maximum_matching_log_difference=knot_difference,
                  queries=[{'coordinates': q, 'response': r} for q, r in zip(points, rows, strict=True)],
                  physical_source_checks_reused=True)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ['scope', 'input_sha256', 'queries']}))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
