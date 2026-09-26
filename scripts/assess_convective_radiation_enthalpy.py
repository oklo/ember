#!/usr/bin/env python3
"""Check an isobaric parcel identity and evaluate its radiation coefficient.

Uses retained EOS replies at the selected stellar nodes. It does not supply
convective parcel correlations, a heat-flux bound, or new stellar evolution.
"""
import argparse
import csv
import gzip
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from audit_full_star_diffusion import digest


def main():
    start = time.monotonic()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    assert not args.output.exists()
    root = Path(__file__).resolve().parents[1]
    profile = root / 'docs/reports/2026-09-11/evolution_latest_profile.csv'
    source_report = root / 'docs/results/native_composition_heat_v2.json'
    source = json.loads(source_report.read_text())
    for p, expected in source['artifacts_sha256'].items():
        assert digest(p) == expected, p
    scratch = Path('/tmp/ember-native-heat-audit-v2')
    queries = gzip.decompress((scratch / 'queries.txt.gz').read_bytes()).decode().splitlines()
    replies = [json.loads(s) for s in gzip.decompress((scratch / 'responses.jsonl.gz').read_bytes()).decode().splitlines()]
    with profile.open() as stream: rows = list(csv.DictReader(stream))
    assert len(rows) == 512 and len(queries) == len(replies)
    # Match Ember's adopted constant, not an independently rounded substitute.
    constants = (root / 'include/ember/constants.hpp').read_text()
    import re
    def literal(name):
        match = re.search(r'\b' + name + r'\s*=\s*([\d.eE+-]+)\s*;', constants)
        assert match, name
        return float(match.group(1))
    assert re.search(r'a_rad\s*=\s*4\.0\s*\*\s*sigma_SB\s*/\s*c\s*;', constants)
    arad = 4 * literal('sigma_SB') / literal('c')
    assert 7e-15 < arad < 8e-15
    Rgas = 8.31446261815324e7
    checks = []
    # An ideal fully ionized H/He mixture with fixed inert metals has C(X)
    # linear in X. Solve density algebraically at fixed TOTAL pressure, then
    # differentiate the bulk material+radiation enthalpy independently.
    for T in (1e4, 1e6, 1e7):
        for x in (.01, .2, .7):
            for beta in (.9999, .9, .5, .1):
                Prad = arad * T**4 / 3
                P = Prad / (1 - beta)
                C = lambda v: 2 * v + .75 * (.98 - v) + .01
                rho = (P - Prad) / (Rgas * T * C(x))
                Pc = rho * Rgas * T * 1.25
                hrad = 4 * Prad / rho * Pc / (P - Prad)
                prediction = 2.5 * Rgas * T * 1.25 + hrad
                def bulk(v):
                    density = (P - Prad) / (Rgas * T * C(v))
                    return 1.5 * Rgas * T * C(v) + arad * T**4 / density + P / density
                epsilon = 1e-4
                measured = (bulk(x + epsilon) - bulk(x - epsilon)) / (2 * epsilon)
                error = abs(measured - prediction) / abs(prediction)
                checks.append(dict(T=T, X=x, material_pressure_fraction=beta,
                                   relative_difference=error, passed=error <= 1e-10))
    nodes = []
    for i, (row, query, e) in enumerate(zip(rows, queries, replies)):
        x, y, T, rho, active0, active1 = map(float, query.split())
        assert [x, y, T, rho] == [float(row[k]) for k in ('X', 'Y3', 'temperature_K', 'density_g_cm3')]
        assert active0 == active1 == 1 and 'error' not in e
        P, _, _, _, _, chiT, chiRho, *_ = e['state']
        Prad = arad * T**4 / 3
        assert 0 < Prad < P
        Pc = rho * T * np.asarray(e['gradient_rho'])
        density_response = -Pc / (P * chiRho)
        hrad = -4 * Prad / rho * density_response
        hmaterial = np.asarray(e['enthalpy'])
        # Independent equivalent expression using total thermal expansion.
        htotal = T * ((chiT / chiRho) * np.asarray(e['gradient_rho']) - np.asarray(e['gradient_T']))
        error = float(np.max(abs(htotal - hmaterial - hrad) / (Rgas * T + abs(hmaterial))))
        checks.append(dict(node=i, identity_normalized_error=error, passed=error <= 1e-12))
        nodes.append(dict(node=i, mass_fraction=float(row['mass_g']) / float(rows[-1]['mass_g']),
                          T=T, rho=rho, radiation_pressure_fraction=Prad / P,
                          radiation_composition_enthalpy=hrad.tolist(),
                          material_composition_enthalpy=hmaterial.tolist(),
                          coefficient_normalized_size=(abs(hrad) / (Rgas * T + abs(hmaterial))).tolist()))
    envelope = nodes[397:]
    maximum = max(envelope, key=lambda n: max(n['coefficient_normalized_size']))
    inputs = [Path(__file__).resolve(), root / 'include/ember/constants.hpp', profile, source_report,
              scratch / 'queries.txt.gz', scratch / 'responses.jsonl.gz']
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  outcome='verified_isobaric_parcel_coefficient' if all(c['passed'] for c in checks) else 'failed_parcel_identity',
                  accepted_for_stellar_evolution=False, new_EOS_queries=0, checks=checks,
                  retained_nodes=nodes, mixed_region=[397, 512],
                  maximum_envelope_coefficient=maximum, elapsed_seconds=time.monotonic() - start,
                  input_sha256={str(p): digest(p) for p in inputs},
                  equation='h_rad,i = -(4 P_rad/rho) (d ln rho/d X_i)_(T,Ptotal) = (4 P_rad/rho) P_Xi/P_lnrho',
                  coefficient_units='erg per gram per unit independent mass fraction',
                  limitations=['This is a thermodynamic coefficient, not a heat-flux or cooling-age bound.',
                               'A parcel interpretation assumes local thermal equilibrium and pressure balance.',
                               'The coefficient multiplies macroscopic composition transport, not automatically the microscopic species rate.',
                               'Convective temperature and composition correlations and boundary movement still require a physical model.',
                               'The identity does not justify homogeneous reservoir convection after crystallization or at an optically thin boundary.'])
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(dict(outcome=result['outcome'],checks=len(checks),maximum_envelope_coefficient=maximum,elapsed_seconds=result['elapsed_seconds'])))
    if not all(c['passed'] for c in checks): raise SystemExit(1)


if __name__ == '__main__':
    main()
