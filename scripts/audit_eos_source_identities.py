#!/usr/bin/env python3
"""Scan completed EOS source planes against the importer's first-law criterion.

Convergence flags, response positivity and thermodynamic consistency are
reported separately. Passing this scan does not accept an interpolated family.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import gzip
import hashlib
import json
import math
from pathlib import Path

from eos_source_coverage import absent_source_rows


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scan(job):
    source, i, expected_x, expected_y = job
    path = source/f'plane-{i:03d}'/'source.json.gz'
    checksum = sha(path); raw = json.loads(gzip.decompress(path.read_bytes()))
    if (raw['hydrogen'], raw['helium3']) != (expected_x, expected_y):
        raise ValueError('source composition differs from the specification')
    absent = absent_source_rows(raw)
    if len(raw['data']) != len(raw['logT'])*len(raw['logQ']):
        raise ValueError('incomplete source cube')
    failures = []; maximum = 0.; nonpositive = 0; nonconverged = 0
    for k, (row, missing) in enumerate(zip(raw['data'], absent, strict=True)):
        if missing:
            continue
        if len(row) != 22 or not all(map(math.isfinite, row)):
            raise ValueError('invalid source row')
        if row[0] != 0:
            nonconverged += 1
            continue
        rho, T, P = row[2:5]; chir, chit, Er, Et, Sr, St, cp, ad, delta = row[7:16]
        defects = [rho*Er/P+chit-1, T*St/Et-1, rho*T*Sr/P+chit]
        worst = max(map(abs, defects)); maximum = max(maximum, worst)
        good = P > 0 and chir > 0 and Et > 0 and cp > 0 and ad > 0 and delta > 0
        nonpositive += not good
        it, iq = divmod(k, len(raw['logQ']))
        t, q = raw['logT'][it], raw['logQ'][iq]
        coordinate_errors = [math.log(T)-t*math.log(10),
                             math.log(rho)-(q+1.5*(t-6))*math.log(10)]
        coordinate_failed = abs(coordinate_errors[0]) > 1e-10 or abs(coordinate_errors[1]) > 1e-9
        if worst > 1e-7 or coordinate_failed:
            failures.append(dict(index=k, temperature_index=it, density_index=iq,
                                 logT=raw['logT'][it], logQ=raw['logQ'][iq],
                                 defects=defects, positive_responses=good,
                                 coordinate_errors=coordinate_errors,
                                 first_law_failed=worst > 1e-7, coordinate_failed=coordinate_failed))
    if sha(path) != checksum:
        raise ValueError('source changed during scan')
    result = dict(plane=i, hydrogen=raw['hydrogen'], helium3=raw['helium3'],
                  source=str(path), source_sha256=checksum,
                  maximum_first_law_defect=maximum, nonpositive_states=nonpositive,
                  nonconverged_states=nonconverged, absent_states=sum(absent), failures=failures)
    print(i, len(failures), maximum, flush=True)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--jobs', type=int, choices=range(1, 5), default=4)
    a = p.parse_args()
    path = a.source/'specification.json'; checksum = sha(path); spec = json.loads(path.read_text())
    jobs = [(a.source, i, x, y) for i, (x, y) in enumerate(
        (x, y) for x in spec['hydrogen'] for y in spec['helium3'])]
    with ProcessPoolExecutor(a.jobs) as pool:
        records = list(pool.map(scan, jobs))
    if sha(path) != checksum:
        raise ValueError('source specification changed')
    report = dict(scope=__doc__, passed=not any(r['failures'] for r in records), criterion=1e-7,
                  source_specification_sha256=checksum, audit_script_sha256=sha(Path(__file__)),
                  coordinate_criteria=dict(log_temperature=1e-10, log_density=1e-9),
                  coordinate_failures=sum(r['coordinate_failed'] for p in records for r in p['failures']),
                  first_law_failures=sum(r['first_law_failed'] for p in records for r in p['failures']),
                  failed_states=sum(len(r['failures']) for r in records),
                  failed_isotherms=sum(len({s['temperature_index'] for s in r['failures']}) for r in records),
                  records=records)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: report[k] for k in ['passed', 'failed_states', 'failed_isotherms']}))
    if not report['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
