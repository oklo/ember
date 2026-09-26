#!/usr/bin/env python3
"""Compare a faster native density inversion with retained native controls."""
import argparse
import gzip
import hashlib
import json
import os
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


def main():
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--scratch', type=Path, required=True)
    args = ap.parse_args()
    work = args.scratch
    work.mkdir()
    out = args.output
    reservation = root / 'docs/research/fable/coordination/reservations' / (work.name + '.json')
    assert not out.exists() and not reservation.exists()
    status = dict(owner='primary', task_ids=['E-COLLISION-INVERSION-001'], status='running',
                  controller_pid=os.getpid(), threads=1, expected_memory_bytes=1000000000,
                  data_cap_bytes=30000000, scratch_paths=[str(work)], jobs=[],
                  start_utc=datetime.now(timezone.utc).isoformat())

    def save():
        for p in (reservation, work / 'receipt.json'):
            p.write_text(json.dumps(status, indent=2) + '\n')

    save()
    identities = {}

    def read(p):
        p = Path(p)
        data = p.read_bytes()
        identities[str(p)] = hashlib.sha256(data).hexdigest()
        return gzip.decompress(data) if p.suffix == '.gz' else data

    table = '/tmp/ember-collision-transport-table-v1.dat'
    old = '/tmp/ember-collision-derivative-build-v1/tests/collision_transport_probe'
    new = '/tmp/ember-collision-fast-inversion-v1/collision_probe'
    for p in (table, old, new, __file__, '/tmp/ember-collision-fast-inversion-v1/manifest.json'):
        read(p)
    reference = Path('/tmp/ember-collision-derivative-audit-v1')
    jobs, failures, maxima, counts = [], [], {}, {}
    try:
        # The scalar set includes exact-zero species and invalid requests; the
        # derivative set also exercises composition and screening derivatives.
        for mode in ('scalar', 'derivatives'):
            q = read(reference / (mode + '.queries.txt.gz'))
            expected = [json.loads(s) for s in read(reference / (mode + '.responses.jsonl.gz')).splitlines()]
            query_values = [list(map(float, line.split())) for line in q.splitlines()]
            responses = None
            for repeat, executable in enumerate((old, new, new, old, old, new)):
                cmd = [executable, table] + (['derivatives'] if mode == 'derivatives' else [])
                start = time.monotonic()
                result = subprocess.run(cmd, input=q, capture_output=True, timeout=30, check=True)
                records = [json.loads(s) for s in result.stdout.splitlines()]
                assert len(records) == len(expected)
                jobs.append(dict(mode=mode, executable=executable, repeat=repeat,
                                 wall_seconds=time.monotonic() - start,
                                 evaluation_seconds=sum(s.get('seconds', 0) for s in records)))
                if executable == new and responses is None:
                    responses = records
                    p = work / (mode + '.responses.jsonl.gz')
                    p.write_bytes(gzip.compress(result.stdout, mtime=0)); read(p)
            counts[mode] = dict(queries=len(expected), errors=sum('error' in s for s in expected))
            for index, (a, b) in enumerate(zip(expected, responses)):
                if 'error' in a or 'error' in b:
                    if ('error' in a) != ('error' in b): failures.append([mode, index, 'domain acceptance changed'])
                    continue
                for key in ('active', 'defined'):
                    if a.get(key) != b.get(key): failures.append([mode, index, key])
                if b['backward_error'] > 1e-12: failures.append([mode, index, 'linear solve residual'])

                def compare(av, bv, prefix, coordinate=None):
                    for key in ('eta', 'ne', 'b', 'energy_scale', 'conductivity', 'enthalpy', 'mobility'):
                        x, y = np.asarray(av[key]), np.asarray(bv[key])
                        # Relative abundance increments keep a trace-species
                        # partial's cancellation error in its physical scale.
                        increment = query_values[index][coordinate] if coordinate in (2,3,4) else 1.
                        if coordinate is not None: x, y = x*increment, y*increment
                        if key == 'mobility':
                            diagonal = np.sqrt(np.diag(a[key]))
                            scale = np.outer(diagonal, diagonal)
                            scale[scale == 0] = 1.
                        else:
                            scale = np.maximum(abs(np.asarray(a[key])), 1.)
                        if coordinate is not None: scale = scale + abs(x)
                        error = float(np.max(abs(x-y) / scale))
                        name = prefix + '/' + key
                        maxima[name] = max(maxima.get(name, 0), error)
                        if error > 2e-8: failures.append([mode, index, name, error])

                compare(a, b, mode)
                if mode == 'derivatives':
                    for k in range(6):
                        if a['defined'][k]: compare(a['partials'][k], b['partials'][k], 'partial_' + str(k), k)
                    if a.get('screening_partials') != b.get('screening_partials'):
                        failures.append([mode, index, 'screening changed'])
        speedups = {}
        for mode in ('scalar', 'derivatives'):
            times = {exe: np.median([s['evaluation_seconds'] for s in jobs
                                    if s['mode'] == mode and s['executable'] == exe]) for exe in (old, new)}
            speedups[mode] = dict(old_seconds=float(times[old]), new_seconds=float(times[new]),
                                 ratio=float(times[old] / times[new]))
        report = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                      outcome='passed' if not failures else 'failed', counts=counts,
                      tolerance=2e-8, normalization='Response scale plus derivative magnitude; relative composition increments; mobility diagonal scales.', maximum_normalized_differences=maxima,
                      speedups=speedups, jobs=jobs, failures=failures, input_sha256=identities,
                      limitations=['Retained independent physical and derivative controls are reused.',
                                   'Timing excludes the stellar solve and its EOS table load.'])
        out.write_text(json.dumps(report, indent=2) + '\n')
        print(json.dumps({k: report[k] for k in ('outcome', 'counts', 'speedups')}, indent=2))
        assert not failures
        status['outcome'] = 'completed'
    finally:
        status.update(status='released', end_utc=datetime.now(timezone.utc).isoformat(), jobs=jobs)
        save()


if __name__ == '__main__':
    main()
