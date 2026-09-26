#!/usr/bin/env python3
"""Record small CPU/history snapshots of explicitly named MESA processes.

Usage: monitor_mesa_comparison.py --configuration FILE --output FILE
The configuration lists native PIDs and run directories; it never discovers,
starts, stops or changes a calculation. Readouts are approximate live samples,
not exact per-model timings. Final /usr/bin/time receipts remain authoritative.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import time


def cpu_seconds(value):
    days, clock = value.split('-', 1) if '-' in value else ('0', value)
    total = 0.
    for field in clock.split(':'):
        total = 60*total + float(field)
    return 86400*int(days) + total


def history_summary(directory):
    path = directory/'LOGS/history.data'
    if not path.exists():
        return dict(stage='initial_model_construction')
    lines = path.read_text().splitlines()
    if len(lines) < 7:
        return dict(stage='history_header')
    keys = lines[5].split()
    rows = []
    for line in lines[6:]:
        values = line.split()
        if len(values) != len(keys):
            continue  # possible partially written final row
        row = dict(zip(keys, map(float, values)))
        # MESA can append a restored model number; discard its replaced tail.
        while rows and rows[-1]['model_number'] >= row['model_number']:
            rows.pop()
        rows.append(row)
    if not rows:
        return dict(stage='history_header')
    fields = ['model_number', 'star_age', 'elapsed_time', 'log_dt', 'log_Teff',
              'log_L', 'log_R', 'log_g', 'center_h1', 'surface_h1',
              'surface_he4', 'total_mass_h1', 'num_zones', 'num_iters',
              'num_retries', 'diffusion_solver_steps', 'diffusion_solver_iters']
    return dict(stage='evolving', saved_models=len(rows),
                last={k: rows[-1][k] for k in fields if k in rows[-1]})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--configuration', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--interval', type=float, default=60.)
    p.add_argument('--samples', type=int, default=1440)
    args = p.parse_args()
    cases = json.loads(args.configuration.read_text())['cases']
    pids = ','.join(str(int(c['pid'])) for c in cases)
    if args.interval < 1 or args.samples < 1:
        raise ValueError('positive sample count and interval >= 1 required')
    with args.output.open('x') as out:
        for _ in range(args.samples):
            start = datetime.now(timezone.utc)
            result = subprocess.run(['/bin/ps', '-p', pids,
                '-o', 'pid=,time=,pcpu=,rss='], capture_output=True, text=True)
            if result.returncode not in (0, 1) or result.stderr.strip():
                raise RuntimeError(result.stderr or 'ps failed')
            counters = {}
            for line in result.stdout.splitlines():
                pid, clock, utilization, rss = line.split()
                counters[int(pid)] = dict(cpu_seconds=cpu_seconds(clock),
                    cpu_percent=float(utilization), resident_bytes=1024*int(rss))
            snapshot = dict(started_utc=start.isoformat(), cases=[])
            for case in cases:
                directory = Path(case['directory'])
                snapshot['cases'].append(dict(name=case['name'], pid=case['pid'],
                    counters=counters.get(case['pid']),
                    history=history_summary(directory),
                    completion_receipt_present=(directory/'completion.json').exists()))
            snapshot['finished_utc'] = datetime.now(timezone.utc).isoformat()
            out.write(json.dumps(snapshot, allow_nan=False)+'\n'); out.flush()
            if not counters:
                break
            time.sleep(args.interval)


if __name__ == '__main__':
    main()
