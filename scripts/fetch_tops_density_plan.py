#!/usr/bin/env python3
"""Run sequential, resumable TOPS density requests from an explicit plan.

Each completed request verifies composition, grid, and unchanged overlap.
This controller records raw-source coverage only; it never installs opacity.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from prepare_nongrey_sources import digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('plan', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    plan = json.loads(a.plan.read_text())
    requests = plan['requests']
    coordinates = [(r['X'], r['Z']) for r in requests]
    directories = [str(Path(r['work']).resolve()) for r in requests]
    if (len(set(coordinates)) != len(coordinates) or len(set(directories)) != len(directories)
            or not requests):
        raise ValueError('empty or duplicate source requests')
    inputs = {str(a.plan.resolve()): digest(a.plan), str(Path(__file__).resolve()): digest(__file__)}
    rows = []
    for request in requests:
        command = [sys.executable, '-B', str(Path(__file__).with_name('fetch_tops_density_extension.py')),
                   request['baseline_manifest'], request['work'], '--hydrogen', str(request['X']),
                   '--metallicity', str(request['Z']), '--maximum-density', str(plan['maximum_density_g_cm3']),
                   '--density-points', str(plan['density_points'])]
        subprocess.run(command, check=True)
        report = Path(request['work'])/'coverage.json'
        rows.append({'X': request['X'], 'Z': request['Z'], 'work': str(Path(request['work']).resolve()),
                     'coverage_sha256': digest(report)})
        if any(digest(f) != checksum for f, checksum in inputs.items()):
            raise ValueError('source plan changed during calculation')
        result = {'scope': __doc__, 'complete': len(rows) == len(requests),
                  'accepted_for_stellar_opacity': False, 'input_sha256': inputs, 'requests': rows}
        a.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = a.output.with_suffix('.pending')
        temporary.write_text(json.dumps(result, indent=2)+'\n')
        temporary.replace(a.output)
        print(f'completed {len(rows)} of {len(requests)} source mixtures', flush=True)


if __name__ == '__main__':
    main()
