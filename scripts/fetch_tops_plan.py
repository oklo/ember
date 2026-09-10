#!/usr/bin/env python3
"""Run a resumable, sequential TOPS source plan; never install runtime tables."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from fetch_tops_composition import fraction_label


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('specification', type=Path)
    p.add_argument('work', type=Path)
    a = p.parse_args()
    spec = json.loads(a.specification.read_text())
    seen = set()
    for job in spec['requests']:
        key = (job['X'], job['Z'])
        if key in seen or job['role'] not in ['candidate_node', 'heldout']:
            raise ValueError('duplicate mixture or unknown validation role')
        if not 0 < job['Z'] < .1 or not 0 <= job['X'] <= 1-job['Z']:
            raise ValueError('invalid source composition')
        seen.add(key)
    a.work.mkdir(parents=True, exist_ok=True)
    saved = a.work/'specification.json'
    if saved.exists() and json.loads(saved.read_text()) != spec:
        raise ValueError('source plan changed; use a new work directory')
    saved.write_text(json.dumps(spec, indent=2)+'\n')
    records = []
    for job in spec['requests']:
        x, z = job['X'], job['Z']
        print(f"request X={x:g} Z={z:g}, role={job['role']}", flush=True)
        subprocess.run([sys.executable, '-B', str(Path(__file__).with_name('fetch_tops_composition.py')),
                        str(a.work), str(x), '--metallicity', str(z)], check=True)
        label = f'tops_gs98_x{fraction_label(x,100)}_z{fraction_label(z,1000)}'
        record = json.loads((a.work/(label+'.receipt.json')).read_text())
        records.append({**record, 'validation_role': job['role']})
        report = {'scope': 'Completed raw source calculations only; no interpolation-accuracy or runtime acceptance',
                  'complete': len(records) == len(spec['requests']), 'planes': records}
        temporary = a.work/'manifest.pending'
        temporary.write_text(json.dumps(report, indent=2)+'\n')
        temporary.replace(a.work/'manifest.json')


if __name__ == '__main__':
    main()
