#!/usr/bin/env python3
"""Check preserved source planes, isotope mapping and opacity responses.

This checks the implementation and source support of an extended family.
It does not estimate physical opacity errors or interpolation error against
an independently computed intermediate-composition opacity table.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_planes(path):
    lines = path.read_text().splitlines()
    nx, nt, nr = map(int, lines[0].split()[:3])
    if len(lines) != 3+nx*(nt+1):
        raise ValueError('invalid table row count')
    return lines[1:3], {lines[3+i*(nt+1)]: lines[4+i*(nt+1):4+i*(nt+1)+nt] for i in range(nx)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['original', 'extended', 'probe', 'output']:
        p.add_argument(name, type=Path)
    p.add_argument('--tracks', type=Path, nargs='+', required=True)
    a = p.parse_args()
    preserved = []
    for z in [10, 20, 30]:
        for kind in ['low', 'high']:
            name = f'tops_gs98_mixture_z{z:03d}_{kind}.dat'
            old_axes, old = source_planes(a.original/name)
            new_axes, new = source_planes(a.extended/name)
            if old_axes != new_axes or any(new.get(k) != v for k, v in old.items()):
                raise ValueError('original source values or density support changed')
            preserved.append({'table': name, 'original_sha256': sha(a.original/name),
                              'extended_sha256': sha(a.extended/name), 'retained_planes': len(old)})
    points = []
    for path in a.tracks:
        track = json.loads(path.read_text())
        if not track['converged'] or not track['metal_inventory'].startswith('GS98'):
            raise ValueError('completed GS98 track required')
        points.extend([r[5], r[6], r[3], r[2]] for r in track['profile'])
    def query(directory, points):
        request = ''.join(' '.join(f'{v:.17g}' for v in r)+'\n' for r in points)
        raw = subprocess.check_output([str(a.probe.resolve()), str(directory)], input=request.encode())
        rows = [list(map(float, line.split())) for line in raw.decode().splitlines()]
        if len(rows) != len(points) or any(len(r) != 11 or r[:4] != q or
                not all(math.isfinite(v) for v in r) or r[4] <= 0 or not r[9] <= r[3] <= r[10]
                for r, q in zip(rows, points, strict=True)):
            raise ValueError('invalid or unsupported runtime state')
        return raw, rows
    before, _ = query(a.original, points)
    after, _ = query(a.extended, points)
    if before != after:
        raise ValueError('opacity extension changed archived stellar profiles')
    states = [(4157., 7.2e-5), (6031., 7.3e-4), (27183., .037),
              (2341000., 13.7), (6783000., 243.1), (9731000., 260.7)]
    trial = [[x, y, T, rho] for x in [.1, .15, .2, .25, .3]
             for y in [0., .005, .06, .12] for T, rho in states]
    _, rows = query(a.extended, trial)
    selected = [r for r in rows if 0 < r[1] < .12]
    h = 1e-6
    differences = []
    for row in selected:
        for axis in range(4):
            for sign in [1, -1]:
                q = row[:4].copy()
                q[axis] = q[axis]+sign*h if axis < 2 else q[axis]*math.exp(sign*h)
                differences.append(q)
    _, varied = query(a.extended, differences)
    error = 0.
    for i, row in enumerate(selected):
        for axis, index in enumerate([7, 8, 5, 6]):
            plus, minus = varied[8*i+2*axis:8*i+2*axis+2]
            measured = (math.log(plus[4])-math.log(minus[4]))/(2*h)
            error = max(error, abs(measured-row[index])/max(1., abs(row[index])))
    rejected = subprocess.run([str(a.probe.resolve()), str(a.extended)],
                              input='.08 0 6783000 243.1\n', text=True, capture_output=True)
    if rejected.returncode == 0 or 'outside table' not in rejected.stderr:
        raise ValueError('missing low-hydrogen source-domain guard')
    report = {'scope': __doc__, 'passed': error < 3e-6,
              'preserved_source_tables': preserved, 'old_profile_queries': len(points),
              'byte_identical_old_profile_output': before == after,
              'old_profile_output_sha256': hashlib.sha256(before).hexdigest(),
              'new_queries': len(rows), 'response_states': len(selected),
              'maximum_scaled_response_error': error, 'hydrogen_extrapolation_rejected': True,
              'input_tracks_sha256': {str(path): sha(path) for path in a.tracks},
              'probe_sha256': sha(a.probe), 'audit_script_sha256': sha(__file__), 'new_rows': rows}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ['new_rows', 'preserved_source_tables']}, indent=2))
    if not report['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
