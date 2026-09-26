#!/usr/bin/env python3
"""Replace only the quadrature rules in a checked collision-transport table.

Scattering splines and ion collision integrals are preserved byte for byte.
The selected 96/128 rules have local and coupled stellar comparisons in
docs/results/collision_quadrature_{comparison,stellar}_v1.json. Other orders
require their own accuracy checks before use in stellar evolution.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_jacobi, roots_legendre


def repack(source, output, jacobi_points, legendre_points):
    if output.exists() or output.with_suffix('.json').exists():
        raise FileExistsError('preserve existing table and source record')
    if any(n < 32 or n > 4096 for n in [jacobi_points, legendre_points]):
        raise ValueError('quadrature point counts must be between 32 and 4096')
    raw = source.read_bytes()
    lines = raw.decode().splitlines(keepends=True)
    if len(lines) < 7 or lines[0].strip() != 'EMBER_COLLISION_TRANSPORT_V1':
        raise ValueError('expected exported collision transport V1 table')
    old_counts = []
    for start in [-6, -3]:
        n = int(lines[start]); old_counts.append(n)
        x, w = [np.fromstring(lines[start+i], sep=' ') for i in [1, 2]]
        if (len(x) != n or len(w) != n or not np.isfinite(x).all()
                or not np.isfinite(w).all() or np.any(np.diff(x) <= 0)
                or np.any(w <= 0) or x[0] <= -1 or x[-1] >= 1):
            raise ValueError('invalid existing quadrature block')
    prefix = ''.join(lines[:-6])
    with output.open('x') as f:
        f.write(prefix)
        for n, jacobi in [(jacobi_points, True), (legendre_points, False)]:
            x, w = roots_jacobi(n, 0, 1.5) if jacobi else roots_legendre(n)
            f.write(str(n)+'\n')
            for values in [x, w]:
                f.write(' '.join(format(v, '.17g') for v in values)+'\n')
    if ''.join(output.read_text().splitlines(keepends=True)[:-6]) != prefix:
        raise AssertionError('scattering table prefix changed')
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    report = dict(outcome='repacked_quadrature_only', source=str(source),
                  source_sha256=hashlib.sha256(raw).hexdigest(), output=str(output),
                  output_sha256=sha(output), script_sha256=sha(Path(__file__)),
                  original_points=old_counts, points=[jacobi_points, legendre_points],
                  preserved_prefix_sha256=hashlib.sha256(prefix.encode()).hexdigest(),
                  new_scattering_integrals=0)
    output.with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--jacobi-points', type=int, required=True)
    parser.add_argument('--legendre-points', type=int, required=True)
    args = parser.parse_args()
    print(json.dumps(repack(args.source, args.output, args.jacobi_points, args.legendre_points)))
