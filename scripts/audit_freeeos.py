#!/usr/bin/env python3
"""Audit direct FreeEOS and the potential at off-grid states and a source join.

Use the retained 4096-point CMS19/grey profile, an external FreeEOS probe,
and the small C++ potential probe. See docs/FREEEOS.md for build commands.
Every source state starts cold; failed or nonfinite states are rejected.
"""
import argparse
import concurrent.futures
import hashlib
import json
import math
from pathlib import Path
import random
import subprocess


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('source_probe', type=Path)
    ap.add_argument('potential_probe', type=Path)
    ap.add_argument('profile', type=Path)
    ap.add_argument('--potential', type=Path, default=Path('data/eos/freeeos300_hhe_x070_potential.dat'))
    a = ap.parse_args()
    model = json.loads(a.profile.read_text())
    if not model['converged'] or model['points'] != 4096 or 'CMS19' not in model['physics']['eos']:
        raise ValueError('use the converged 4096-point CMS19/grey reference profile')
    cols = model['columns']
    ir, it = cols.index('density_g_cm3'), cols.index('temperature_K')
    points = [(math.log(row[ir]), math.log(row[it])) for row in model['profile'][::8]]
    rng = random.Random(12)
    for _ in range(400):
        t, q = rng.uniform(3.31, 6.89), rng.uniform(-4.8, 1.9)
        points.append(((q+1.5*(t-6))*math.log(10), t*math.log(10)))
    eps = [.7/1.00782503, .3/4.00260325]+[0.]*18
    header = ' '.join(map(str, eps))+'\n3 1 -2\n'

    def source(point):
        lr, lt = point
        p = subprocess.run([str(a.source_probe.resolve())],
                           input=header+f'{lr:.17g} {lt:.17g}\n',
                           text=True, capture_output=True, check=True, timeout=60)
        r = [float(v) for v in p.stdout.split()]
        if len(r) != 22 or r[0] != 0 or not all(math.isfinite(v) for v in r):
            raise ValueError(f'failed source evaluation at {point}: {r[:6]}')
        if abs(math.log(r[2])-lr) > 1e-9 or abs(math.log(r[3])-lt) > 1e-10:
            raise ValueError('source returned the wrong state')
        return r

    # Fundamental-theorem check crossing the H2/H2+ fit join at 9000 K.
    # Radiation is included in both F/T and its derivative and cancels.
    n = 200
    join = [((.882+1.5*(3.95+i*.0125/n-6))*math.log(10),
             (3.95+i*.0125/n)*math.log(10)) for i in range(n+1)]
    with concurrent.futures.ThreadPoolExecutor(4) as pool:
        direct = list(pool.map(source, points))
        joined = list(pool.map(source, join))
    p = subprocess.run([str(a.potential_probe.resolve()), str(a.potential.resolve())],
                       input=''.join(f'{lr:.17g} {lt:.17g}\n' for lr, lt in points),
                       text=True, capture_output=True, check=True, timeout=60)
    interpolated = [[float(v) for v in row.split()] for row in p.stdout.splitlines()]
    if len(interpolated) != len(points) or any(len(r) != 17 or r[0] != 0 or
            not all(math.isfinite(v) for v in r) for r in interpolated):
        raise ValueError('potential probe did not return every requested state')
    errors = {}
    for name, column, index in [('P',4,1),('E',5,2),('S',6,3),
                                 ('cv',10,6),('cp',13,7),('grad_ad',14,8),('delta',15,9)]:
        values = [r[column]/r[3] if name == 'cv' else r[column] for r in direct]
        diffs = [abs(r[index]/v-1) for r, v in zip(interpolated, values)]
        k = max(range(len(diffs)), key=diffs.__getitem__)
        errors[name] = dict(max_fractional_difference=diffs[k],
                            temperature_K=direct[k][3], density_g_cm3=direct[k][2])
    local = max(max(abs(r[2]*r[9]/r[4]+r[8]-1),
                    abs(r[3]*r[12]/r[10]-1),
                    abs(r[2]*r[3]*r[11]/r[4]+r[8])) for r in direct+joined)
    potential_change = (joined[-1][5]/joined[-1][3]-joined[-1][6]
                        -joined[0][5]/joined[0][3]+joined[0][6])
    slopes = [-r[5]/r[3]+1.5*r[4]/(r[2]*r[3]) for r in joined]
    integral = (.0125*math.log(10)/n)/3 * (slopes[0]+slopes[-1]
               +4*sum(slopes[1:-1:2])+2*sum(slopes[2:-1:2]))
    output = dict(source='FreeEOS 3.0.0 EOS1 (3,1,-2); H=.7 He=.3; all states queried cold',
        potential_sha256=hashlib.sha256(a.potential.read_bytes()).hexdigest(),
        sample=dict(count=len(points), profile='Every eighth point of the 4096-point CMS19/grey reference',
                    random='400 Python Random(12) states, logT uniform 3.31..6.89, logQ uniform -4.8..1.9'),
        max_direct_local_identity_residual=local, source_comparison=errors,
        molecular_join=dict(logQ=.882, logT_interval=[3.95,3.9625], simpson_intervals=n,
            potential_change_erg_g_K=potential_change, integrated_derivative_erg_g_K=integral,
            difference_erg_g_K=potential_change-integral),
        qualification='Sampled source agreement, not uniform physical error bounds; source-fit joins are regularized by the C2 potential')
    print(json.dumps(output, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
