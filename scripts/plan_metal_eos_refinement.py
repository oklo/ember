#!/usr/bin/env python3
"""Choose candidate EOS composition nodes from direct-source curvature.

This inexpensive pilot tests six thermal states. It does not validate a
runtime potential family or replace fresh heldout and profile audits.
"""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess

from metal_eos_composition import mixture

STATES = [(4157., 7.2e-5), (6031., 7.3e-4), (27183., .037),
          (2341000., 13.7), (6783000., 243.1), (9731000., 260.7)]
NAMES = ['P', 'E', 'cv', 'cp', 'grad_ad']
# Leave room for thermal interpolation and isotope-axis errors in the
# subsequent runtime audit, which uses the original acceptance thresholds.
LIMITS = [.0005, .001, .0015, .0015, .0015]


def observables(jet, temperature, density):
    pressure, energy, cv, pressure_rho, pressure_temperature = jet
    cp = cv + pressure_temperature**2 / (pressure_rho*density*temperature)
    grad_ad = pressure*pressure_temperature/(pressure_rho*density*temperature*cp)
    return [pressure, energy, cv, cp, grad_ad]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source_probe', type=Path)
    p.add_argument('work', type=Path)
    p.add_argument('--hydrogen', type=float, nargs='+',
                   default=[0, .001, .01, .025, .05, .075, .1])
    a = p.parse_args()
    if len(a.hydrogen) < 2 or a.hydrogen != sorted(set(a.hydrogen)):
        raise ValueError('candidate hydrogen axis must increase')
    for x in a.hydrogen:
        mixture(x, .12)
    a.work.mkdir(parents=True, exist_ok=False)
    probe_sha = hashlib.sha256(a.source_probe.read_bytes()).hexdigest()
    cache = {}

    def source(x, y):
        key = x, y
        if key not in cache:
            m = mixture(x, y)
            scale = m['source_mass_scale']
            request = ' '.join(map(str, m['eps']))+'\n3 1 -2\n'+''.join(
                f'{math.log(scale*rho):.17g} {math.log(T):.17g}\n' for T, rho in STATES)
            result = subprocess.run([str(a.source_probe.resolve())], input=request,
                                    text=True, capture_output=True, check=True, timeout=600)
            rows = [list(map(float, line.split())) for line in result.stdout.splitlines()]
            if len(rows) != len(STATES) or any(len(r) != 22 or r[0] != 0
                    or not all(math.isfinite(v) for v in r) for r in rows):
                raise ValueError('invalid direct-source response')
            for row, (T, rho) in zip(rows, STATES, strict=True):
                if abs(row[2]/(scale*rho)-1) > 1e-9 or abs(row[3]/T-1) > 1e-10:
                    raise ValueError('source did not reach requested state')
            # These five quantities are linear in the free-energy jet at
            # fixed T,rho. Derive cp and grad_ad after mixing; they are not
            # quantities the runtime independently interpolates.
            jets = [[r[4], r[5]*scale, r[12]*scale, r[4]*r[7], r[4]*r[8]] for r in rows]
            cache[key] = {'x': x, 'y': y, 'input': request, 'stdout': result.stdout,
                          'stderr': result.stderr, 'jets': jets}
        return cache[key]['jets']

    xs = a.hydrogen.copy()
    rounds = []
    for iteration in range(8):
        additions, rejected = [], []
        for lo, hi in zip(xs, xs[1:]):
            peak, worst = 0., None
            for fraction in [1/3, .5, 2/3]:
                x = round(lo+(hi-lo)*fraction, 12)
                for y in [0, .005, .06, .12]:
                    truth = source(x, y)
                    corners = [source(xx, yy) for xx, yy in
                               [(lo, 0), (hi, 0), (lo, .12), (hi, .12)]]
                    u, v = (x-lo)/(hi-lo), y/.12
                    weights = [(1-u)*(1-v), u*(1-v), (1-u)*v, u*v]
                    for i, (T, rho) in enumerate(STATES):
                        jet = [sum(w*c[i][k] for w, c in zip(weights, corners)) for k in range(5)]
                        errors = [abs(predicted/exact-1)/limit for predicted, exact, limit in
                                  zip(observables(jet, T, rho), observables(truth[i], T, rho), LIMITS)]
                        if max(errors) > peak:
                            peak, worst = max(errors), [x, y, T, rho, errors]
            if peak > 1:
                additions.append(round((lo+hi)/2, 12))
                rejected.append({'interval': [lo, hi], 'limit_fraction': peak, 'worst': worst})
        rounds.append(rejected)
        print('round', iteration, 'new nodes', additions, flush=True)
        if not additions:
            break
        xs = sorted(xs+additions)
        if len(xs) > 40:
            raise ValueError('excessive refinement; inspect thermal and isotope dependence')
    else:
        raise ValueError('refinement did not converge')
    raw = gzip.compress(json.dumps(list(cache.values()), allow_nan=False).encode(), mtime=0)
    (a.work/'source.json.gz').write_bytes(raw)
    report = {'scope': __doc__, 'hydrogen': xs,
              'new_hydrogen': [x for x in xs if x not in a.hydrogen], 'states': STATES,
              'limits': dict(zip(NAMES, LIMITS)), 'rounds': rounds,
              'source_queries': len(cache)*len(STATES), 'source_probe_sha256': probe_sha,
              'raw_source_sha256': hashlib.sha256(raw).hexdigest(),
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    if hashlib.sha256(a.source_probe.read_bytes()).hexdigest() != probe_sha:
        raise ValueError('source executable changed during planning')
    (a.work/'plan.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'rounds'}, indent=2))


if __name__ == '__main__':
    main()
