#!/usr/bin/env python3
"""Compare a GS98 EOS family with fresh FreeEOS states and local identities.

These are interpolation and consistency checks of the specified source, not
physical uncertainty estimates. Entropy zero points and isotope entropy are
excluded from the direct source comparison; entropy derivatives are checked.
"""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import shlex
import subprocess

from metal_eos_composition import mixture


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def query(command, request, width):
    result = subprocess.run(command, input=request, text=True, capture_output=True, check=True)
    rows = [list(map(float, line.split())) for line in result.stdout.splitlines()]
    if any(len(r) != width or not all(math.isfinite(v) for v in r) for r in rows):
        raise ValueError('invalid probe response')
    return rows, result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['family', 'ember_probe', 'source_probe', 'output']:
        p.add_argument(name, type=Path)
    p.add_argument('--hydrogen', type=float, nargs='+', required=True)
    a = p.parse_args()
    lines = a.family.read_text().splitlines()
    if lines[0] != 'EMBER_METAL_HELMHOLTZ 1':
        raise ValueError('expected metal EOS family')
    xs = list(map(float, lines[1].split()[2:]))
    ys = list(map(float, lines[2].split()[2:]))
    if ys != [0., .12] or any(not xs[0] <= x <= xs[-1] for x in a.hydrogen):
        raise ValueError('audit mixtures outside supported family')
    inputs = {str(a.family.resolve()): sha(a.family)}
    for line in lines[3:]:
        plane = a.family.parent / shlex.split(line)[0]
        inputs[str(plane.resolve())] = sha(plane)
    states = [(4157., 7.2e-5), (6031., 7.3e-4), (27183., .037),
              (2341000., 13.7), (6783000., 243.1), (9731000., 260.7)]
    records, raw_source = [], []
    for x in a.hydrogen:
        for y in [0., .005, .06, .12]:
            m = mixture(x, y)
            scale = m['source_mass_scale']
            request = ' '.join(map(str, m['eps'])) + '\n3 1 -2\n' + ''.join(
                f'{math.log(scale * rho):.17g} {math.log(T):.17g}\n' for T, rho in states)
            source, result = query([str(a.source_probe.resolve())], request, 22)
            raw_source.append({'XH': x, 'X3': y, 'input': request,
                               'stdout': result.stdout, 'stderr': result.stderr})
            if len(source) != len(states):
                raise ValueError('incomplete source response')
            for (T, rho), r in zip(states, source, strict=True):
                if r[0] != 0 or abs(r[2] / (scale * rho) - 1) > 1e-9 or abs(r[3] / T - 1) > 1e-10:
                    raise ValueError('source did not reach requested state')
                records.append({'query': [x, y, T, rho], 'source': {
                    'P': r[4], 'E': r[5] * scale, 'cv': r[12] * scale,
                    'cp': r[13] * scale, 'grad_ad': r[14]}})
    a.output.parent.mkdir(parents=True, exist_ok=True)
    raw_path = a.output.with_suffix('.source.json.gz')
    raw_path.write_bytes(gzip.compress(json.dumps(raw_source, allow_nan=False).encode(), mtime=0))
    command = [str(a.ember_probe.resolve()), str(a.family.resolve())]
    def evaluate(points):
        request = ''.join(' '.join(f'{v:.17g}' for v in row) + '\n' for row in points)
        rows, _ = query(command, request, 18)
        if len(rows) != len(points) or any(row[:4] != point for row, point in zip(rows, points, strict=True)):
            raise ValueError('Ember probe changed requested states')
        return rows
    original = evaluate([r['query'] for r in records])
    errors = dict.fromkeys(['P', 'E', 'cv', 'cp', 'grad_ad'], 0.)
    inversion = 0.
    for record, row in zip(records, original, strict=True):
        record['relative_difference'] = {}
        for name, index in [('P', 4), ('E', 5), ('cv', 10), ('cp', 9), ('grad_ad', 11)]:
            d = row[index] / record['source'][name] - 1
            record['relative_difference'][name] = d
            errors[name] = max(errors[name], abs(d))
        inversion = max(inversion, abs(row[13] / row[3] - 1))
    # Centered differences stay inside composition intervals; derivatives
    # at piecewise-linear composition knots are intentionally not tested.
    selected = [r for r in original if r[0] not in xs and 0 < r[1] < .12]
    first_law, responses = 0., 0.
    h = 1e-6
    points = []
    for row in selected:
        for axis in [0, 1, 2, 3]:
            for sign in [1, -1]:
                q = row[:4].copy()
                q[axis] = q[axis] + sign*h if axis < 2 else q[axis]*math.exp(sign*h)
                points.append(q)
    varied = evaluate(points) if points else []
    for i, s in enumerate(selected):
        x, y, T, rho, P, E, S, chi_r, chi_t, cp, cv, ad, delta = s[:13]
        for k in range(4):
            plus, minus = varied[8*i+2*k:8*i+2*k+2]
            dP = (plus[4]-minus[4])/(2*h)
            dE = (plus[5]-minus[5])/(2*h)
            if k < 2:
                responses = max(responses, abs(dP-s[14+k])/P, abs(dE-s[16+k])/(T*cv))
            elif k == 2:
                responses = max(responses, abs(dP/P-chi_t), abs(dE/(T*cv)-1))
            else:
                responses = max(responses, abs(dP/P-chi_r))
                dS = (plus[6]-minus[6])/(2*h)
                first_law = max(first_law, abs(dE-P/rho*(1-chi_t))/(T*cv),
                                abs(dS+P*chi_t/(rho*T))/cv, abs(P*delta/(rho*T*cp*ad)-1))
    limits = {'P': .001, 'E': .002, 'cv': .003, 'cp': .003, 'grad_ad': .003}
    passed = (all(errors[k] < limits[k] for k in errors) and inversion < 1e-8
              and first_law < 1e-5 and responses < 1e-5 and bool(selected))
    if any(sha(path) != digest for path, digest in inputs.items()):
        raise ValueError('EOS inputs changed during audit')
    report = {'scope': __doc__, 'passed': passed, 'source_queries': len(records),
              'identity_states': len(selected), 'maximum_relative_source_differences': errors,
              'source_comparison_limits': limits, 'maximum_density_inversion_error': inversion,
              'maximum_first_law_error': first_law, 'maximum_response_error': responses,
              'family_files_sha256': inputs, 'ember_probe_sha256': sha(a.ember_probe),
              'source_probe_sha256': sha(a.source_probe), 'raw_source_sha256': sha(raw_path),
              'audit_script_sha256': sha(__file__), 'records': records}
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ['records', 'family_files_sha256']}, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
