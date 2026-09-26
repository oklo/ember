#!/usr/bin/env python3
"""Check continuous hot spectral opacity with Ember's actual interpolation.

Reuse verified scalar means, compute only new source states, write separate
single-composition tables, and query the compiled runtime. Withhold temperature
and density states independently. No stellar input is replaced.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess

from audit_tops_spectral_means import source, digest
from audit_tops_electron_dispersion import electron_moments, mean, KEV, KB, HBAR, ME


def collect(plan, prior_paths):
    prior, inputs = {}, {}
    for path in prior_paths:
        report = json.loads(path.read_text())
        for name, h in report['input_sha256'].items():
            if digest(Path(name)) != h:
                raise ValueError('prior calculation dependency changed: '+name)
            inputs[name] = h
        inputs[str(path.resolve())] = digest(path)
        for row in report['records']:
            key = row['temperature_keV'], row['density_atomic_g_cm3']
            if key in prior and row['rosseland_atomic_cm2_g'] != prior[key]['rosseland_atomic_cm2_g']:
                raise ValueError('prior means disagree')
            prior[key] = row
    records, reused = {}, 0
    for directory in plan['reused_sources']+[v['work'] for v in plan['requests']]:
        root = Path(directory)
        spectra, grey = source(root)
        receipt = json.loads((root/'receipt.json').read_text())
        if receipt['X'] != 0 or receipt['Z'] != .02:
            raise ValueError('pilot expects X=0, Z=0.02')
        source_paths = [root/n for n in ('receipt.json', 'request.json', 'source.txt')]
        composition = (root/'source.txt').read_text().split(
            'No. Fraction Mass Fraction  At. No.  Chem. Sym.  Mat ID.\n')[1].split('Temperature grid')[0]
        helium = [v.split() for v in composition.splitlines() if len(v.split()) == 5 and v.split()[3] == 'He']
        if len(helium) != 1:
            raise ValueError('missing helium fraction')
        number, mass = map(float, helium[0][:2])
        atomic_mass = 4.002602*number/mass
        for key, data in spectra.items():
            if key in records:
                raise ValueError('duplicate source state')
            if key in prior:
                if any(inputs.get(str(path.resolve())) != digest(path) for path in source_paths):
                    raise ValueError('reused state has different source identity')
                row = dict(prior[key]); row.pop('role', None)
                records[key] = row; reused += 1
                continue
            t, rho = key
            ne = rho/(atomic_mass*1.66053906660e-24)*grey[key]['free_electrons_per_ion']
            u0 = HBAR*math.sqrt(4*math.pi*4.80320471257e-10**2*ne/ME)/(t*KEV)
            electron = electron_moments(t*KEV/KB, ne)
            cutoff = u0*math.sqrt(electron['plasma_frequency_squared_ratio'])
            values = [float(mean(data, t, cutoff, electron['vstar_squared'], order)) for order in (8, 16)]
            error = abs(values[0]/values[1]-1)
            if error > 1e-7:
                raise ValueError('spectral quadrature failed')
            records[key] = {'temperature_keV': t, 'density_atomic_g_cm3': rho,
                'rosseland_atomic_cm2_g': values[1],
                'reported_conditional_rosseland_atomic_cm2_g': grey[key]['rosseland'],
                'quadrature_relative_change': error, 'electron': electron,
                'source': str(root), 'free_electron_density_cm3': ne}
        inputs.update({str(path.resolve()): digest(path) for path in source_paths})
    if any(key not in records for key in prior):
        raise ValueError('plan omits a previously calculated state')
    return records, inputs, reused


def table(path, records, tt, rr):
    lines = [f'1 {len(tt)} {len(rr)} X=0 Z=0.02 continuous spectral pilot; atomic mass basis',
             ' '.join(format(math.log10(r), '.17g') for r in rr),
             ' '.join(format(math.log10(t*KEV/KB), '.17g') for t in tt), '0 0.02']
    for t in tt:
        lines.append(' '.join(format(math.log10(records[t, r]['rosseland_atomic_cm2_g']), '.17g') for r in rr))
    path.write_text('\n'.join(lines)+'\n')


def evaluate(probe, path, points, label):
    request = path.parent/(label+'.request.txt')
    request.write_text(''.join(f'0 {t*KEV/KB:.17g} {r:.17g}\n' for t, r in points))
    result = subprocess.run([str(probe), str(path), str(path)], input=request.read_text(),
                            text=True, capture_output=True)
    (path.parent/(label+'.stdout')).write_text(result.stdout)
    (path.parent/(label+'.stderr')).write_text(result.stderr)
    result.check_returncode()
    rows = [list(map(float, line.split())) for line in result.stdout.splitlines()]
    if len(rows) != len(points) or any(len(row) != 3 or not all(map(math.isfinite, row)) or row[0] <= 0 for row in rows):
        raise ValueError('invalid runtime response')
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'prior', 'probe', 'scratch', 'output'):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    if a.output.exists() or a.scratch.exists():
        raise FileExistsError('use new output and scratch paths')
    a.scratch.mkdir(parents=True)
    spec = json.loads(a.plan.read_text())
    records, inputs, reused = collect(spec, [a.prior])
    inputs.update({str(path.resolve()): digest(path) for path in [a.plan, a.probe, Path(__file__)]})
    for name in ('audit_tops_spectral_means.py', 'audit_tops_electron_dispersion.py', 'import_tops_composition.py'):
        path = Path(__file__).with_name(name)
        if str(path.resolve()) in inputs and inputs[str(path.resolve())] != digest(path):
            raise ValueError('physical integrator differs from the previous calculation')
        inputs[str(path.resolve())] = digest(path)
    for name in ('scripts/tops_table_probe.cpp', 'src/composition.cpp', 'src/opacity_table.cpp',
                 'src/interp.cpp', 'src/opacity_blend.cpp', 'include/ember/composition.hpp',
                 'include/ember/opacity_table.hpp', 'include/ember/interp.hpp', 'include/ember/opacity_blend.hpp'):
        path = Path(name); inputs[str(path.resolve())] = digest(path)
    tt, rr = spec['node_temperatures_keV'], spec['node_densities_g_cm3']
    withheld = spec['temperature_check_keV']
    full, reduced = a.scratch/'full.dat', a.scratch/'temperature-withheld.dat'
    table(full, records, tt, rr)
    table(reduced, records, [t for t in tt if t != withheld], rr)
    node_points = [(t, r) for t in tt for r in rr]
    node_values = evaluate(a.probe, full, node_points, 'nodes')
    node_error = max(abs(row[0]/records[key]['rosseland_atomic_cm2_g']-1)
                     for key, row in zip(node_points, node_values, strict=True))
    if node_error > 1e-12:
        raise ValueError('runtime does not reproduce tabulated nodes')
    density_points = [(t, r) for t in spec['density_check_temperatures_keV']
                      for r in spec['density_check_densities_g_cm3']]
    temperature_points = [(withheld, r) for r in rr+spec['density_check_densities_g_cm3']]
    comparisons = []
    for label, path, points in [('density', full, density_points), ('temperature', reduced, temperature_points)]:
        values = evaluate(a.probe, path, points, label)
        for key, row in zip(points, values, strict=True):
            category = label if label == 'density' or key[1] in rr else 'temperature_and_density'
            comparisons.append({'kind': category, 'temperature_keV': key[0], 'density_atomic_g_cm3': key[1],
                'reference_opacity': records[key]['rosseland_atomic_cm2_g'], 'runtime_opacity': row[0],
                'relative_error': row[0]/records[key]['rosseland_atomic_cm2_g']-1,
                'runtime_dlnk_dlnT': row[1], 'runtime_dlnk_dlnRho': row[2]})
    # Check derivatives at interior off-node locations, where a symmetric
    # finite difference does not cross an interpolation interval boundary.
    derivative_points = [(math.sqrt(x*y), r) for x, y in zip(tt[1:3], tt[2:4])
                         for r in spec['density_check_densities_g_cm3']]
    queries = []
    for t, r in derivative_points:
        queries.append((t, r))
        for axis in (0, 1):
            for h in (2e-5, 1e-5):
                for sign in (-1, 1):
                    q = [t, r]; q[axis] *= math.exp(sign*h); queries.append(tuple(q))
    values = evaluate(a.probe, full, queries, 'derivatives')
    derivative_checks = []
    for index, key in enumerate(derivative_points):
        rows = values[9*index:9*(index+1)]
        for axis in (0, 1):
            pairs = rows[1+4*axis:5+4*axis]
            fd = [(math.log(pairs[2*i+1][0])-math.log(pairs[2*i][0]))/(2*h)
                  for i, h in enumerate((2e-5, 1e-5))]
            estimate = (4*fd[1]-fd[0])/3
            actual = rows[0][axis+1]
            derivative_checks.append({'temperature_keV': key[0], 'density_atomic_g_cm3': key[1],
                'axis': 'temperature' if axis == 0 else 'density', 'runtime_derivative': actual,
                'finite_difference_derivative': estimate,
                'scaled_error': abs(actual-estimate)/max(1., abs(actual)),
                'finite_difference_scaled_change': abs(fd[0]-fd[1])/max(1., abs(actual))})
    maximum = max(abs(row['relative_error']) for row in comparisons)
    derivative_error = max(row['scaled_error'] for row in derivative_checks)
    numerical_ok = node_error <= 1e-12 and derivative_error <= 1e-5
    inputs.update({str(path.resolve()): digest(path) for path in (full, reduced)})
    if any(digest(Path(name)) != h for name, h in inputs.items()):
        raise ValueError('input changed during runtime checks')
    report = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
              'reused_scalar_states': reused, 'new_scalar_states': len(records)-reused,
              'records': list(records.values()), 'comparisons': comparisons,
              'maximum_relative_interpolation_error': maximum,
              'radiative_interpolation_relative_criterion': spec['radiative_interpolation_relative_criterion'],
              'radiative_interpolation_passed': maximum <= spec['radiative_interpolation_relative_criterion'],
              'runtime_node_relative_error': node_error, 'derivative_checks': derivative_checks,
              'maximum_scaled_derivative_error': derivative_error, 'runtime_numerical_checks_passed': numerical_ok,
              'assumptions': ['Single atomic-mass X=0, Z=0.02 plane; no composition interpolation tested.',
                  'The full table includes the middle temperature; a separate table excludes it for independent checks.',
                  'Physical prescription and hot-state limitations are those of the verified finite-temperature spectral calculation.'],
              'input_sha256': inputs}
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('scope', 'records', 'comparisons', 'derivative_checks', 'assumptions', 'input_sha256')}), flush=True)
    if not numerical_ok:
        raise ValueError('runtime numerical control failed')


if __name__ == '__main__':
    main()
