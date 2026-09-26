#!/usr/bin/env python3
"""Independent C++ transport controls and comparison with retained native results."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np

from audit_conservative_material_transport import AnalyticMixture
from conservative_material_transport_v2 import implicit_step
from conservative_material_transport import positive_solve
from native_material_transport import RGAS


def numbers(*arrays):
    return ' '.join(f'{float(v):.17g}' for a in arrays for v in np.asarray(a).ravel())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--native', action='store_true')
    parser.add_argument('--reuse', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    assert args.work.is_dir() and not args.output.exists()
    start = time.process_time()
    inputs, checks, records = {}, [], []

    def source(path):
        path = Path(path)
        raw = path.read_bytes()
        inputs[str(path)] = hashlib.sha256(raw).hexdigest()
        return raw

    def report(name):
        return json.loads(source(root/'docs/results'/f'{name}.json'))

    def check(name, error, tolerance, **meta):
        checks.append(dict(name=name, error=float(error), tolerance=tolerance,
                           passed=bool(np.isfinite(error) and abs(error) <= tolerance), **meta))

    probe = Path('/tmp/ember-material-transport-cpp-build-v1/tests/material_transport_probe')
    source(probe)
    for p in ('include/ember/material_transport.hpp', 'src/material_transport.cpp',
              'src/material_thermodynamics.cpp', 'tests/material_test_eos.hpp',
              'scripts/material_transport_probe.cpp', 'scripts/audit_material_transport_cpp_v2.py',
              'scripts/conservative_material_transport.py', 'scripts/conservative_material_transport_v2.py',
              'scripts/audit_conservative_material_transport.py'):
        source(root/p)
    command = [str(probe)]
    if args.native:
        family = Path('/private/tmp/ember-smooth-eos-refined-family-v1/freeeos300_gs98_z020.dat')
        retained = report('smooth_eos_refined_retained_sources_v1')
        source(family)
        assert inputs[str(family)] == retained['inputs_sha256'][str(family)]
        command.append(str(family))
    stderr = (args.work/'probe.stderr').open('x')
    raw = (args.work/'queries.jsonl').open('x')
    proc = None
    cache = {}
    reused = 0
    if args.reuse:
        for line in source(args.reuse).decode().splitlines():
            entry = json.loads(line)
            if entry['output']['ok']:
                cache[entry['input']] = entry['output']

    def query(row, expected=True):
        nonlocal proc, reused
        if row in cache:
            reused += 1
            value = cache[row]
            raw.write(json.dumps(dict(input=row, output=value, reused=True))+'\n');raw.flush()
            return value
        if proc is None:
            proc = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                    stderr=stderr, text=True, bufsize=1)
        proc.stdin.write(row+'\n');proc.stdin.flush()
        reply = proc.stdout.readline()
        if not reply:
            raise RuntimeError('C++ probe returned no reply')
        value = json.loads(reply)
        raw.write(json.dumps(dict(input=row, output=value))+'\n');raw.flush()
        if expected and not value['ok']:
            raise RuntimeError(value['error'])
        return value

    def advance(mode, old, guess, mass, rho, k, dt, e0=1., s0=1., mixture=None, changing=False):
        count, size = old.shape
        if np.shape(k) != (count-1, size, size):
            raise ValueError('audit driver requires one matrix per face')
        prefix = f'{mode} {count} {size} {dt:.17g} 2e-11 {e0:.17g} {s0:.17g} {int(changing)} '
        if mixture:
            prefix += numbers(mixture.a, mixture.z, mixture.binding, mixture.degeneracy)+' '
        return query(prefix+numbers(rho, mass, old, guess, k))

    def conservative_checks(value, **meta):
        check('integrated species and energy conservation',
              np.max(abs(np.asarray(value['integrated_conservation_error']))), 2e-9, **meta)
        check('backward Euler entropy inequality',
              min(value['entropy_change']-value['backward_euler_entropy_bound'], 0), 2e-9, **meta)
        check('converged residual', value['residual'], 2e-11, **meta)

    try:
        if not args.native:
            rng = np.random.default_rng(9451)
            for size in (1, 2, 3):
                for count in (1, 7, 31):
                    a = rng.normal(size=(count, size, size))
                    a = a@a.swapaxes(-1, -2)+np.eye(size)
                    b = rng.normal(size=(count-1, size, size))
                    b = b@b.swapaxes(-1, -2)+np.eye(size)
                    rhs = rng.normal(size=(count, size))
                    for stiffness in (0., 1e-7, 1., 1e4):
                        k = stiffness*b
                        full = np.zeros((count*size, count*size))
                        for i in range(count):
                            s = slice(i*size, (i+1)*size)
                            full[s, s] += a[i]
                            if i+1 < count:
                                t = slice((i+1)*size, (i+2)*size)
                                full[s, s] += k[i];full[t, t] += k[i]
                                full[s, t] -= k[i];full[t, s] -= k[i]
                        reference = np.linalg.solve(full, rhs.ravel()).reshape(count, size)
                        answer = query(f'chain {count} {size} '+numbers(a, k, rhs))
                        error = np.max(abs(np.asarray(answer['answer'])-reference))/np.max(abs(reference))
                        check('independently assembled dense solution', error, 2e-9,
                              components=size, cells=count, stiffness=stiffness)
            # Unequal component scales and stiff noncommuting conductances.
            count, size = 12, 3
            scales = np.array([1e-3, 1., 1e3])
            a = rng.normal(size=(count, size, size));a = a@a.swapaxes(-1, -2)+np.eye(size)
            k = rng.normal(size=(count-1, size, size));k = k@k.swapaxes(-1, -2)+np.eye(size)
            a *= scales[:, None]*scales[None, :];k *= scales[:, None]*scales[None, :]
            uniform = np.array([.3, -.8, 1.1])/scales
            rhs = np.einsum('nij,j->ni', a, uniform)
            for stiffness in (1., 1e6, 1e12):
                answer = query(f'chain {count} {size} '+numbers(a, stiffness*k, rhs))
                check('scaled stiff uniform mode', np.max(abs((np.asarray(answer['answer'])-uniform)*scales)),
                      2e-9, stiffness=stiffness)
            # A species/heat coupled nonlinear control, checked against Python.
            count = 12
            mass = np.linspace(1, 2, count);mass /= sum(mass)
            rho = np.linspace(.7, 1.3, count)
            x = np.empty((count, 3));x[:6] = [.65, 0., .35];x[6:] = [0., .15, .85]
            old = np.c_[x[:, :2], np.log(np.linspace(.8, 1.4, count))]
            k = np.broadcast_to([[1., .2, .1], [.2, .6, -.15], [.1, -.15, .9]], (11, 3, 3)).copy()
            for degeneracy in (0., 2.):
                eos = AnalyticMixture([1, 3, 4], [1, 2, 2], [1., -.4, .2], rho, degeneracy)
                for dt in (1e-4, .05, 100.):
                    for mixing in (.01, .1):
                        interior = (1-mixing)*x+mixing*(mass@x)
                        guess = np.c_[interior[:, :2], old[:, -1]]
                        reference = implicit_step(old_primitives=old, initial_guess=guess, mass=mass,
                                                  conductance=k, dt=dt, thermodynamics=eos)
                        result = advance('analytic', old, guess, mass, rho, k, dt, mixture=eos)
                        meta = dict(degeneracy=degeneracy, dt=dt, mixing=mixing)
                        check('nonlinear Python comparison', np.max(abs(np.asarray(result['primitives'])-reference['primitives'])),
                              2e-9, **meta)
                        conservative_checks(result, **meta)
                        records.append(dict(kind='analytic', **meta, result=result))
            # Independent backward-Euler Fick solution, including initially pure cells.
            for count in (2, 12, 32):
                mass = np.linspace(1, 2, count);mass /= sum(mass)
                rho = np.ones(count)
                x = np.zeros(count);x[:count//2] = 1
                old = np.c_[x, np.zeros(count)]
                guess = np.c_[.99*x+.01*(mass@x), np.zeros(count)]
                eos = AnalyticMixture([1, 1], [1, 1], [0, 0], rho)
                k = np.broadcast_to(np.eye(2), (count-1, 2, 2)).copy()
                laplacian = np.zeros((count, count))
                for i in range(count-1):
                    laplacian[i:i+2, i:i+2] += [[1, -1], [-1, 1]]
                for dt in (1e-4, .01, 1., 100.):
                    reference = np.linalg.solve(np.diag(mass)+dt*laplacian, mass*x)
                    result = advance('analytic', old, guess, mass, rho, k, dt, mixture=eos, changing=True)
                    check('independent changing-mobility Fick solution',
                          np.max(abs(np.asarray(result['primitives'])[:, 0]-reference)), 2e-9,
                          cells=count, dt=dt)
                    conservative_checks(result, cells=count, dt=dt)
                    records.append(dict(kind='Fick', cells=count, dt=dt, result=result))
        else:
            patch = report('stellar_transport_patch_v1')
            transport = report('stellar_diffusion_pair_transport_v1')
            profile_path = root/'docs/reports/2026-09-11/evolution_latest_profile.csv'
            source(profile_path)
            assert inputs[str(profile_path)] == patch['input_sha256'][str(profile_path)]
            profile = list(csv.DictReader(profile_path.open()))
            zones = patch['zones']
            radius = np.array([float(profile[i]['radius_cm']) for i in zones])
            old = np.asarray(patch['initial_primitives']);rho = patch['fixed_density_g_cm3']
            mass = patch['normalized_cell_masses'];e0 = patch['energy_scale_erg_g']
            conductances = {}
            for record in patch['records']:
                screening = record['screening']
                if screening not in conductances:
                    matrices = []
                    for i in zones:
                        row = transport['records'][i]
                        case = next(c for c in row['cases'] if c['screening'] == screening)
                        g = np.eye(3);g[2, :2] = np.array(row['exchange_enthalpy_erg_g'])/e0
                        g[2, 2] = case['energy_scale_erg_g']/e0
                        m = RGAS*g@np.array(case['mobility_scaled_heat'])@g.T
                        matrices.append((m+m.T)/2)
                    matrices = np.asarray(matrices)
                    area = 4*np.pi*((radius[:-1]+radius[1:])/2)**2
                    conductances[screening] = .5*(matrices[:-1]+matrices[1:])*(area/np.diff(radius)/patch['closed_patch_mass_g'])[:, None, None]
                result = advance('native', old, old, mass, rho, conductances[screening],
                                 record['years']*365.25*86400, e0, RGAS)
                meta = dict(screening=screening, years=record['years'])
                check('retained native stellar segment primitives',
                      np.max(abs(np.asarray(result['primitives'])-record['final_primitives'])), 2e-10, **meta)
                reference_flux = np.asarray(record['face_species_and_scaled_energy_rate'])/patch['closed_patch_mass_g']
                check('retained native stellar segment face flux',
                      np.max(abs((np.asarray(result['face_flux'])-reference_flux)/np.max(abs(reference_flux), axis=0))),
                      2e-7, **meta)
                conservative_checks(result, **meta)
                records.append(dict(kind='native_segment', **meta, result=result))
            zero = report('native_transport_zero_species_v1')
            old = np.asarray(zero['old_primitives']);mass = np.asarray(zero['normalized_masses'])
            x = np.c_[old[:, :2], .98-old[:, :2].sum(axis=1)]
            for record in zero['records']:
                mixing = record['initial_guess_mixing']
                interior = (1-mixing)*x+mixing*(mass@x)
                guess = np.c_[interior[:, :2], old[:, -1]]
                result = advance('native', old, guess, mass, np.full(4, 1000.), np.broadcast_to(zero['prescribed_conductance'], (3, 3, 3)),
                                 1., RGAS*1e7, RGAS)
                check('retained native zero-species solution',
                      np.max(abs(np.asarray(result['primitives'])-record['final_primitives'])), 2e-10, mixing=mixing)
                conservative_checks(result, mixing=mixing)
                records.append(dict(kind='native_zero_species', mixing=mixing, result=result))
            # Compare the native adapter directly with retained EOS replies.
            raw_path = Path(patch['native_query_file']['path'])
            source(raw_path);assert inputs[str(raw_path)] == patch['native_query_file']['sha256']
            for index, line in enumerate(raw_path.read_text().splitlines()[::17]):
                saved = json.loads(line)
                x, y, T, rho, chemical = map(float, saved['input'].split())
                assert saved['result']['ok'] and chemical == 1
                v = np.asarray(saved['result']['values']);e0 = RGAS*1e7
                value = query('point '+numbers(x, y, T, rho, e0, RGAS, 1))
                cv = v[3];ec = v[19:21]
                h = np.empty((3, 3));h[:2, :2] = v[24:28].reshape(2, 2)+np.outer(ec, ec)/(cv*T*T)
                h[:2, 2] = h[2, :2] = -e0*ec/(cv*T*T);h[2, 2] = e0*e0/(cv*T*T)
                expected = positive_solve(h/RGAS, np.eye(3))
                check('native entropy capacity from retained source response',
                      np.max(abs(np.asarray(value['capacity'])-expected))/np.max(abs(expected)), 2e-10, point=index)
                check('native energy from retained source response', value['conserved'][2]/(v[1]/e0)-1, 2e-12, point=index)
    finally:
        code = 0
        if proc is not None:
            proc.stdin.close()
            code = proc.wait(timeout=30)
        raw.close();stderr.close()
        if code:
            raise RuntimeError(f'C++ probe exit {code}')
    passed = all(c['passed'] for c in checks)
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  outcome='passed' if passed else 'failed', native=args.native,
                  accepted_for_stellar_evolution=False, reused_complete_probe_replies=reused, checks=checks, records=records,
                  input_sha256=inputs, python_cpu_seconds=time.process_time()-start,
                  raw_queries=dict(path=str(args.work/'queries.jsonl'),
                                   sha256=hashlib.sha256((args.work/'queries.jsonl').read_bytes()).hexdigest()),
                  new_FreeEOS_source_calculations=0, new_collision_integrals=0,
                  limitations=['Closed fixed-density segments only; no track age advance.',
                               'Native kinetic coefficients remain frozen for the retained segment comparison.',
                               'Coupling burning, pressure work, radiation and convection remains required.'])
    with args.output.open('x') as out:
        json.dump(result, out, indent=2, allow_nan=False);out.write('\n')
    print(json.dumps(dict(outcome=result['outcome'], checks=len(checks),
                         failed=[c for c in checks if not c['passed']]), indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
