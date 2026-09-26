#!/usr/bin/env python3
"""Independent block algebra and fixed-structure native burning/diffusion controls.

The native experiment seals the nonconvective hot interior at its last node.
Temperature, density, geometry and initial collision coefficients stay fixed.
It is not a new stellar model or an accuracy test of the collision physics.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.linalg import solve
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('probe', type=Path)
    p.add_argument('work', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--native', action='store_true')
    args = p.parse_args()
    assert args.work.is_dir() and not args.output.exists()
    root = Path(__file__).resolve().parents[1]
    inputs = {}

    def read(path):
        path = Path(path); raw = path.read_bytes()
        inputs[str(path)] = hashlib.sha256(raw).hexdigest()
        return raw

    for path in [args.probe, Path(__file__), root/'src/species_transport.cpp',
                 root/'include/ember/species_transport.hpp', root/'scripts/species_transport_probe.cpp',
                 root/'tests/test_species_transport.cpp', root/'src/material_flux.cpp',
                 root/'src/evolution.cpp', root/'src/nuclear_pp.cpp', root/'src/eos_smooth_mixture.cpp']:
        read(path)
    family = Path('/tmp/ember-smooth-eos-refined-family-v1/freeeos300_gs98_z020.dat')
    if args.native:
        read(family)
        assert inputs[str(family)] == '7d6bfa0b26264834307dad651852bc5ff693f91f3861b4934e582dee047fd7f5'
    checks, records = [], []

    def check(name, error, tolerance, **meta):
        checks.append(dict(name=name, error=float(error), tolerance=tolerance,
                           passed=bool(np.isfinite(error) and abs(error) <= tolerance), **meta))

    with (args.work/'probe.stderr').open('x') as log:
        proc = subprocess.Popen([str(args.probe), *([str(family)] if args.native else [])],
                                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=log, text=True)

        def query(mode, *parts):
            values = [mode]
            for part in parts:
                values.extend(format(float(x), '.17g') for x in np.asarray(part).ravel())
            line = ' '.join(values)+'\n'
            proc.stdin.write(line); proc.stdin.flush()
            reply = proc.stdout.readline()
            if not reply:
                raise RuntimeError(f'probe closed unexpectedly: {proc.poll()}')
            result = json.loads(reply)
            with (args.work/'requests.jsonl').open('a') as stream:
                stream.write(json.dumps(dict(input_sha256=hashlib.sha256(line.encode()).hexdigest(), result=result))+'\n')
            return result

        try:
            if not args.native:
                rng = np.random.default_rng(582)
                for n in (1, 2, 7, 40):
                    for stiffness in (0., 1e-5, 1., 100.):
                        e = np.array([np.diag(np.exp(rng.uniform(-2, 2, 2)))+.05*rng.normal(size=(2, 2)) for _ in range(n)])
                        h = []
                        for _ in range(n):
                            x = rng.normal(size=(2, 2)); h.append(np.eye(2)+x@x.T)
                        h = np.asarray(h)
                        k = []
                        for _ in range(n-1):
                            x = rng.normal(size=(2, 2)); k.append(stiffness*(np.eye(2)+x@x.T))
                        k = np.asarray(k).reshape(n-1, 2, 2)
                        a, b = k@h[:-1], k@h[1:]
                        rhs, d = rng.normal(size=(n, 2)), rng.normal(size=(n-1, 2))
                        matrix = np.zeros((2*n, 2*n)); dense_rhs = rhs.copy()
                        for i in range(n): matrix[2*i:2*i+2, 2*i:2*i+2] = e[i]
                        for i in range(n-1):
                            lo, hi = slice(2*i, 2*i+2), slice(2*i+2, 2*i+4)
                            matrix[lo, lo] += a[i]; matrix[hi, hi] += b[i]
                            matrix[lo, hi] -= b[i]; matrix[hi, lo] -= a[i]
                            dense_rhs[i] -= d[i]; dense_rhs[i+1] += d[i]
                        expected = solve(matrix, dense_rhs.ravel()).reshape(n, 2)
                        result = query('chain', n, e, a, b, rhs, d)
                        assert result['ok'], result
                        found = np.asarray(result['answer'])
                        check('general nonsymmetric block chain against dense solve', np.max(abs(found-expected))/np.max(abs(expected)),
                              2e-10, cells=n, stiffness=stiffness)
                        check('independent full matrix residual', np.max(abs(matrix@found.ravel()-dense_rhs.ravel()))/
                              (np.linalg.norm(matrix, np.inf)*np.max(abs(found))+np.max(abs(dense_rhs))),
                              2e-12, cells=n, stiffness=stiffness)
                for stiffness in (1., 1e4, 1e8, 1e12):
                    n = 100; mass = np.geomspace(1e-3, 1, n)
                    e = mass[:, None, None]*np.array([[1., -.1], [.2, 1.3]])
                    k = stiffness*np.geomspace(.5, 2, n-1)[:, None, None]*np.array([[1., .3], [-.2, .8]])
                    expected = np.array([.15, .002])
                    rhs = e@expected
                    result = query('chain', n, e, k, k, rhs, np.zeros((n-1, 2)))
                    assert result['ok'], result
                    check('stiff uniform mode with nonuniform masses', np.max(abs(np.asarray(result['answer'])-expected)),
                          2e-12, cells=n, stiffness=stiffness)
            else:
                data = json.loads(read(root/'docs/results/stellar_diffusion_pair_transport_v1.json'))
                assert data['outcome'] == 'completed_conditional_pair_transport'
                profile_path = root/'docs/reports/2026-09-11/evolution_latest_profile.csv'; read(profile_path)
                profile = list(csv.DictReader(profile_path.open()))
                zones = [r['zone'] for r in data['records'] if not r['convective']]
                assert zones == list(range(397))
                grid = np.array([[float(profile[i][key]) for key in ['mass_g', 'density_g_cm3', 'temperature_K', 'X', 'Y3']] for i in zones])
                n = len(zones); m, rho, T = grid[:, :3].T; old = grid[:, 3:].copy()
                weights = np.zeros(n); weights[0] = m[0]; weights[:-1] += .5*np.diff(m); weights[1:] += .5*np.diff(m)
                radius = np.array([float(profile[i]['radius_cm']) for i in zones])
                geometry = (4*np.pi*(.5*(radius[:-1]+radius[1:]))**2)**2*.5*(rho[:-1]+rho[1:])/np.diff(m)
                s0 = 1.380649e-16/1.66053906660e-24; e0 = s0*1e7
                for screening in [c['screening'] for c in data['records'][0]['cases']]:
                    full = []
                    for i in zones:
                        row = data['records'][i]; c = next(c for c in row['cases'] if c['screening'] == screening)
                        reduced = np.asarray(c['mobility_scaled_heat'])
                        assert np.max(abs(reduced-reduced.T))/np.max(abs(reduced)) < 1e-12
                        g = np.eye(3); g[2, :2] = np.asarray(row['exchange_enthalpy_erg_g'])/e0
                        g[2, 2] = c['energy_scale_erg_g']/e0
                        v = s0*g@((reduced+reduced.T)/2)@g.T
                        full.append((v+v.T)/2)
                    full = np.asarray(full)
                    face = .5*(full[:-1]+full[1:])*geometry[:, None, None]
                    for years in (1000., 1e5, 1e6, 2e7):
                        dt = years*365.25*86400
                        result = query('native', [n, n, dt, 1e-13, e0, s0, 1], np.arange(1, n+1), grid, face)
                        record = dict(screening=screening, years=years, result=result)
                        records.append(record)
                        (args.work/'progress.json').write_text(json.dumps(dict(records=records, checks=checks), indent=2)+'\n')
                        print(json.dumps(dict(screening=screening, years=years, **{k:result[k] for k in
                            ('ok', 'error', 'iterations', 'residual', 'elapsed_seconds') if k in result})), flush=True)
                        if not result['ok']:
                            check('native nonlinear convergence', 1., 0., screening=screening, years=years)
                            continue
                        final = np.asarray(result['composition']); rates = np.asarray(result['nuclear_source'])
                        found_flux = np.asarray(result['boundary_flux']); phi = np.asarray(result['chemical_potential'])
                        delta = np.column_stack((np.diff(phi, axis=0)/s0, e0*np.diff(T)/(s0*T[:-1]*T[1:])))
                        expected_flux = -np.einsum('nij,nj->ni', face, delta)[:, :2]
                        check('native face flux independent matrix application', np.max(abs(found_flux-expected_flux))/np.max(abs(expected_flux)),
                              2e-10, screening=screening, years=years)
                        integrated = np.sum(weights[:, None]*(final-old-dt*rates), axis=0)/sum(weights)
                        check('integrated species change equals nuclear burning', np.max(abs(integrated)), 2e-13,
                              screening=screening, years=years)
                        residual = weights[:, None]*(final-old-dt*rates)
                        residual[:-1] += dt*found_flux; residual[1:] -= dt*found_flux
                        raw_residual = float(np.max(abs(residual/weights[:, None])))
                        hessian = np.asarray(result['chemical_hessian'])
                        reaction_jac = np.asarray(result['nuclear_jacobian'])
                        jac = lil_matrix((2*n, 2*n))
                        for i in range(n):
                            q = slice(2*i, 2*i+2)
                            jac[q, q] = np.eye(2)-dt*reaction_jac[i]
                        for i in range(n-1):
                            lo, hi = slice(2*i, 2*i+2), slice(2*i+2, 2*i+4)
                            a = dt*face[i, :2, :2]@hessian[i]/s0
                            b = dt*face[i, :2, :2]@hessian[i+1]/s0
                            jac[lo, lo] += a/weights[i]; jac[lo, hi] -= b/weights[i]
                            jac[hi, lo] -= a/weights[i+1]; jac[hi, hi] += b/weights[i+1]
                        correction = spsolve(jac.tocsc(), -(residual/weights[:, None]).ravel())
                        check('independent Newton abundance correction', np.max(abs(correction)), 1e-13,
                              screening=screening, years=years)
                        check('correction norm agrees with C++', np.max(abs(correction))-result['abundance_correction'],
                              1e-14, screening=screening, years=years)
                        record['independent_maximum_local_residual'] = raw_residual
                        record['independent_maximum_abundance_correction'] = float(np.max(abs(correction)))
                        check('nonnegative species including He4', min(np.min(final), np.min(.98-final.sum(axis=1)), 0.), 0.,
                              screening=screening, years=years)
                        record['maximum_abundance_change'] = float(np.max(abs(final-old)))
                        if years == 2e7:
                            endpoints = {1:final}
                            refinement = []
                            for steps in (2, 4):
                                trial = grid.copy(); replies = []
                                for step in range(steps):
                                    reply = query('native', [n, n, dt/steps, 1e-13, e0, s0, 1],
                                                  np.arange(1, n+1), trial, face)
                                    assert reply['ok'], reply
                                    replies.append(reply)
                                    trial[:, 3:] = reply['composition']
                                endpoints[steps] = trial[:, 3:].copy()
                                refinement.append(dict(steps=steps, maximum_change=float(np.max(abs(trial[:, 3:]-old))),
                                                       replies=replies))
                            e12 = float(np.max(abs(endpoints[1]-endpoints[2])))
                            e24 = float(np.max(abs(endpoints[2]-endpoints[4])))
                            record['time_refinement'] = refinement
                            record['one_vs_two_steps_abundance_difference'] = e12
                            record['two_vs_four_steps_abundance_difference'] = e24
                            check('time subdivision reduces endpoint difference', max(e24-e12, 0.), 0.,
                                  screening=screening, years=years)
                            print(json.dumps(dict(screening=screening, years=years,
                                                  one_vs_two=e12, two_vs_four=e24)), flush=True)
        finally:
            proc.stdin.close()
            proc.wait(timeout=20)
            assert proc.returncode == 0, proc.returncode
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  outcome='passed' if all(c['passed'] for c in checks) else 'failed',
                  scope=__doc__, mode='native' if args.native else 'analytic',
                  accepted_for_stellar_evolution=False, checks=checks, records=records,
                  input_sha256=inputs, new_FreeEOS_source_calculations=0, new_collision_integrals=0)
    if args.native:
        report.update(zones=zones, initial_states=grid.tolist(), cell_masses_g=weights.tolist(),
                      energy_scale_erg_g=e0, entropy_scale=s0,
                      limitations=['Fixed temperatures, densities, geometry and initial mobility; no structure or energy update.',
                                   'Sealed outer boundary at the last nonconvective hot node; no envelope exchange.',
                                   'Successful implicit steps alone do not establish evolutionary time accuracy.'])
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps(dict(outcome=report['outcome'], checks=len(checks), failed=[c for c in checks if not c['passed']])), flush=True)
    if report['outcome'] != 'passed': raise SystemExit(1)


if __name__ == '__main__': main()
