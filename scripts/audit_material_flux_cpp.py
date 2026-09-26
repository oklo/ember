#!/usr/bin/env python3
"""Check C++ material flux conventions and evaluate saved hot stellar faces.

All kinetic coefficients and EOS replies are retained inputs. The local-force
controls encode a point gradient as a two-point difference with the same
geometric-mean temperature. Actual stellar faces are reported separately.
Neither calculation evolves the structure or selects new conductivity.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from native_material_transport import RGAS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    assert args.work.is_dir() and not args.output.exists()
    root = Path(__file__).resolve().parents[1]
    inputs = {}

    def source(path):
        p = Path(path);raw = p.read_bytes()
        inputs[str(p)] = hashlib.sha256(raw).hexdigest()
        return raw

    def report(name):
        return json.loads(source(root/'docs/results'/f'{name}.json'))

    transport = report('stellar_diffusion_pair_transport_v1')
    forces = report('stellar_diffusion_flux_comparison_3890gyr_v1')
    retained = report('smooth_eos_refined_retained_sources_v1')
    assert transport['outcome'] == 'completed_conditional_pair_transport'
    assert inputs[str(root/'docs/results/stellar_diffusion_flux_comparison_3890gyr_v1.json')] == transport['input_sha256'][str(root/'docs/results/stellar_diffusion_flux_comparison_3890gyr_v1.json')]
    profile_path = root/'docs/reports/2026-09-11/evolution_latest_profile.csv'
    source(profile_path)
    profile = list(csv.DictReader(profile_path.open()))
    old_input = Path('/tmp/ember-smooth-eos-refined-physical-run-v1/retained/retained.input')
    old_output = old_input.with_suffix('.output')
    for p in (old_input, old_output):
        source(p);assert inputs[str(p)] == retained['artifacts_sha256'][str(p)]
    cache = {}
    for row, reply in zip(old_input.read_text().splitlines(), old_output.read_text().splitlines(), strict=True):
        x, y, T, rho, chemical = map(float, row.split());d = json.loads(reply)
        if chemical == 1 and d['ok']:
            cache[x, y, T, rho] = d['values'][22:24]
    probe = Path('/tmp/ember-material-flux-build-v1/tests/material_flux_probe')
    for p in (probe, Path(__file__), root/'include/ember/material_flux.hpp',
              root/'src/material_flux.cpp', root/'src/material_transport.cpp',
              root/'scripts/material_flux_probe.cpp'):
        source(p)
    stderr = (args.work/'probe.stderr').open('x')
    raw = (args.work/'queries.jsonl').open('x')
    proc = subprocess.Popen([str(probe)], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            stderr=stderr, text=True, bufsize=1)
    e0, s0 = RGAS*1e7, RGAS
    checks, local_records, face_records = [], [], []
    maximum_input_asymmetry = 0.
    start = time.process_time()

    def check(name, error, tolerance, **meta):
        checks.append(dict(name=name, error=float(error), tolerance=tolerance,
                          passed=bool(np.isfinite(error) and abs(error) <= tolerance), **meta))

    def query(reduced, h, el, Ta, Tb, rho, area, dm, lo, hi):
        nonlocal maximum_input_asymmetry
        reduced = np.asarray(reduced)
        diagonal = np.sqrt(np.diag(reduced))
        asymmetry = float(np.max(abs(reduced-reduced.T)/diagonal[:, None]/diagonal[None, :]))
        maximum_input_asymmetry = max(maximum_input_asymmetry, asymmetry)
        if not np.isfinite(asymmetry) or asymmetry > 1e-12:
            raise ValueError('retained response is not reciprocal within roundoff')
        reduced = (reduced+reduced.T)/2
        arrays = (reduced, h, el, e0, s0, Ta, Tb, rho, area, dm, lo, hi)
        row = ' '.join(f'{float(v):.17g}' for a in arrays for v in np.asarray(a).ravel())
        proc.stdin.write(row+'\n');proc.stdin.flush()
        reply = proc.stdout.readline()
        if not reply:
            raise RuntimeError('material flux probe stopped before replying')
        value = json.loads(reply)
        raw.write(json.dumps(dict(input=row, output=value))+'\n');raw.flush()
        if not value['ok']:
            raise RuntimeError(value['error'])
        return value

    def full_matrix(record, case):
        g = np.eye(3);g[2, :2] = np.asarray(record['exchange_enthalpy_erg_g'])/e0
        g[2, 2] = case['energy_scale_erg_g']/e0
        full = s0*g@np.asarray(case['mobility_scaled_heat'])@g.T
        return (full+full.T)/2, g

    try:
        for record, driving in zip(transport['records'], forces['records'], strict=True):
            zone = record['zone'];assert driving['zone'] == zone
            T = record['temperature_K'];d = driving['gradients']['1']
            for case in record['cases']:
                full, g = full_matrix(record, case)
                reduced_force = np.r_[-np.asarray(d['chemical_acceleration_cm_s2'])/T,
                                      -case['energy_scale_erg_g']*d['dlnT_dr']/T]
                gradient = -np.linalg.solve(g.T, reduced_force)/s0
                # Same geometric-mean T and exact inverse-T difference. This
                # represents one retained local gradient, not an actual cell.
                if gradient[2] == 0:
                    Ta = Tb = T;length = 1.
                else:
                    contrast = np.copysign(.001, gradient[2])
                    Ta, Tb = T*np.exp(-contrast/2), T*np.exp(contrast/2)
                    length = e0/s0*2*np.sinh(contrast/2)/T/gradient[2]
                answer = query(case['mobility_scaled_heat'], record['exchange_enthalpy_erg_g'],
                               case['energy_scale_erg_g'], Ta, Tb, 1., 1., length,
                               [0., 0.], s0*length*gradient[:2])
                expected_h = np.asarray(record['exchange_enthalpy_erg_g'])+case['transport_enthalpy_erg_g']
                reference = case['solutions']['1:saved_density']
                mass = np.asarray(reference['independent_mass_flux_g_cm2_s'])
                energy = reference['full_material_energy_flux_erg_cm2_s']
                meta = dict(zone=zone, screening=case['screening'])
                check('transported EOS plus kinetic enthalpy',
                      np.max(abs(np.asarray(answer['enthalpy'])-expected_h))/max(np.max(abs(expected_h)), e0),
                      2e-10, **meta)
                check('zero-species-flux conductivity', answer['conductivity']/case['conductivity_zero_element_flux_cgs']-1,
                      2e-10, **meta)
                check('retained local species flux',
                      np.max(abs(np.asarray(answer['species_rate'])-mass))/max(np.max(abs(mass)), 1e-300),
                      2e-9, **meta)
                check('retained local full energy flux', (answer['material_luminosity']-energy)/max(abs(energy), 1.),
                      2e-9, **meta)
                check('single energy-flux decomposition',
                      (answer['material_luminosity']-answer['carried_luminosity']-answer['conductive_luminosity'])/
                      max(abs(answer['material_luminosity']), abs(answer['carried_luminosity']), abs(answer['conductive_luminosity']), 1.),
                      2e-10, **meta)
                local_records.append(dict(**meta, result=answer))
        by_zone = {r['zone']: r for r in transport['records']}
        for i in range(len(transport['records'])-1):
            a, b = by_zone[i], by_zone[i+1]
            pa, pb = profile[i], profile[i+1]
            Ta, Tb = float(pa['temperature_K']), float(pb['temperature_K'])
            rho = .5*(float(pa['density_g_cm3'])+float(pb['density_g_cm3']))
            radius = .5*(float(pa['radius_cm'])+float(pb['radius_cm']))
            area = 4*np.pi*radius**2
            dm = float(pb['mass_g'])-float(pa['mass_g'])
            phi = [cache[tuple(float(p[k]) for k in ['X','Y3','temperature_K','density_g_cm3'])] for p in [pa, pb]]
            stellar_luminosity = .5*(float(pa['luminosity_erg_s'])+float(pb['luminosity_erg_s']))
            for ca in a['cases']:
                cb = next(c for c in b['cases'] if c['screening'] == ca['screening'])
                full = .5*(full_matrix(a, ca)[0]+full_matrix(b, cb)[0])
                # Identity heat transformation sends the already full matrix.
                answer = query(full/s0, [0., 0.], e0, Ta, Tb, rho, area, dm, *phi)
                meta = dict(face=i, screening=ca['screening'])
                check('physical face energy decomposition',
                      (answer['material_luminosity']-answer['carried_luminosity']-answer['conductive_luminosity'])/
                      max(abs(answer['material_luminosity']), abs(answer['carried_luminosity']), abs(answer['conductive_luminosity']), 1.),
                      2e-10, **meta)
                check('nonnegative face entropy production', min(answer['entropy_production'], 0), 0., **meta)
                face_records.append(dict(**meta, both_nodes_nonconvective=not(a['convective'] or b['convective']),
                                         stellar_luminosity_erg_s=stellar_luminosity, result=answer))
    finally:
        proc.stdin.close();code = proc.wait(timeout=30);raw.close();stderr.close()
        if code:
            raise RuntimeError(f'material flux probe exited {code}')
    summaries = []
    for screening in {r['screening'] for r in face_records}:
        selected = [r for r in face_records if r['screening'] == screening and r['both_nodes_nonconvective']]
        carried = np.array([r['result']['carried_luminosity']/r['stellar_luminosity_erg_s'] for r in selected])
        thermal = np.array([r['result']['conductive_luminosity']/r['stellar_luminosity_erg_s'] for r in selected])
        summaries.append(dict(screening=screening, nonconvective_faces=len(selected),
                              carried_over_saved_luminosity_range=[float(min(carried)), float(max(carried))],
                              conductive_over_saved_luminosity_range=[float(min(thermal)), float(max(thermal))]))
    passed = all(c['passed'] for c in checks)
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(), outcome='passed' if passed else 'failed',
                  accepted_for_stellar_evolution=False, input_sha256=inputs, checks=checks,
                  local_coefficient_sets=len(local_records), physical_face_screening_sets=len(face_records),
                  maximum_normalized_import_symmetry_difference=maximum_input_asymmetry,
                  local_records=local_records, face_records=face_records, summaries=summaries,
                  new_EOS_queries=0, new_collision_integrals=0, python_cpu_seconds=time.process_time()-start,
                  raw_queries=dict(path=str(args.work/'queries.jsonl'),sha256=hashlib.sha256((args.work/'queries.jsonl').read_bytes()).hexdigest()),
                  limitations=['Saved local kinetic coefficients with fixed ionization and stationary metals.',
                               'Arithmetic mobility average across faces; no coefficient evolution or step integration.',
                               'Reported luminosity ratios compare this conditional transport with the saved structure, not a relaxed new model.'])
    with args.output.open('x') as out:
        json.dump(result, out, indent=2, allow_nan=False);out.write('\n')
    print(json.dumps(dict(outcome=result['outcome'], checks=len(checks), summaries=summaries,
                         failed=[c for c in checks if not c['passed']][:20]), indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
