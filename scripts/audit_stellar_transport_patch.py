#!/usr/bin/env python3
"""Closed-boundary transport checks on a segment of the saved stellar structure.

Density and geometry stay fixed, with no burning, radiation flux or convection.
This is a numerical material-transport experiment, not a stellar continuation.
The initial kinetic conductances are frozen during these short comparisons.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from conservative_material_transport_v2 import implicit_step
from native_material_transport import NativeMaterialThermodynamics, RGAS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    assert args.work.is_dir() and not args.output.exists()
    root = Path(__file__).resolve().parents[1]
    inputs = {}
    def read(path):
        path = Path(path);raw = path.read_bytes()
        inputs[str(path)] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)
    for name in ('conservative_material_transport_v1', 'conservative_material_transport_v2'):
        assert read(root/'docs/results'/f'{name}.json')['outcome'] == 'passed'
    transport = read(root/'docs/results/stellar_diffusion_pair_transport_v1.json')
    assert transport['outcome'] == 'completed_conditional_pair_transport'
    profile_path = root/'docs/reports/2026-09-11/evolution_latest_profile.csv'
    inputs[str(profile_path)] = hashlib.sha256(profile_path.read_bytes()).hexdigest()
    profile = list(csv.DictReader(profile_path.open()))
    zones = list(range(292, 308))
    selected = [profile[i] for i in zones]
    radius = np.array([float(p['radius_cm']) for p in selected])
    density = np.array([float(p['density_g_cm3']) for p in selected])
    old = np.array([[float(p['X']), float(p['Y3']), np.log(float(p['temperature_K']))] for p in selected])
    star_mass = float(profile[-1]['mass_g'])
    mass = star_mass*np.array([transport['records'][i]['mass_fraction'] for i in zones])
    patch_mass = float(sum(mass))
    weights = mass/patch_mass
    energy_scale = RGAS*1e7
    thermo = NativeMaterialThermodynamics(root, args.work, density, energy_scale)
    inputs.update(thermo.input_sha256)
    for name in ('audit_stellar_transport_patch.py', 'native_material_transport.py',
                 'conservative_material_transport.py', 'conservative_material_transport_v2.py'):
        p = root/'scripts'/name;inputs[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
    checks, records = [], []
    def check(name, error, tolerance=2e-9, **meta):
        checks.append(dict(name=name, error=float(error), tolerance=tolerance,
                          passed=bool(np.isfinite(error) and abs(error) <= tolerance), **meta))
    start = time.process_time()
    try:
        for screening in [c['screening'] for c in transport['records'][0]['cases']]:
            matrices = []
            for zone in zones:
                row = transport['records'][zone]
                assert not row['convective']
                case = next(c for c in row['cases'] if c['screening'] == screening)
                conversion = np.eye(3)
                conversion[2, :2] = np.asarray(row['exchange_enthalpy_erg_g'])/energy_scale
                conversion[2, 2] = case['energy_scale_erg_g']/energy_scale
                matrix = RGAS*conversion@np.asarray(case['mobility_scaled_heat'])@conversion.T
                matrices.append((matrix+matrix.T)/2)
            matrices = np.asarray(matrices)
            area = 4*np.pi*((radius[:-1]+radius[1:])/2)**2
            conductance = .5*(matrices[:-1]+matrices[1:])*(area/np.diff(radius)/patch_mass)[:, None, None]
            for years in (10., 100., 1000.):
                dt = years*365.25*86400
                result = implicit_step(old_primitives=old, initial_guess=old, mass=weights,
                    conductance=conductance, dt=dt, thermodynamics=thermo)
                check('closed-patch species and material energy', np.max(abs(result['integrated_conservation_error'])),
                      screening=screening, years=years)
                check('closed-patch entropy inequality', min(result['entropy_change']-result['backward_euler_entropy_bound'], 0.),
                      screening=screening, years=years)
                final = result['primitives']
                record = dict(screening=screening, years=years, iterations=result['iterations'],
                    residual=result['residual'], integrated_conservation_error=result['integrated_conservation_error'].tolist(),
                    entropy_change_over_Rgas=result['entropy_change'],
                    backward_euler_entropy_bound_over_Rgas=result['backward_euler_entropy_bound'],
                    maximum_absolute_abundance_change=float(np.max(abs(final[:, :2]-old[:, :2]))),
                    maximum_fractional_temperature_change=float(np.max(abs(np.expm1(final[:, 2]-old[:, 2])))),
                    final_primitives=final.tolist(), face_species_and_scaled_energy_rate=(result['face_flux']*patch_mass).tolist())
                if years == 100.:
                    first = implicit_step(old_primitives=old, initial_guess=old, mass=weights,
                        conductance=conductance, dt=dt/2, thermodynamics=thermo)
                    second = implicit_step(old_primitives=first['primitives'], initial_guess=first['primitives'],
                        mass=weights, conductance=conductance, dt=dt/2, thermodynamics=thermo)
                    record['one_vs_two_steps_maximum_abundance_difference'] = float(np.max(abs(final[:, :2]-second['primitives'][:, :2])))
                    record['one_vs_two_steps_maximum_log_temperature_difference'] = float(np.max(abs(final[:, 2]-second['primitives'][:, 2])))
                records.append(record)
                (args.work/'progress.json').write_text(json.dumps(dict(records=records, checks=checks,
                    new_native_queries=thermo.new_queries), indent=2)+'\n')
                print(json.dumps({k:record[k] for k in ['screening','years','iterations','residual',
                    'maximum_absolute_abundance_change','maximum_fractional_temperature_change']}), flush=True)
    finally:
        thermo.close()
    check('native chemical/energy Maxwell identity', thermo.maximum_maxwell_error, 1e-9)
    passed = all(c['passed'] for c in checks)
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(),
        outcome='passed_closed_patch_controls' if passed else 'failed', accepted_for_stellar_evolution=False,
        zones=zones, closed_patch_mass_g=patch_mass, energy_scale_erg_g=energy_scale,
        initial_primitives=old.tolist(), fixed_density_g_cm3=density.tolist(),
        normalized_cell_masses=weights.tolist(), records=records, checks=checks,
        new_native_table_queries=thermo.new_queries, reused_native_queries=thermo.reused_queries,
        new_FreeEOS_source_queries=0, new_collision_integrals=0,
        python_cpu_seconds=time.process_time()-start, input_sha256=inputs,
        native_query_file=dict(path=str(args.work/'native_queries.jsonl'),
            sha256=hashlib.sha256((args.work/'native_queries.jsonl').read_bytes()).hexdigest()),
        limitations=['Sealed 16-cell interior segment with fixed density and geometry; not an evolutionary checkpoint.',
                     'Initial kinetic coefficients are frozen. The adaptive coefficient iteration is checked separately in the analytic Fick limit.',
                     'No nuclear burning, radiative flux, hydrostatic work, convective exchange or atmosphere is evolved.',
                     'The stipulated collision physics and its nonideal thermodynamic matching remain conditional.'])
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False);stream.write('\n')
    print(json.dumps(dict(outcome=report['outcome'], checks=len(checks),
        failed=[c for c in checks if not c['passed']], new_native_queries=thermo.new_queries), indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
