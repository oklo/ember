#!/usr/bin/env python3
"""Compare electron heat approximations on the saved hot stellar structure.

Both responses use prescribed fixed-ion elastic Born scattering. This is not
an accepted stellar conductivity, diffusion velocity or error bound on Ember's
selected conduction. Retain all completed zone calculations for reuse.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import math
import multiprocessing
from pathlib import Path
import re
import time

import numpy as np

from electron_ion_heat import current_basis, born_response, physical_prefactor
from electron_ion_born import KB, ME, HBAR
from ion_collision_integrals import MU

CLASSICAL_ETA = -25.
_classical_basis = None


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_new(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def zone_comparison(record, saved, metals):
    global _classical_basis
    started = time.process_time()
    if _classical_basis is None:
        _classical_basis = current_basis(CLASSICAL_ETA)
    if record['zone'] != saved['zone']:
        raise ValueError('saved profiles use different zones')
    T, rho = record['temperature_K'], record['density_baryonic_g_cm3']
    if (T != saved['temperature_K'] or rho != saved['density_g_cm3'] or
            record['mass_fraction'] != saved['mass_fraction'] or
            record['convective'] != saved['convective']):
        raise ValueError('saved thermodynamic state differs')
    fractions = np.array([record['X'], record['Y3'], record['Y4']] +
                         [(1-record['X']-record['Y3']-record['Y4'])*m[3] for m in metals])
    masses = np.array([1., 3., 4.] + [m[1] for m in metals])
    charges = np.array([1., 2., 2.] + [m[0] for m in metals])
    n = rho/MU*fractions/masses
    ne = float(n@charges)
    if np.any(n <= 0) or abs(ne/record['electron_density_cm3']-1) > 1e-12:
        raise ValueError('electron population convention differs')
    eta = saved['eta_nonrelativistic']
    basis = current_basis(eta)
    ne_current = (2*ME*KB*T)**1.5/(3*math.pi**2*HBAR**3)*basis.normalization
    density_error = abs(ne_current/ne-1)
    if density_error > 5e-10:
        raise ValueError('independent current normalization disagrees with saved density')
    charge_weight = n*charges**2
    total_weight = float(np.sum(charge_weight))
    unit_prefactor = physical_prefactor(1., 1.)
    prefactor = unit_prefactor*total_weight
    scale = basis.normalization/_classical_basis.normalization
    cases = []
    for j, old in enumerate(saved['modes'][0]['cases']):
        label, length = old['screening'], old['screening_length_cm']
        B = 8*ME*KB*T*length**2/HBAR**2
        fermi = born_response(eta, B, basis=basis)
        classical = born_response(CLASSICAL_ETA, B, basis=_classical_basis)
        I, L = np.array(fermi['collision_moments']), np.array(fermi['exact_elastic_mobility'])
        S = fermi['zero_drift_heat_response_exact']
        kappa = ne**2*KB**2*T*S/prefactor
        ratio = fermi['two_mode_to_exact_zero_drift_heat']
        classical_ratio = classical['zero_drift_heat_response_exact']/scale/S
        if not (0 < ratio <= 1+2e-10 and kappa > 0 and classical_ratio > 0):
            raise ValueError('nonpositive or nonvariational heat comparison')
        drag_checks = []
        for i, mode in enumerate(saved['modes']):
            prior = mode['cases'][j]
            if prior['screening'] != label or prior['screening_length_cm'] != length:
                raise ValueError('prescribed screening differs between saved modes')
            c = charges[i]*n[i]/ne
            weight = charge_weight[i]*(1-c)**2 + float(np.sum(np.delete(charge_weight, i)))*c*c
            recovered = unit_prefactor*I[0, 0]*weight
            error = abs(recovered/prior['electron_resistance_g_cm3_s']-1)
            if error > 1e-9:
                raise ValueError('heat calculation does not recover the saved electron drag')
            drag_checks.append(dict(isotope=mode['isotope'], relative_difference=error))
        cases.append(dict(screening=label, screening_length_cm=length,
                          b_thermal=B, dimensionless_collision_moments=I.tolist(),
                          dimensionless_full_energy_response=L.tolist(),
                          collision_prefactor_cgs=prefactor,
                          full_energy_kappa_erg_cm_s_K=kappa,
                          one_heat_variable_over_full_energy_kappa=ratio,
                          classical_over_fermi_full_energy_kappa=classical_ratio,
                          classical_over_fermi_displaced_drag=float(classical['collision_moments'][0][0]*scale/I[0, 0]),
                          recovered_drag_checks=drag_checks))
    return dict(zone=record['zone'], mass_fraction=record['mass_fraction'],
                enclosed_mass_fraction=record['enclosed_mass_fraction'],
                convective=record['convective'], temperature_K=T,
                density_baryonic_g_cm3=rho, electron_density_cm3=ne,
                eta_nonrelativistic=eta,
                temperature_over_fermi_temperature=record['temperature_over_fermi_temperature'],
                current_basis_density_relative_difference=density_error,
                current_variance=basis.variance,
                classical_density_normalization=scale, cases=cases,
                cpu_seconds=time.process_time()-started)


def summarize(rows):
    result = {}
    selections = dict(all_hot=rows,
                      nonconvective_hot=[r for r in rows if not r['convective']],
                      nonconvective_T_over_TF_below_0p3=[r for r in rows if not r['convective'] and r['temperature_over_fermi_temperature'] < .3])
    keys = ['one_heat_variable_over_full_energy_kappa',
            'classical_over_fermi_full_energy_kappa',
            'classical_over_fermi_displaced_drag', 'full_energy_kappa_erg_cm_s_K']
    for name, selected in selections.items():
        weight = np.array([r['mass_fraction'] for r in selected])
        item = dict(zones=len(selected), mass_fraction=float(np.sum(weight)), screening={})
        for j in range(2):
            s = {}
            for key in keys:
                values = np.array([r['cases'][j][key] for r in selected])
                s[key] = dict(minimum=float(np.min(values)), maximum=float(np.max(values)),
                              mass_weighted_mean=float(np.average(values, weights=weight)))
            item['screening'][selected[0]['cases'][j]['screening']] = s
        result[name] = item
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if args.output.exists() or not 1 <= args.workers <= 6:
        raise ValueError('new output and one to six workers required')
    paths = [Path(__file__), Path('scripts/electron_ion_heat.py'),
             Path('docs/results/electron_ion_heat_v2.json'),
             Path('docs/results/stellar_electron_drag_3890gyr_v1.json'),
             Path('docs/results/diffusion_material_regime_3890gyr_v1.json'),
             Path('include/ember/gs98_mixture.hpp')]
    inputs = {str(p.resolve()): digest(p) for p in paths}
    audit, drag, regime = [json.loads(p.read_text()) for p in paths[2:5]]
    if audit['outcome'] != 'passed' or not all(c['passed'] for c in audit['checks']):
        raise ValueError('electron heat checks have not passed')
    for document, key in [(audit, 'inputs_sha256'), (drag, 'input_sha256'), (regime, 'input_sha256')]:
        for name, h in document[key].items():
            if digest(name) != h:
                raise ValueError('saved calculation source changed: '+name)
            inputs[str(Path(name).resolve())] = h
    metals = [tuple(map(float, s.split(','))) for s in re.findall(r'^\s*\{([^{}]+)\}, //', paths[5].read_text(), re.M)]
    if len(metals) != 19 or len(regime['records']) != len(drag['records']):
        raise ValueError('mixture or profile dimensions changed')
    args.work.mkdir(parents=True, exist_ok=True)
    manifest = args.work/'inputs.json'
    if manifest.exists():
        if json.loads(manifest.read_text()) != inputs:
            raise ValueError('retained zone calculations have different inputs')
    else:
        write_new(manifest, inputs)
    rows, pending = [], []
    for record, saved in zip(regime['records'], drag['records'], strict=True):
        path = args.work/f"zone-{record['zone']:04d}.json"
        if path.exists():
            prior = json.loads(path.read_text())
            if prior['inputs_sha256'] != inputs or prior['result']['zone'] != record['zone']:
                raise ValueError('retained zone calculation differs')
            rows.append(prior['result'])
        else:
            pending.append((record, saved, path))
    reused = len(rows)
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context('spawn')) as pool:
        futures = {pool.submit(zone_comparison, record, saved, metals): path for record, saved, path in pending}
        for future in as_completed(futures):
            row = future.result()
            write_new(futures[future], dict(inputs_sha256=inputs, result=row))
            rows.append(row)
            if sum(p.stat().st_size for p in args.work.iterdir() if p.is_file()) > 10000000:
                raise ValueError('zone artifact cap exceeded')
    rows.sort(key=lambda r: r['zone'])
    summary = summarize(rows)
    if abs(summary['all_hot']['mass_fraction']/regime['selected_hot_mass_fraction']-1) > 1e-12:
        raise ValueError('summed saved-zone mass differs')
    result = dict(scope=__doc__, created_utc=datetime.now(timezone.utc).isoformat(),
                  outcome='completed_model_comparison', accepted_for_stellar_evolution=False,
                  hot_zones=len(rows), summary=summary, records=rows,
                  saved_zone_calculations_reused=reused,
                  classical_comparison=dict(shape_eta=CLASSICAL_ETA,
                                            density='renormalized to the actual ne at the same T and screening length',
                                            occupation_weight_relative_bound=2*math.exp(CLASSICAL_ETA)),
                  maximum_current_density_relative_difference=max(r['current_basis_density_relative_difference'] for r in rows),
                  maximum_retained_drag_relative_difference=max(c['relative_difference'] for r in rows for case in r['cases'] for c in case['recovered_drag_checks']),
                  limitations=['Fixed-ion nonrelativistic elastic Born scattering with fully stripped ions.',
                               'Two prescribed screening choices; no accepted correlation or dynamic-screening prescription.',
                               'No electron-electron collisions, recoil, moving-ion heat coupling or relativistic corrections.',
                               'This does not compare with or replace Ember conductivity and does not give stellar diffusion velocities.'],
                  inputs_sha256=inputs)
    for name, h in inputs.items():
        if digest(name) != h:
            raise ValueError('input changed during the calculation')
    write_new(args.output, result)
    print(json.dumps({k:result[k] for k in ['outcome','hot_zones','saved_zone_calculations_reused','summary','maximum_current_density_relative_difference','maximum_retained_drag_relative_difference']}), flush=True)


if __name__ == '__main__':
    main()
