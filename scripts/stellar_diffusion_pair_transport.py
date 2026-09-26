#!/usr/bin/env python3
"""Apply the checked pair interpolator to coupled transport in the saved hot star.

Retains the EOS forces and exact elastic comparison. No pair integrals or EOS
queries are repeated. Coefficients remain conditional on the specified kinetic
model; the output is not a conservative abundance or structure update.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import time

import numpy as np

from audit_electron_pair_interpolation import relative_operator_error
from diffusion_burgers import KB, MU
from electron_electron_energy_modes_v2 import PolynomialBasis
from electron_mixture_collision import reduce_electron_collision, symmetric
from electron_pair_table import ElectronPairTable
from ion_collision_integrals import E2
from ion_collision_table import ExtendedIonTable
from ion_electron_transport import coupled_elastic_operator, heat_decomposition


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    assert not args.output.exists()
    root = Path(__file__).resolve().parents[1]
    inputs = {}
    def read(name):
        path = root/'docs/results'/name
        raw = path.read_bytes()
        inputs[str(path)] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)
    acceptance = read('electron_pair_interpolation_v2.json')
    assert acceptance['outcome'] == 'passed_independent_controls'
    data = read('electron_pair_table_sources_v2.json')
    forces = read('stellar_diffusion_flux_comparison_3890gyr_v1.json')
    material = read('diffusion_material_regime_3890gyr_v1.json')
    heat = read('stellar_electron_heat_3890gyr_v1.json')
    old_mixture = read('electron_mixture_collision_v1.json')
    extension = read('ion_collision_table_extension_v1.json')
    assert hashlib.sha256(Path(extension['table']).read_bytes()).hexdigest() == extension['table_sha256']
    inputs[extension['table']] = extension['table_sha256']
    ion_table = ExtendedIonTable(extension['table'])
    pair_table = ElectronPairTable(root/'docs/results/electron_pair_table_sources_v2.json')
    assert pair_table.report_sha256 == acceptance['input_sha256'][str(root/'docs/results/electron_pair_table_sources_v2.json')]
    metal_path = root/'include/ember/gs98_mixture.hpp'
    inputs[str(metal_path)] = hashlib.sha256(metal_path.read_bytes()).hexdigest()
    metals = [tuple(map(float, s.split(','))) for s in re.findall(r'^\s*\{([^{}]+)\}, //', metal_path.read_text(), re.M)]
    mass = np.array([1., 3., 4.]+[m[1] for m in metals])
    charge = np.array([1., 2., 2.]+[m[0] for m in metals])
    reduced_mass = MU*mass[:, None]*mass[None, :]/(mass[:, None]+mass[None, :])
    charge_product = np.outer(charge, charge)*E2
    for name in ['stellar_diffusion_pair_transport.py', 'electron_pair_table.py',
                 'electron_electron_energy_modes_v2.py', 'electron_mixture_collision.py',
                 'ion_electron_transport.py', 'ion_collision_table.py', 'electron_elastic_mixture.py']:
        path = root/'scripts'/name
        inputs[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    controls = {r['zone']: r for r in old_mixture['records']}
    maxima = dict(baryon_flux=0., electric_current=0., entropy_balance=0., heat_identity=0.,
                  negative_entropy=0., old_elastic_heat_reproduction=0., retained_mixture_response_error=0.,
                  retained_mixture_conductivity_difference=0., eight_to_ten_heat_difference=0.)
    def check(name, error):
        if not np.isfinite(error):
            raise ValueError('nonfinite check: '+name)
        maxima[name] = max(maxima[name], abs(float(error)))
    records = []
    start = time.process_time()
    for mr, hr, fr in zip(material['records'], heat['records'], forces['records'], strict=True):
        zone = mr['zone']
        assert hr['zone'] == fr['zone'] == zone
        t, rho = mr['temperature_K'], mr['density_baryonic_g_cm3']
        assert (t, rho) == (fr['temperature_K'], fr['density_g_cm3'])
        ztotal = 1-mr['X']-mr['Y3']-mr['Y4']
        x = np.array([mr['X'], mr['Y3'], mr['Y4']]+[ztotal*m[3] for m in metals])
        ni = rho/MU*x/mass
        ne = float(ni@charge)
        basis = PolynomialBasis(hr['eta_nonrelativistic'], 9)
        pair_matrices = pair_table.matrices(hr['eta_nonrelativistic'], [c['b_thermal'] for c in hr['cases']])
        cases = []
        for sample, ee in zip(hr['cases'], pair_matrices, strict=True):
            label = sample['screening']
            old = next(c for c in fr['cases'] if c['screening'] == label)
            moments = ion_table.moments(charge_product/(KB*t*sample['screening_length_cm']))
            omega11 = np.sqrt(2*math.pi/reduced_mass)*charge_product**2/(KB*t)**1.5*moments[:, :, 0]
            ki = 16/3*np.outer(ni, ni)*reduced_mass*omega11
            ei = symmetric(basis.ion_matrix(sample['b_thermal']))
            reduced = reduce_electron_collision(electron_ion_matrix=ei, electron_electron_matrix=ee,
                                                ion_density=ni, ion_charges=charge)
            common = dict(density=rho, temperature=t, mass_fractions=x, mass_numbers=mass,
                charges=charge, independent_ions=[0, 1], reference_ion=2, ion_resistance=ki,
                ion_z=1-.4*moments[:, :, 1]/moments[:, :, 0],
                ion_zprime=2.5-2*moments[:, :, 1]/moments[:, :, 0]+.4*moments[:, :, 2]/moments[:, :, 0],
                ion_zdoubleprime=moments[:, :, 3]/moments[:, :, 0], electron_moments=ei[:2, :2],
                energy_scale=KB/MU*t)
            op = coupled_elastic_operator(**common, electron_full_response=reduced.retained_response)
            htransport, kappa = heat_decomposition(op.mobility, temperature=t, energy_scale=op.energy_scale)
            reduced8 = reduce_electron_collision(electron_ion_matrix=ei[:8, :8], electron_electron_matrix=ee[:8, :8],
                                                 ion_density=ni, ion_charges=charge)
            op8 = coupled_elastic_operator(**common, electron_full_response=reduced8.retained_response)
            _, k8 = heat_decomposition(op8.mobility, temperature=t, energy_scale=op8.energy_scale)
            check('eight_to_ten_heat_difference', kappa/k8-1)
            if zone in controls and label == controls[zone]['screening']:
                control = controls[zone]
                check('retained_mixture_response_error', relative_operator_error(reduced.retained_response,
                      np.asarray(control['retained_electron_response'])))
                check('retained_mixture_conductivity_difference', kappa/control['conductivity_zero_element_flux_cgs']-1)
                exact = coupled_elastic_operator(**{**common, 'electron_moments': sample['dimensionless_collision_moments']},
                                                electron_full_response=sample['dimensionless_full_energy_response'])
                _, exact_k = heat_decomposition(exact.mobility, temperature=t, energy_scale=exact.energy_scale)
                check('old_elastic_heat_reproduction', exact_k/old['conductivity_zero_element_flux_cgs']-1)
            solutions = {}
            for stride in ('1', '2'):
                gradient = fr['gradients'][stride]
                for forcing, key in [('saved_density', 'chemical_acceleration_cm_s2'),
                                     ('hydrostatic_density', 'hydrostatic_chemical_acceleration_cm_s2')]:
                    force = np.r_[-np.asarray(gradient[key])/t, -op.energy_scale*gradient['dlnT_dr']/t]
                    sol = op.solve(force)
                    noheat = op.solve(np.r_[force[:-1], 0.])
                    speed = max(float(np.max(abs(sol['velocity']))), np.finfo(float).tiny)
                    entropy_scale = max(abs(sol['entropy_from_forces']), np.finfo(float).tiny)
                    check('baryon_flux', x@sol['velocity'][:-1]/speed)
                    check('electric_current', (ni*charge@sol['velocity'][:-1]-ne*sol['velocity'][-1])/(2*ne*speed))
                    check('entropy_balance', (sol['entropy_from_forces']-sol['entropy_from_collisions'])/entropy_scale)
                    check('negative_entropy', min(sol['entropy_from_forces'], 0.)/entropy_scale)
                    reconstructed = htransport@sol['mass_flux']-kappa*t*gradient['dlnT_dr']
                    check('heat_identity', (reconstructed-sol['reduced_heat_flux'])/
                          max(abs(reconstructed), abs(sol['reduced_heat_flux']), 1.))
                    solutions[f'{stride}:{forcing}'] = dict(velocity_cm_s=sol['velocity'][:3].tolist(),
                        thermal_velocity_contribution_cm_s=(sol['velocity']-noheat['velocity'])[:3].tolist(),
                        independent_mass_flux_g_cm2_s=sol['mass_flux'].tolist(),
                        reduced_heat_flux_erg_cm2_s=sol['reduced_heat_flux'],
                        full_material_energy_flux_erg_cm2_s=float(sol['reduced_heat_flux']+
                            np.asarray(fr['exchange_enthalpy_erg_g'])@sol['mass_flux']),
                        pressure_scale_crossing_years=[float(fr['pressure_scale_height_cm']/abs(v)/(365.25*86400))
                                                       if v else None for v in sol['velocity'][:3]])
            cases.append(dict(screening=label, electron_electron_prefactor_ratio=reduced.electron_electron_prefactor_ratio,
                mobility_scaled_heat=op.mobility.tolist(), energy_scale_erg_g=op.energy_scale,
                conductivity_zero_element_flux_cgs=kappa,
                conductivity_over_exact_elastic=kappa/old['conductivity_zero_element_flux_cgs'],
                transport_enthalpy_erg_g=htransport.tolist(), solutions=solutions))
        records.append(dict(zone=zone, mass_fraction=fr['mass_fraction'], enclosed_mass_fraction=fr['enclosed_mass_fraction'],
            convective=fr['convective'], temperature_K=t, density_g_cm3=rho,
            exchange_enthalpy_erg_g=fr['exchange_enthalpy_erg_g'], cases=cases))
    summaries = {}
    for label in [c['screening'] for c in records[0]['cases']]:
        rows = [r for r in records if not r['convective']]
        weights = np.array([r['mass_fraction'] for r in rows])
        weights /= sum(weights)
        cases = [next(c for c in r['cases'] if c['screening'] == label) for r in rows]
        v = np.array([c['solutions']['1:saved_density']['velocity_cm_s'] for c in cases])
        thermal = np.array([c['solutions']['1:saved_density']['thermal_velocity_contribution_cm_s'] for c in cases])
        old_cases = [next(c for c in forces['records'][r['zone']]['cases'] if c['screening'] == label) for r in rows]
        old_v = np.array([c['solutions']['1:saved_density']['velocity_cm_s'] for c in old_cases])
        ratios = [c['conductivity_over_exact_elastic'] for c in cases]
        summaries[label] = dict(nonconvective_hot_layers=len(rows),
            conductivity_over_exact_elastic_minimum=min(ratios), conductivity_over_exact_elastic_maximum=max(ratios),
            mass_weighted_absolute_velocity_cm_s=(weights@abs(v)).tolist(),
            mass_weighted_absolute_velocity_change_over_elastic=(weights@abs(v-old_v)/(weights@abs(old_v))).tolist(),
            mass_weighted_absolute_thermal_over_nonthermal=(weights@abs(thermal)/(weights@abs(v-thermal))).tolist())
    tolerances = {key: 2e-10 for key in maxima}
    for key in ['retained_mixture_response_error', 'retained_mixture_conductivity_difference', 'eight_to_ten_heat_difference']:
        tolerances[key] = .005
    passed = all(maxima[k] <= tolerances[k] for k in maxima)
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(),
        outcome='completed_conditional_pair_transport' if passed else 'failed_numerical_checks',
        accepted_for_stellar_evolution=False, records=records, summaries=summaries,
        maximum_check_errors=maxima, check_tolerances=tolerances,
        input_sha256=inputs, new_EOS_queries=0, new_pair_integrals=0,
        new_electron_ion_moment_queries=2*len(records), cpu_seconds=time.process_time()-start,
        limitations=['Fully stripped stationary metals, common effective static screening, classical ions and nonrelativistic Born/Pauli electrons.',
                     'Mixture correlations, changing ionization and matching kinetic forces to the nonideal EOS remain conditional.',
                     'Pressure-scale drift times are local diagnostics, not stellar lifetimes.',
                     'No conservative abundance, convective coupling, heat update or stellar structure evolution.'])
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({k: result[k] for k in ['outcome', 'summaries', 'maximum_check_errors', 'cpu_seconds']}, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
