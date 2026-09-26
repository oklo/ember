#!/usr/bin/env python3
"""Apply a conditional coupled transport model to the saved stellar gradients.

All EOS replies, electron integrals and ion scattering sources are reused.
This does not evolve composition or select a conductivity. Electron-electron
collisions, correlation corrections and changing ionization remain absent.
"""
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np
from audit_smooth_eos_runtime import ARAD
from diffusion_burgers import KB, MU
from diffusion_composition_forces import chemical_acceleration
from diffusion_thermal_transport import exchange_enthalpy, HeatFluxConvention
from ion_collision_integrals import E2
from ion_collision_table import ExtendedIonTable
from ion_electron_transport import coupled_elastic_operator, heat_decomposition


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def derivative(radius, values, index, stride):
    """Local quadratic derivative, using an even central reflection if needed."""
    if index < stride:
        xs = np.array([-radius[index], radius[index], radius[index+stride]])
        ys = np.array([values[index], values[index], values[index+stride]])
    else:
        xs = radius[[index-stride, index, index+stride]]
        ys = values[[index-stride, index, index+stride]]
    t = xs-radius[index]
    weights = np.array([-(t[1]+t[2])/((t[0]-t[1])*(t[0]-t[2])),
                        -(t[0]+t[2])/((t[1]-t[0])*(t[1]-t[2])),
                        -(t[0]+t[1])/((t[2]-t[0])*(t[2]-t[1]))])
    # Constant offsets cancel analytically rather than through large products.
    return weights @ (ys-ys[1])


def main():
    inputs = {}
    def read(path):
        path = Path(path); inputs[str(path)] = sha(path)
        return json.loads(path.read_text())
    retained_path = Path('docs/results/smooth_eos_refined_retained_sources_v1.json')
    retained = read(retained_path)
    assert retained['outcome'] == 'source_value_comparisons_passed'
    cached_input = Path('/tmp/ember-smooth-eos-refined-physical-run-v1/retained/retained.input')
    cached_output = cached_input.with_suffix('.output')
    for p in [cached_input, cached_output]:
        assert sha(p) == retained['artifacts_sha256'][str(p)]
        inputs[str(p)] = sha(p)
    cache = {}
    for text, line in zip(cached_input.read_text().splitlines(),
                          cached_output.read_text().splitlines(), strict=True):
        key = tuple(map(float, text.split())); reply = json.loads(line)
        if reply['ok'] and key[-1] == 1:
            cache[key[:4]] = np.array(reply['values'])
    profile = Path('docs/reports/2026-09-11/evolution_latest_profile.csv')
    inputs[str(profile)] = sha(profile)
    data = list(csv.DictReader(profile.open()))
    vectors = {k: np.array([float(r[k]) for r in data]) for k in data[0]}
    radius, rho, T = [vectors[k] for k in ['radius_cm','density_g_cm3','temperature_K']]
    composition = np.c_[vectors['X'], vectors['Y3']]
    eos = np.array([cache[tuple(float(r[k]) for k in ['X','Y3','temperature_K','density_g_cm3'])]
                    for r in data])
    heat_path = Path('docs/results/stellar_electron_heat_3890gyr_v1.json')
    heat = read(heat_path)
    regime = read('docs/results/diffusion_material_regime_3890gyr_v1.json')
    extension = read('docs/results/ion_collision_table_extension_v1.json')
    previous = read('docs/results/ion_electron_transport_v2.json')
    assert previous['outcome'] == 'passed' and extension['status'] == 'pass'
    assert sha(extension['table']) == extension['table_sha256']
    table = ExtendedIonTable(extension['table']); inputs[extension['table']] = sha(extension['table'])
    for report in [heat, regime, extension, previous]:
        for key in ['input_sha256', 'inputs_sha256']:
            for path, digest in report.get(key, {}).items():
                assert sha(path) == digest, path
                inputs[path] = digest
    metals_path = Path('include/ember/gs98_mixture.hpp')
    metals = [tuple(map(float,s.split(','))) for s in re.findall(
        r'^\s*\{([^{}]+)\}, //', metals_path.read_text(), re.M)]
    mass = np.array([1.,3.,4.]+[m[1] for m in metals])
    charge = np.array([1.,2.,2.]+[m[0] for m in metals])
    constants = Path('include/ember/constants.hpp')
    G = float(re.search(r'\bG\s*=\s*([\d.eE+-]+)',constants.read_text()).group(1))
    seconds_per_year = 365.25*86400
    paths = [Path(__file__), constants, metals_path, Path('scripts/audit_smooth_eos_runtime.py')]
    inputs.update({str(p):sha(p) for p in paths})
    checks = []
    def check(name, error, tolerance=2e-10, **meta):
        checks.append(dict(name=name,error=float(error),tolerance=tolerance,
                           passed=bool(np.isfinite(error) and abs(error)<=tolerance),**meta))
    rr = np.geomspace(.01,2.,101)
    for stride in [1,2]:
        for i in [0,1,12,50,97]:
            got = derivative(rr, 3+2*rr**2, i, stride)
            check('even quadratic gradient',abs(got-4*rr[i])/(4*rr[i]),1e-10,index=i,stride=stride)
    records = []
    previous_cases = {(r['zone'],r['screening']):r for r in previous['physical_cases']}
    for hr, material in zip(heat['records'], regime['records'], strict=True):
        i = hr['zone']; assert i == material['zone'] and i+2 < len(radius)
        v = eos[i]; temp = T[i]; density = rho[i]
        prho = v[0]*v[6]; prad = ARAD*temp**4/3
        delta = (v[0]*v[5]-4*prad)/prho
        hessian = v[24:28].reshape(2,2); gt=v[28:30]; gr=v[30:32]
        enthalpy = exchange_enthalpy(temperature=temp,gradient_phi_lnT=gt,
                                    gradient_phi_lnrho=gr,material_delta=delta)
        check('physical enthalpy identities',np.max(abs(enthalpy-(v[19:21]+delta*v[17:19]/density)))/
              (KB/MU*temp+np.max(abs(enthalpy))),zone=i)
        gravity = G*vectors['mass_g'][i]/radius[i]**2
        hp = v[0]/(density*gravity)
        gradients = {}
        for stride in [1,2]:
            dlnt = derivative(radius,np.log(T),i,stride)
            dlnrho = derivative(radius,np.log(rho),i,stride)
            dc = derivative(radius,composition,i,stride)
            pressure_gradient = prho*dlnrho+v[0]*v[5]*dlnt+v[17:19]@dc
            # Alternate density gradient imposes hydrostatic balance with the
            # actual enclosed mass and the same EOS. Retain both explicitly.
            hydro_dlnrho = (-density*gravity-v[0]*v[5]*dlnt-v[17:19]@dc)/prho
            accel = chemical_acceleration(temperature=temp,hessian_phi=hessian,
                gradient_phi_lnrho=gr,material_delta=delta,composition_gradient=dc,
                logarithmic_density_gradient=dlnrho,logarithmic_temperature_gradient=dlnt)
            hydro_accel = chemical_acceleration(temperature=temp,hessian_phi=hessian,
                gradient_phi_lnrho=gr,material_delta=delta,composition_gradient=dc,
                logarithmic_density_gradient=hydro_dlnrho,logarithmic_temperature_gradient=dlnt)
            gradients[str(stride)] = dict(dlnT_dr=float(dlnt),dlnrho_dr=float(dlnrho),
                dcomposition_dr=dc.tolist(),chemical_acceleration_cm_s2=accel.tolist(),
                hydrostatic_chemical_acceleration_cm_s2=hydro_accel.tolist(),
                relative_hydrostatic_residual=float((pressure_gradient+density*gravity)/(density*gravity)))
        total_metals = 1-material['X']-material['Y3']-material['Y4']
        x = np.array([material['X'],material['Y3'],material['Y4']]+[total_metals*m[3] for m in metals])
        ni = density/MU*x/mass; ne=float(ni@charge)
        cases = []
        for sample in hr['cases']:
            label = sample['screening']; lam=sample['screening_length_cm']
            strength = np.outer(charge,charge)*E2/(KB*temp*lam)
            moments = table.moments(strength)
            reduced = MU*(mass[:,None]*mass[None,:]/(mass[:,None]+mass[None,:]))
            omega11 = np.sqrt(2*math.pi/reduced)*(np.outer(charge,charge)*E2)**2/(KB*temp)**1.5*moments[:,:,0]
            ki = 16/3*np.outer(ni,ni)*reduced*omega11
            op = coupled_elastic_operator(density=density,temperature=temp,mass_fractions=x,
                mass_numbers=mass,charges=charge,independent_ions=[0,1],reference_ion=2,
                ion_resistance=ki,ion_z=1-.4*moments[:,:,1]/moments[:,:,0],
                ion_zprime=2.5-2*moments[:,:,1]/moments[:,:,0]+.4*moments[:,:,2]/moments[:,:,0],
                ion_zdoubleprime=moments[:,:,3]/moments[:,:,0],
                electron_moments=sample['dimensionless_collision_moments'],
                electron_full_response=sample['dimensionless_full_energy_response'],energy_scale=KB/MU*temp)
            htransport,kappa=heat_decomposition(op.mobility,temperature=temp,energy_scale=op.energy_scale)
            if (i,label) in previous_cases:
                check('retained conductivity comparison',kappa/previous_cases[i,label]['conductivity_zero_element_flux_cgs']-1,
                      zone=i,screening=label)
            solutions = {}
            for stride in [1,2]:
                d=gradients[str(stride)];dlnt=d['dlnT_dr']
                for forcing in ['saved_density','hydrostatic_density']:
                    name='chemical_acceleration_cm_s2' if forcing=='saved_density' else 'hydrostatic_chemical_acceleration_cm_s2'
                    f=np.r_[-np.array(d[name])/temp,-op.energy_scale*dlnt/temp]
                    sol=op.solve(f); noheat=op.solve(np.r_[f[:-1],0.])
                    speed=max(abs(sol['velocity']));normal=max(abs(sol['entropy_from_forces']),np.finfo(float).tiny)
                    meta=dict(zone=i,screening=label,stride=stride,forcing=forcing)
                    check('baryonic mass conservation',x@sol['velocity'][:-1]/speed,1e-13,**meta)
                    check('electric current conservation',(ni*charge@sol['velocity'][:-1]-ne*sol['velocity'][-1])/(2*ne*speed),1e-13,**meta)
                    check('entropy balance',(sol['entropy_from_forces']-sol['entropy_from_collisions'])/normal,**meta)
                    check('positive dissipation',min(sol['entropy_from_forces'],0.)/normal,**meta)
                    qprime=htransport@sol['mass_flux']-kappa*temp*dlnt
                    check('direct heat reconstruction',(qprime-sol['reduced_heat_flux'])/max(abs(qprime),abs(sol['reduced_heat_flux']),1.),**meta)
                    convention=HeatFluxConvention(enthalpy,op.energy_scale)
                    if forcing=='saved_density':
                        direct=-hessian@np.array(d['dcomposition_dr'])-gr*d['dlnrho_dr']-gt*dlnt
                        full=convention.full_forces(f)
                        check('matched EOS thermal force',np.max(abs(full[:-1]-direct))/max(np.max(abs(direct)),1e-300),**meta)
                    velocity=sol['velocity'][[0,1,2]];thermal=velocity-noheat['velocity'][[0,1,2]]
                    solutions[f'{stride}:{forcing}']=dict(velocity_cm_s=velocity.tolist(),
                        velocity_without_thermal_force_cm_s=noheat['velocity'][[0,1,2]].tolist(),
                        thermal_velocity_contribution_cm_s=thermal.tolist(),
                        independent_mass_flux_g_cm2_s=sol['mass_flux'].tolist(),
                        reduced_heat_flux_erg_cm2_s=sol['reduced_heat_flux'],
                        full_material_energy_flux_erg_cm2_s=float(sol['reduced_heat_flux']+enthalpy@sol['mass_flux']),
                        entropy_production_erg_cm3_s_K=sol['entropy_from_forces'],
                        pressure_scale_crossing_years=[float(hp/abs(w)/seconds_per_year) if w else None for w in velocity])
            cases.append(dict(screening=label,energy_scale_erg_g=op.energy_scale,
                mobility_scaled_heat=op.mobility.tolist(),conductivity_zero_element_flux_cgs=kappa,
                transport_enthalpy_erg_g=htransport.tolist(),solutions=solutions))
        records.append(dict(zone=i,mass_fraction=hr['mass_fraction'],convective=hr['convective'],
            enclosed_mass_fraction=hr['enclosed_mass_fraction'],temperature_K=float(temp),density_g_cm3=float(density),
            gravity_cm_s2=float(gravity),pressure_scale_height_cm=float(hp),
            exchange_enthalpy_erg_g=enthalpy.tolist(),gradients=gradients,cases=cases))
    summaries={}
    for screening in [s['screening'] for s in records[0]['cases']]:
        rows=[r for r in records if not r['convective']]
        weights=np.array([r['mass_fraction'] for r in rows]);weights/=sum(weights)
        def values(key,quantity):
            return np.array([next(c for c in r['cases'] if c['screening']==screening)['solutions'][key][quantity] for r in rows])
        fine=values('1:saved_density','velocity_cm_s');coarse=values('2:saved_density','velocity_cm_s')
        hydro=values('1:hydrostatic_density','velocity_cm_s');thermal=values('1:saved_density','thermal_velocity_contribution_cm_s')
        base=values('1:saved_density','velocity_without_thermal_force_cm_s')
        summaries[screening]=dict(nonconvective_hot_layers=len(rows),
            stellar_mass_fraction=float(sum(r['mass_fraction'] for r in rows)),
            isotope_order=['H1','He3','He4'],
            minimum_velocity_cm_s=np.min(fine,axis=0).tolist(),maximum_velocity_cm_s=np.max(fine,axis=0).tolist(),
            mass_weighted_mean_absolute_velocity_cm_s=(weights@abs(fine)).tolist(),
            mass_weighted_absolute_thermal_over_nonthermal=(weights@abs(thermal)/(weights@abs(base))).tolist(),
            mass_fraction_of_hot_nonconvective_region_with_thermal_over_ten_percent=(weights@(abs(thermal)>.1*abs(base))).tolist(),
            mass_fraction_of_hot_nonconvective_region_with_direction_reversal=(weights@(fine*base<0)).tolist(),
            sampling_change_relative_to_mass_weighted_absolute_velocity=(weights@abs(fine-coarse)/(weights@abs(fine))).tolist(),
            hydrostatic_change_relative_to_mass_weighted_absolute_velocity=(weights@abs(fine-hydro)/(weights@abs(fine))).tolist())
    assert all(sha(p)==h for p,h in inputs.items())
    passed=all(c['passed'] for c in checks)
    report=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),
        outcome='completed_conditional_comparison' if passed else 'failed_numerical_checks',
        accepted_for_stellar_evolution=False,new_EOS_queries=0,new_scattering_calculations=0,
        cached_stellar_chemical_replies_reused=len(eos),saved_electron_layer_screen_pairs_reused=844,
        checks=checks,summaries=summaries,records=records,inputs_sha256=inputs,
        limitations=['No electron-electron collisions, higher-order recoil or effective-potential mixture correlations.',
            'Fully stripped stationary metals and prescribed common Born screening.',
            'Local velocities are conditional on matching the kinetic reduced heat variable to the EOS partial-enthalpy convention.',
            'No conservative abundance update, convective coupling or re-equilibrated stellar structure.',
            'Pressure-scale crossing times are local drift diagnostics, not remaining stellar lifetimes.',
            'Two radial stencils and hydrostatic reconstruction test one saved profile, not evolutionary mesh convergence.'])
    target=Path('docs/results/stellar_diffusion_flux_comparison_3890gyr_v1.json')
    with target.open('x') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(outcome=report['outcome'],checks=len(checks),
        failures=[c for c in checks if not c['passed']],summaries=summaries)))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
