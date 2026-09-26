#!/usr/bin/env python3
"""Apply retained electron pair integrals to three actual stellar mixtures.

Independently assemble all ion/electron energy modes and compare their solve
with the reduced drift/heat matrix. No EOS or collision integrals are repeated.
"""
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np

from diffusion_burgers import KB,MU
from electron_ion_heat import physical_prefactor
from electron_mixture_collision import reduce_electron_collision,positive_solve,symmetric
from ion_collision_integrals import E2
from ion_collision_table import ExtendedIonTable
from ion_electron_transport import coupled_elastic_operator,heat_decomposition

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT/'docs/results/electron_mixture_collision_v1.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    assert not OUTPUT.exists()
    inputs = {}
    def read(name,expected=None):
        path=ROOT/'docs/results'/name
        sha=digest(path)
        if expected is not None:assert sha==expected
        inputs[str(path)]=sha
        return json.loads(path.read_text())
    energy=read('electron_electron_energy_refinement_v1.json','11414f3aa812fe07300d0d1bb527efc89c8187dcb04841588109db9bf4248b4d')
    forces=read('stellar_diffusion_flux_comparison_3890gyr_v1.json','19c39a5276330de8dd3a70b3345cef286efe2193135ecef7930c2899aef833cd')
    heat=read('stellar_electron_heat_3890gyr_v1.json')
    material=read('diffusion_material_regime_3890gyr_v1.json')
    extension=read('ion_collision_table_extension_v1.json')
    assert digest(extension['table'])==extension['table_sha256']
    inputs[extension['table']]=digest(extension['table'])
    table=ExtendedIonTable(extension['table'])
    assert energy['outcome']=='passed' and forces['outcome']=='completed_conditional_comparison'
    for name in ['electron_mixture_collision.py','ion_electron_transport.py','electron_elastic_mixture.py',
                 'electron_ion_heat.py','ion_collision_table.py','diffusion_composition_forces.py',
                 'ion_collision_integrals.py','diffusion_burgers.py']:
        p=ROOT/'scripts'/name;inputs[str(p)]=digest(p)
    inputs[str(Path(__file__))]=digest(__file__)
    path=ROOT/'include/ember/gs98_mixture.hpp';inputs[str(path)]=digest(path)
    metals=[tuple(map(float,s.split(','))) for s in re.findall(r'^\s*\{([^{}]+)\}, //',path.read_text(),re.M)]
    mass=np.array([1.,3.,4.]+[m[1] for m in metals])
    charge=np.array([1.,2.,2.]+[m[0] for m in metals])
    checks=[];records=[]
    def check(name,error,tolerance=2e-10,**meta):
        checks.append(dict(name=name,error=float(error),tolerance=tolerance,
                           passed=bool(np.isfinite(error) and abs(error)<=tolerance),**meta))
    def matrix_error(a,b):
        scale=np.sqrt(np.diag(b))
        return float(np.max(abs(a-b)/scale[:,None]/scale[None,:]))
    for zone in [0,300,396]:
        mr,hr,fr=material['records'][zone],heat['records'][zone],forces['records'][zone]
        assert mr['zone']==hr['zone']==fr['zone']==zone
        calc=energy['calculations'][f'zone_{zone}']
        sample=next(c for c in hr['cases'] if c['screening']=='total_effective_static_screening')
        previous=next(c for c in fr['cases'] if c['screening']==sample['screening'])
        assert calc['eta']==hr['eta_nonrelativistic'] and calc['b_thermal']==sample['b_thermal']
        t,rho=mr['temperature_K'],mr['density_baryonic_g_cm3']
        assert (t,rho)==(fr['temperature_K'],fr['density_g_cm3'])
        ztotal=1-mr['X']-mr['Y3']-mr['Y4']
        fractions=np.array([mr['X'],mr['Y3'],mr['Y4']]+[ztotal*m[3] for m in metals])
        ni=rho/MU*fractions/mass;ne=float(ni@charge)
        strength=np.outer(charge,charge)*E2/(KB*t*sample['screening_length_cm'])
        moments=table.moments(strength)
        reduced_mass=MU*(mass[:,None]*mass[None,:]/(mass[:,None]+mass[None,:]))
        omega11=np.sqrt(2*math.pi/reduced_mass)*(np.outer(charge,charge)*E2)**2/(KB*t)**1.5*moments[:,:,0]
        ki=16/3*np.outer(ni,ni)*reduced_mass*omega11
        ei=symmetric(calc['electron_ion_matrix']);ee=symmetric(calc['collision_matrix'])
        common=dict(density=rho,temperature=t,mass_fractions=fractions,mass_numbers=mass,charges=charge,
                    independent_ions=[0,1],reference_ion=2,ion_resistance=ki,
                    ion_z=1-.4*moments[:,:,1]/moments[:,:,0],
                    ion_zprime=2.5-2*moments[:,:,1]/moments[:,:,0]+.4*moments[:,:,2]/moments[:,:,0],
                    ion_zdoubleprime=moments[:,:,3]/moments[:,:,0],electron_moments=ei[:2,:2],
                    energy_scale=KB/MU*t)
        history=[];operators={}
        for count in [2,4,6,8,9,10]:
            reduced=reduce_electron_collision(electron_ion_matrix=ei[:count,:count],
                       electron_electron_matrix=ee[:count,:count],ion_density=ni,ion_charges=charge)
            op=coupled_elastic_operator(**common,electron_full_response=reduced.retained_response)
            h,kappa=heat_decomposition(op.mobility,temperature=t,energy_scale=op.energy_scale)
            history.append(dict(modes=count,conductivity_zero_element_flux_cgs=kappa,
                                transport_enthalpy_erg_g=h.tolist(),mobility=op.mobility.tolist()))
            operators[count]=(op,reduced)
        op,reduced=operators[10]
        check('eight-to-ten-mode conductivity',history[-1]['conductivity_zero_element_flux_cgs']/history[-3]['conductivity_zero_element_flux_cgs']-1,.005,zone=zone)
        for first,second in zip(history,history[1:]):
            check('more energy modes increase conductivity',min(second['conductivity_zero_element_flux_cgs']/first['conductivity_zero_element_flux_cgs']-1,0),zone=zone,modes=second['modes'])
        check('electron-electron prefactor uses actual mixture',reduced.electron_electron_prefactor_ratio/(ne/float(ni@(charge*charge)))-1,zone=zone)
        # Independent assembly sums the full electron-ion collision cost for
        # each ion velocity. It does not use the reduced mean/variance split.
        nc=op.velocity_map.shape[1];count=len(mass);retained_count=len(op.collision_matrix)
        size=retained_count+8
        full=np.zeros((size,size));full[:retained_count,:retained_count]=op.ion_collision_matrix
        full[nc+np.arange(count),nc+np.arange(count)] += op.brownian_heat_diagonal
        for density,z,b in zip(ni,charge,op.velocity_map[:-1],strict=True):
            mapping=np.zeros((10,size));mapping[0,:nc]=op.velocity_map[-1]-b
            mapping[1,nc+count]=1;mapping[2:,retained_count:]=np.eye(8)
            full += physical_prefactor(float(density),float(z))*mapping.T@ei@mapping
        mapping=np.zeros((10,size));mapping[1,nc+count]=1;mapping[2:,retained_count:]=np.eye(8)
        full += physical_prefactor(ne,1.)*mapping.T@ee@mapping
        schur=symmetric(full[:retained_count,:retained_count]-full[:retained_count,retained_count:]@
                        positive_solve(full[retained_count:,retained_count:],full[retained_count:,:retained_count]))
        check('full independent assembly versus reduction',matrix_error(schur,op.collision_matrix),zone=zone)
        g=np.zeros((len(op.mobility),size));g[:,:retained_count]=op.flux_map
        full_response=t*positive_solve(full,g.T)
        full_mobility=g@full_response
        check('full solve versus reduced mobility',matrix_error(full_mobility,op.mobility),zone=zone)
        off_reduced=reduce_electron_collision(electron_ion_matrix=ei,electron_electron_matrix=ee,
                         ion_density=ni,ion_charges=charge,electron_electron_scale=0.)
        off=coupled_elastic_operator(**common,electron_full_response=off_reduced.retained_response)
        exact=coupled_elastic_operator(**{**common,'electron_moments':sample['dimensionless_collision_moments']},
                          electron_full_response=sample['dimensionless_full_energy_response'])
        _,koff=heat_decomposition(off.mobility,temperature=t,energy_scale=off.energy_scale)
        _,kexact=heat_decomposition(exact.mobility,temperature=t,energy_scale=exact.energy_scale)
        h,kappa=heat_decomposition(op.mobility,temperature=t,energy_scale=op.energy_scale)
        check('retained no-pair conductivity',kexact/previous['conductivity_zero_element_flux_cgs']-1,zone=zone)
        check('pair collisions reduce mobility',min(np.linalg.eigvalsh((off.mobility-op.mobility)/
              np.sqrt(np.diag(off.mobility))[:,None]/np.sqrt(np.diag(off.mobility))[None,:]).min(),0),zone=zone)
        check('finite energy space below exact elastic conductivity',max(koff/kexact-1,0),zone=zone)
        # Rephasing eliminated modes must change neither retained current.
        signs=np.diag([1.,1.]+[-1. if i%2 else 1. for i in range(8)])
        rephased=reduce_electron_collision(electron_ion_matrix=signs@ei@signs,electron_electron_matrix=signs@ee@signs,
                                         ion_density=ni,ion_charges=charge)
        check('independent higher-mode sign convention',matrix_error(rephased.retained_response,reduced.retained_response),zone=zone)
        solutions={}
        for stride in ['1','2']:
            gradient=fr['gradients'][stride]
            for forcing,key in [('saved_density','chemical_acceleration_cm_s2'),('hydrostatic_density','hydrostatic_chemical_acceleration_cm_s2')]:
                force=np.r_[-np.asarray(gradient[key])/t,-op.energy_scale*gradient['dlnT_dr']/t]
                sol=op.solve(force);base=off.solve(force);noheat=op.solve(np.r_[force[:-1],0.])
                variables=full_response@force;full_flux=g@variables
                scale=max(np.max(abs(sol['velocity'])),np.finfo(float).tiny)
                meta=dict(zone=zone,stride=stride,forcing=forcing)
                check('baryonic mass flux constraint',fractions@sol['velocity'][:-1]/scale,2e-13,**meta)
                check('electric current constraint',(ni*charge@sol['velocity'][:-1]-ne*sol['velocity'][-1])/(2*ne*scale),2e-13,**meta)
                check('entropy balance',sol['entropy_from_collisions']/sol['entropy_from_forces']-1,**meta)
                flux=np.r_[sol['mass_flux'],sol['reduced_heat_flux']/op.energy_scale]
                norm=np.sqrt(np.diag(op.mobility))
                check('full solve currents at stellar forces',np.max(abs(full_flux-flux)/norm)/max(np.max(abs(flux)/norm),1e-300),**meta)
                q=h@sol['mass_flux']-kappa*t*gradient['dlnT_dr']
                check('matched diffusion and heat flux',q/sol['reduced_heat_flux']-1,**meta)
                solutions[f'{stride}:{forcing}']=dict(velocity_cm_s=sol['velocity'][:3].tolist(),
                       ten_mode_no_pair_velocity_cm_s=base['velocity'][:3].tolist(),
                       thermal_velocity_contribution_cm_s=(sol['velocity']-noheat['velocity'])[:3].tolist(),
                       independent_mass_flux_g_cm2_s=sol['mass_flux'].tolist(),
                       reduced_heat_flux_erg_cm2_s=sol['reduced_heat_flux'],
                       full_material_energy_flux_erg_cm2_s=float(sol['reduced_heat_flux']+np.asarray(fr['exchange_enthalpy_erg_g'])@sol['mass_flux']),
                       pressure_scale_crossing_years=[float(fr['pressure_scale_height_cm']/abs(v)/(365.25*86400)) if v else None for v in sol['velocity'][:3]])
        records.append(dict(zone=zone,temperature_K=t,density_baryonic_g_cm3=rho,X=mr['X'],
                   eta=calc['eta'],b_thermal=calc['b_thermal'],screening=sample['screening'],
                   electron_electron_prefactor_ratio=reduced.electron_electron_prefactor_ratio,
                   retained_electron_response=reduced.retained_response.tolist(),
                   conductivity_zero_element_flux_cgs=kappa,conductivity_over_no_pair_ten_modes=kappa/koff,
                   conductivity_over_no_pair_exact=kappa/kexact,no_pair_ten_modes_over_exact=koff/kexact,
                   transport_enthalpy_erg_g=h.tolist(),mode_convergence=history,solutions=solutions,
                   scaled_condition_number=op.scaled_condition_number))
    passed=all(c['passed'] for c in checks)
    report=dict(outcome='completed_conditional_mixture_comparison' if passed else 'failed',
                accepted_for_stellar_evolution=False,created_utc=datetime.now(timezone.utc).isoformat(),
                checks=checks,records=records,input_sha256=inputs,new_EOS_queries=0,new_collision_integrals=0,
                limitations=['Prescribed common static Born screening, nonrelativistic electrons, classical ions and leading Brownian ion heat term.',
                             'Fully stripped, stationary metals exchange heat; ionization, higher-order recoil and mixture correlations remain unresolved.',
                             'EOS forces and kinetic heat convention are conditionally matched as in the retained stellar-force comparison.',
                             'Three sampled zones and one screening prescription do not establish full-track coverage.',
                             'No selected conductivity replacement, conservative abundance or energy update, or structural evolution.'])
    with OUTPUT.open('x') as stream:json.dump(report,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(dict(outcome=report['outcome'],checks=len(checks),failed=[c for c in checks if not c['passed']],
               summaries=[{k:r[k] for k in ['zone','electron_electron_prefactor_ratio','conductivity_over_no_pair_ten_modes',
                                           'conductivity_over_no_pair_exact','no_pair_ten_modes_over_exact']} for r in records],
               report_sha256=digest(OUTPUT)),indent=2))
    if not passed:raise SystemExit(1)


if __name__ == '__main__':
    main()
