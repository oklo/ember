#!/usr/bin/env python3
"""Check Brownian ion heat relaxation, the classical limit and saved-layer responses.

Physical-layer outputs remain comparisons for a specified collision model.
No stellar driving gradient, settling velocity or evolved abundance is inferred.
"""
import hashlib
import json
import math
from pathlib import Path
import re
from datetime import datetime, timezone

import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from scipy.interpolate import CubicSpline
from diffusion_burgers import KB, MU, ME
from ion_collision_integrals import E2
from diffusion_reciprocal_transport import maxwellian_operator
from diffusion_thermal_transport import HeatFluxConvention
from ion_electron_transport import coupled_elastic_operator, heat_decomposition
from electron_ion_heat import physical_prefactor


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    checks, limits, physical = [], [], []

    def check(name, actual, expected, tolerance=2e-10, **meta):
        a, e = np.asarray(actual), np.asarray(expected)
        scale = float(np.max(abs(e)))
        if scale == 0: scale = 1.
        error = float(np.max(abs(a-e))/scale)
        checks.append(dict(name=name,relative_error=error,tolerance=tolerance,
                           passed=bool(np.isfinite(error) and error <= tolerance),**meta))

    # Independent velocity-space Gaussian integrals, not a Burgers matrix.
    c,w = hermegauss(8)
    grid = np.stack(np.meshgrid(c,c,c,indexing='ij'),axis=-1).reshape(-1,3)
    weight = np.einsum('i,j,k->ijk',w,w,w).ravel()/(2*math.pi)**1.5
    square = np.sum(grid*grid,axis=1)
    drift = grid[:,2]
    heat = drift*(square-5)/5
    grad_drift = np.zeros_like(grid); grad_drift[:,2] = 1.
    grad_heat = 2*drift[:,None]*grid/5
    grad_heat[:,2] += (square-5)/5
    check('Gaussian drift current',weight@(drift*drift),1.)
    check('Gaussian heat has no number current',weight@(drift*heat),0.)
    check('Gaussian residual heat normalization',weight@(drift*(square/2-2.5)*heat),1.)
    check('Brownian heat relaxation relative to drag',weight@np.sum(grad_heat**2,axis=1),1.2)
    check('Brownian drift/heat cross cost',weight@np.sum(grad_drift*grad_heat,axis=1),0.)

    # Infinite-ion-mass limit of independently checked finite-mass Burgers
    # equations, with the actual Coulomb electron heat moment ratios.
    I = np.array([[1.,-.6],[-.6,.52]])
    L = np.array([[6.,9.],[9.,37.5]])/math.gamma(2.5)**2
    for trace in [.15,1e-12]:
        x = np.array([trace,.05,.95-trace]);z = np.array([1.,2.,2.])
        previous = None
        for mass_scale in [1.,16.,256.]:
            mass = np.array([1.,3.,4.])*mass_scale
            density,temperature = 1300.,1.1e7
            ni = density/MU*x/mass
            ke = physical_prefactor(1.,1.)*ni*z*z
            ki = np.sqrt(np.outer(ke,ke))*4
            cz,zp,zpp = [np.full((3,3),v) for v in [.35,1.2,1.4]]
            common = dict(density=density,temperature=temperature,mass_fractions=x,
                mass_numbers=mass,charges=z,independent_ions=[0,1],reference_ion=2,
                ion_resistance=ki,ion_z=cz,ion_zprime=zp,ion_zdoubleprime=zpp,
                electron_moments=I,electron_full_response=L,energy_scale=KB/MU*temperature)
            op = coupled_elastic_operator(**common,electron_relaxation='one_heat_variable')
            k = np.zeros((4,4));k[:3,:3]=ki;k[:3,3]=ke;k[3,:3]=ke
            zz,pp,p2 = [np.full((4,4),v) for v in [.6,1.3,2.]]
            zz[:3,:3]=cz;pp[:3,:3]=zp;p2[:3,:3]=zpp
            old = maxwellian_operator(density=density,temperature=temperature,
                mass_fractions=x,mass_numbers=mass,charges=z,independent_ions=[0,1],
                reference_ion=2,resistance=k,collision_z=zz,collision_zprime=pp,
                collision_zdoubleprime=p2,energy_scale=common['energy_scale'])
            scale = np.sqrt(np.diag(op.collision_matrix))
            error = float(np.max(abs(old.collision_matrix-op.collision_matrix)/scale[:,None]/scale[None,:]))
            bound = 10*ME/MU/min(mass)
            check('finite-mass Burgers approaches Brownian elastic form',error,0.,bound,
                  hydrogen=trace,mass_scale=mass_scale)
            if previous is not None:
                check('mass-ratio error decreases',max(0.,error/previous-.07),0.,1e-12,
                      hydrogen=trace,mass_scale=mass_scale)
            previous=error
            limits.append(dict(hydrogen=trace,mass_scale=mass_scale,
                               normalized_matrix_difference=error,maximum_mass_ratio=ME/MU/min(mass)))

    paths = [Path(__file__),Path('scripts/ion_electron_transport.py'),
        Path('scripts/electron_elastic_mixture.py'),Path('scripts/electron_ion_heat.py'),
        Path('scripts/diffusion_composition_forces.py'),Path('scripts/diffusion_burgers.py'),
        Path('scripts/diffusion_reciprocal_transport.py'),Path('scripts/diffusion_thermal_transport.py'),
        Path('docs/results/diffusion_reciprocal_classical_v1.json'),
        Path('docs/results/electron_elastic_mixture_v1.json'),
        Path('docs/results/stellar_electron_heat_3890gyr_v1.json'),
        Path('docs/results/diffusion_material_regime_3890gyr_v1.json'),
        Path('docs/results/yukawa_collision_table_v1.json'),Path('include/ember/gs98_mixture.hpp')]
    inputs = {str(p):sha(p) for p in paths}
    classical,elastic,heat_profile,regime,table_report = [json.loads(p.read_text()) for p in paths[8:13]]
    for report in [classical,elastic,heat_profile,regime,table_report]:
        for key in ['inputs_sha256','input_sha256']:
            for p,h in report.get(key,{}).items():
                assert sha(p) == h,p
                inputs[p]=h
    assert classical['outcome']=='passed' and elastic['outcome']=='passed'
    table_path = Path(table_report['table'])
    assert sha(table_path)==table_report['table_sha256'] and table_report['status']=='pass'
    inputs[str(table_path)]=sha(table_path)
    table = json.loads(table_path.read_text())
    spline = CubicSpline(table['log10_strength'],np.log(table['dimensionless_integrals']),axis=0,extrapolate=False)
    metals = [tuple(map(float,s.split(','))) for s in re.findall(r'^\s*\{([^{}]+)\}, //',paths[-1].read_text(),re.M)]
    mass = np.array([1.,3.,4.]+[m[1] for m in metals])
    z = np.array([1.,2.,2.]+[m[0] for m in metals])
    unsupported=[]
    for index in [0,100,234,300,396,421]:
        row,hr = regime['records'][index],heat_profile['records'][index]
        assert row['zone']==hr['zone']
        T,rho = row['temperature_K'],row['density_baryonic_g_cm3']
        metal_total=1-row['X']-row['Y3']-row['Y4']
        x=np.array([row['X'],row['Y3'],row['Y4']]+[metal_total*m[3] for m in metals])
        ni=rho/MU*x/mass;ne=float(ni@z)
        for sample in hr['cases']:
            label=sample['screening'];lam=sample['screening_length_cm']
            strength=np.outer(z,z)*E2/(KB*T*lam)
            moments=np.exp(spline(np.log10(strength)))
            if not np.all(np.isfinite(moments)):
                unsupported.append(dict(zone=index,screening=label,strength_min=float(strength.min()),strength_max=float(strength.max())))
                continue
            reduced=MU*(mass[:,None]*mass[None,:]/(mass[:,None]+mass[None,:]))
            omega11=np.sqrt(2*math.pi/reduced)*(np.outer(z,z)*E2)**2/(KB*T)**1.5*moments[:,:,0]
            ki=16/3*np.outer(ni,ni)*reduced*omega11
            cz=1-.4*moments[:,:,1]/moments[:,:,0]
            zp=2.5-2*moments[:,:,1]/moments[:,:,0]+.4*moments[:,:,2]/moments[:,:,0]
            zpp=moments[:,:,3]/moments[:,:,0]
            common=dict(density=rho,temperature=T,mass_fractions=x,mass_numbers=mass,
                charges=z,independent_ions=[0,1],reference_ion=2,ion_resistance=ki,
                ion_z=cz,ion_zprime=zp,ion_zdoubleprime=zpp,
                electron_moments=sample['dimensionless_collision_moments'],
                electron_full_response=sample['dimensionless_full_energy_response'],energy_scale=KB/MU*T)
            op=coupled_elastic_operator(**common)
            trial=coupled_elastic_operator(**common,electron_relaxation='one_heat_variable')
            h,kappa=heat_decomposition(op.mobility,temperature=T,energy_scale=op.energy_scale)
            ht,kt=heat_decomposition(trial.mobility,temperature=T,energy_scale=op.energy_scale)
            meta=dict(zone=index,screening=label)
            scale=np.sqrt(np.diag(op.mobility));a=op.mobility/scale[:,None]/scale[None,:]
            check('coupled reciprocity',a,a.T,**meta)
            diag=np.sqrt(np.diag(trial.collision_matrix))
            difference=(trial.collision_matrix-op.collision_matrix)/diag[:,None]/diag[None,:]
            check('energy relaxation reduces collision cost',max(0.,-float(np.linalg.eigvalsh(difference).min())),0.,2e-11,**meta)
            check('full energy increases zero-element-flux conductivity',max(0.,kt/kappa-1),0.,2e-11,**meta)
            enthalpy=np.array([.3,-.2])*op.energy_scale
            convention=HeatFluxConvention(enthalpy,op.energy_scale)
            for f in [np.array([1.,0.,0.]),np.array([0.,1.,0.]),np.array([0.,0.,1.]),np.array([.4,-.7,.2])]:
                force=f/scale
                sol=op.solve(force)
                speed=float(np.max(abs(sol['velocity'])))
                if speed:
                    check('zero baryonic flux',abs(x@sol['velocity'][:-1])/speed,0.,1e-13,**meta)
                    current=ni*z@sol['velocity'][:-1]-ne*sol['velocity'][-1]
                    check('zero electron current',abs(current)/(2*ne*speed),0.,1e-13,**meta)
                check('stationary metal velocities',sol['velocity'][3:-1],0.,**meta)
                check('coupled entropy balance',sol['entropy_from_forces'],sol['entropy_from_collisions'],**meta)
                check('positive coupled dissipation',max(0.,-sol['entropy_from_forces']),0.,**meta)
                gradT=-force[-1]*T*T/op.energy_scale
                direct=h@sol['mass_flux']-kappa*gradT
                check('heat decomposition versus collision solution',direct,sol['reduced_heat_flux'],**meta)
                reduced=np.r_[sol['mass_flux'],direct/op.energy_scale]
                full=convention.full_fluxes(reduced)
                check('full material energy convention',convention.full_matrix(op.mobility)@convention.full_forces(force),full,**meta)
            rescaled=coupled_elastic_operator(**{**common,'energy_scale':op.energy_scale*13})
            hs,ks=heat_decomposition(rescaled.mobility,temperature=T,energy_scale=rescaled.energy_scale)
            check('conductivity independent of numerical energy scale',ks,kappa,**meta)
            check('transport enthalpy independent of numerical energy scale',hs,h,**meta)
            physical.append(dict(**meta,temperature_K=T,density_g_cm3=rho,
                strength_min=float(strength.min()),strength_max=float(strength.max()),
                conductivity_zero_element_flux_cgs=kappa,one_heat_over_full_conductivity=kt/kappa,
                transport_enthalpy_erg_g=h.tolist(),scaled_condition=op.scaled_condition_number,
                scaled_backward_error=op.scaled_backward_error,
                ion_heat_electron_relaxation_over_ion_diagonal=(op.brownian_heat_diagonal/np.diag(op.ion_collision_matrix)[2:-1]).tolist()))
    guards=0
    for matrix,t,e in [(np.eye(3),0.,1.),(np.eye(3),1.,0.),
                       (np.array([[1.,2.],[2.,1.]]),1.,1.),
                       (np.array([[1.,.1],[0.,1.]]),1.,1.)]:
        try:heat_decomposition(matrix,temperature=t,energy_scale=e)
        except ValueError:guards+=1
    check('heat decomposition invalid input guards',guards,4,0.)
    passed=all(c['passed'] for c in checks) and not unsupported and len(physical)==12
    for p,h in inputs.items():assert sha(p)==h,p
    report=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),
        outcome='passed' if passed else 'failed',accepted_for_stellar_evolution=False,
        checks=checks,classical_mass_ratio_limits=limits,physical_cases=physical,
        unsupported_physical_cases=unsupported,inputs_sha256=inputs,
        limitations=['No electron-electron collisions or higher-order ion recoil.',
                     'Fixed fully stripped metals and prescribed common Born screening.',
                     'Classical ion moments from the checked screened-potential table; no additional mixture-correlation correction.',
                     'Forces are prescribed comparison directions, not saved stellar chemical gradients.',
                     'No stellar velocity, selected conductivity replacement or abundance evolution.'])
    with Path('docs/results/ion_electron_transport_v1.json').open('x') as stream:
        json.dump(report,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(dict(outcome=report['outcome'],checks=len(checks),physical_cases=len(physical),
        failures=[c for c in checks if not c['passed']],unsupported=unsupported)))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
