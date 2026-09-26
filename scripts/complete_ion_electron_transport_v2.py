#!/usr/bin/env python3
"""Complete only the previously unsupported ion/electron transport cases.

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
from ion_collision_table import ExtendedIonTable
from diffusion_burgers import KB, MU, ME
from ion_collision_integrals import E2
from diffusion_thermal_transport import HeatFluxConvention
from ion_electron_transport import coupled_elastic_operator, heat_decomposition
from electron_ion_heat import physical_prefactor


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    previous_path=Path('docs/results/ion_electron_transport_v1.json')
    previous=json.loads(previous_path.read_text())
    assert previous['outcome']=='failed' and len(previous['unsupported_physical_cases'])==6
    assert all(c['passed'] for c in previous['checks'])
    for path,digest in previous['inputs_sha256'].items():assert sha(path)==digest,path
    checks=list(previous['checks'])
    limits=list(previous['classical_mass_ratio_limits'])
    physical=list(previous['physical_cases'])
    completed={(c['zone'],c['screening']) for c in physical}
    extension_path=Path('docs/results/ion_collision_table_extension_v1.json')
    extension=json.loads(extension_path.read_text())
    assert extension['status']=='pass' and extension['old_queries_byte_identical']
    for path,digest in extension['input_sha256'].items():assert sha(path)==digest,path
    assert sha(extension['table'])==extension['table_sha256']
    extended=ExtendedIonTable(extension['table'])

    def check(name, actual, expected, tolerance=2e-10, **meta):
        a, e = np.asarray(actual), np.asarray(expected)
        scale = float(np.max(abs(e)))
        if scale == 0: scale = 1.
        error = float(np.max(abs(a-e))/scale)
        checks.append(dict(name=name,relative_error=error,tolerance=tolerance,
                           passed=bool(np.isfinite(error) and error <= tolerance),**meta))

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
    for p in [previous_path,extension_path,Path(extension['table']),Path('scripts/ion_collision_table.py')]:inputs[str(p)]=sha(p)
    inputs.update(previous['inputs_sha256'])
    inputs.update(extension['input_sha256'])
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
    spline = extended.log_moments
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
            if (index,label) in completed:continue
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
    passed=all(c['passed'] for c in checks) and not unsupported and len(physical)==12
    for p,h in inputs.items():assert sha(p)==h,p
    report=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),
        outcome='passed' if passed else 'failed',accepted_for_stellar_evolution=False,
        checks=checks,classical_mass_ratio_limits=limits,physical_cases=physical,
        previous_checks_reused=len(previous['checks']),previous_physical_cases_reused=len(completed),
        new_physical_cases=len(physical)-len(completed),
        unsupported_physical_cases=unsupported,inputs_sha256=inputs,
        limitations=['No electron-electron collisions or higher-order ion recoil.',
                     'Fixed fully stripped metals and prescribed common Born screening.',
                     'Classical ion moments from the checked screened-potential table; no additional mixture-correlation correction.',
                     'Forces are prescribed comparison directions, not saved stellar chemical gradients.',
                     'No stellar velocity, selected conductivity replacement or abundance evolution.'])
    with Path('docs/results/ion_electron_transport_v2.json').open('x') as stream:
        json.dump(report,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(dict(outcome=report['outcome'],checks=len(checks),physical_cases=len(physical),
        failures=[c for c in checks if not c['passed']],unsupported=unsupported)))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
