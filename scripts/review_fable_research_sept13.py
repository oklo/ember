#!/usr/bin/env python3
"""Independent scalar review of retained Fable encounter and thermal results.

No SPH simulation, particle potential sum, or stellar evolution is rerun.
"""
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
FABLE=ROOT/'docs/research/fable'
RESULTS=FABLE/'results'
OUTPUT=ROOT/'docs/results/fable_research_review_sept13_v1.json'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    if OUTPUT.exists():raise FileExistsError(OUTPUT)
    inputs={str(Path(__file__)):sha(__file__)}
    def read(p):
        inputs[str(p)]=sha(p);return json.loads(p.read_text())
    batches=[]
    for p in sorted(RESULTS.glob('*/READY.json')):
        d=read(p);bad=[]
        for name,h in d.get('files_sha256',{}).items():
            q=p.parent/name
            if not q.exists() or sha(q)!=h:bad.append(name)
            else:inputs[str(q)]=h
        historical=[]
        for name,h in d.get('source_sha256',{}).items():
            q=FABLE/name
            if not q.exists() or sha(q)!=h:
                historical.append(dict(source=name,expected_sha256=h,current_sha256=sha(q) if q.exists() else None))
        batches.append(dict(batch=p.parent.name,verified_files=len(d.get('files_sha256',{})),
            mismatches=bad,historical_source_differences=historical))
    assert all(not b['mismatches'] for b in batches)
    # Recalculate energy changes from the archived summaries, independently
    # of FINDINGS tables or prose. The particle membership remains the
    # worker's kinetic-plus-gravitational criterion, not a new gas-binding test.
    cases={}
    ids=['001','006','007','008','014','016','017','018','019','020','021','022','023','025','026']
    for case in ids:
        matches=sorted(RESULTS.glob(f'fable-batch-[123]/F-SPH-{case}_*/analysis_v1.json'))
        assert len(matches)==1,(case,matches)
        d=read(matches[0]);a,b=d['snapshots'][0],d['snapshots'][-1]
        loss=a['orbit']['E_orb_pointmass']-b['orbit']['E_orb_pointmass']
        gains=[b[f'remnant{i}']['E_body']-a[f'remnant{i}']['E_body'] for i in [1,2]]
        cases[case]=dict(source=str(matches[0]),initial_orbital_energy=a['orbit']['E_orb_pointmass'],
            final_orbital_energy=b['orbit']['E_orb_pointmass'],orbital_loss=loss,
            radial_velocity_out=b['orbit']['radial_velocity'],stellar_gains=gains,
            closure_relative_to_orbital_loss=(sum(gains)-loss)/loss,
            star1_heat=b['remnant1']['U_thermal'],star1_mass_initial=a['remnant1']['M'],
            star1_mass_final=b['remnant1']['M'],units_energy_J=1e41)
    controls=[];ref=cases['006']
    for case in ['016','020','021','022']:
        row=cases[case]
        controls.append(dict(case=case,
            outgoing_orbit_relative_difference=row['final_orbital_energy']/ref['final_orbital_energy']-1,
            mean_deposit_relative_difference=sum(row['stellar_gains'])/sum(ref['stellar_gains'])-1))
    mode_file=RESULTS/'tidal_linear_n15_equal_mass_v1.json';linear=read(mode_file)
    saved={r['q']:r['dE_total_GM2_over_R'] for r in linear['rows'] if r['w']==.1}
    slope=math.log(saved[2.]/saved[1.5])/math.log(2./1.5)
    sensitivity=dict(q_interval=[1.5,2.],w=.1,mean_log_energy_slope_vs_q=slope,
        inferred_physical_energy_radius_exponent_at_fixed_pericenter=-1-slope,
        meaning='Finite-interval scale sensitivity from the retained fluid model, not a local derivative or an elastic-star error bound.')
    # Re-evaluate the actual sealed thermal code and its known latent-energy
    # interval. It is discontinuous in U(T); a scalar root cannot represent
    # the missing phase fraction.
    common=FABLE/'sph/fable_common.py';inputs[str(common)]=sha(common)
    sys.path.insert(0,str(common.parent))
    thermal_path=RESULTS/'fable-batch-3-thermal-correction-v1/thermal_response_v2.py'
    spec=importlib.util.spec_from_file_location('sealed_fable_thermal_review',thermal_path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    remnant=module.Remnant();us=remnant.U(remnant.T_m);latent=remnant.N_ion*remnant.L_ion
    gap=[]
    for fraction in [.1,.5,.9]:
        u=us+fraction*latent;t=remnant.T_of_U(u)
        gap.append(dict(required_liquid_fraction=fraction,target_energy_J=u,returned_temperature_K=t,
            relative_energy_residual=(remnant.U(t)-u)/u,
            phase_fraction_representation_energy_J=us+fraction*latent))
    assert max(abs(r['relative_energy_residual']) for r in gap)>.2
    # Baiko & Chugunov (2022), equation 40, with the ion Bohr radius of
    # their section 2. Only evaluate inside the published fit's domain.
    hbar,mu,e,eps0,kb=1.054571817e-34,1.66053906660e-27,1.602176634e-19,8.8541878128e-12,1.380649e-23
    A,Z=4.,2.;ab=hbar*hbar*4*math.pi*eps0/(A*mu*(Z*e)**2)
    quantum=[]
    for rho in [6782.7,40630.]:
        ni=rho*1000/(A*mu);a=(3/(4*math.pi*ni))**(1/3);rs=a/ab
        assert 500<=rs<=120000
        gamma=175.7-1300/rs+4.1e6/rs**2
        classical=(Z*e)**2/(4*math.pi*eps0*a*kb*175)
        quantum.append(dict(rho_cgs=rho,ion_rs=rs,Gamma_m=gamma,
            T_m_classical_175_K=classical,T_m_fit_K=classical*175/gamma,
            relative_melting_temperature_change=175/gamma-1))
    rate_file=RESULTS/'fable-batch-3-thermal-correction-v1/encounter_rates_v2.json';rates=read(rate_file)
    row=next(r for r in rates['rows'] if r['n_pc3']==.1 and r['sigma_kms']==10 and r['q_max']==3.2)
    G,M,pc,year=6.67430e-11,.2*1.98847e30,3.0857e16,3.15576e7
    r0=row['r0_m'];s=math.sqrt(2)*1e4;density=.1/pc**3;focus=2*G*M/r0
    mean_v=math.sqrt(8/math.pi)*s;inv_v=math.sqrt(2/math.pi)/s
    rate=density*math.pi*r0*r0*(mean_v+focus*inv_v)*year
    assert abs(rate/row['rate_per_yr']-1)<1e-12
    event_mean=(3*s*s+focus)/(mean_v+focus*inv_v)
    rate_review=dict(all_passages_waiting_time_yr=1/rate,events_by_1e14yr=rate*1e14,
        field_mean_speed_m_s=mean_v,encounter_weighted_mean_speed_m_s=event_mean,
        measured_local_WD_density_pc3=.00447,
        same_model_waiting_time_at_measured_local_WD_density_yr=1/rate*.1/.00447,
        limitation='The published density counts local WD systems; a future mass-resolved perturber population is a separate scenario. This rate includes collisions and captures.')
    corrections=dict(q3_orbital_noise_over_q2_loss=abs(cases['023']['orbital_loss']/cases['017']['orbital_loss']),
        q15_star1_gain_J=cases['008']['stellar_gains'][0]*1e41,
        stripping_origin1_particles_mass_msun=8915/49986*.1)
    source_refs={
        'mode_damping':'https://academic.oup.com/mnras/article/357/3/834/1077414',
        'quantum_melting_fit':'https://academic.oup.com/mnras/article/510/2/2628/6464199',
        'local_WD_density':'https://doi.org/10.1093/mnras/stab2672',
        'additional_collision_literature':'https://arxiv.org/abs/2308.14449'}
    output=dict(scope=__doc__,outcome='survey_useful_with_required_interpretation_and_thermal_corrections',
        accepted_for_stellar_evolution=False,input_sha256=inputs,batches=batches,cases=cases,
        clean_survivor_controls=controls,fluid_radius_sensitivity=sensitivity,
        thermal=dict(U100_J=remnant.U(100),Tp_K=remnant.T_p,Tm_K=remnant.T_m,
            latent_heat_J=latent,energy_gap_cases=gap),quantum_melting_comparison=quantum,
        encounter_rates=rate_review,prose_arithmetic_corrections=corrections,primary_sources=source_refs,
        limitations=['No particle simulation is rerun; numerical controls apply to the cases tested and do not validate a realistic cold remnant.',
            'Quantum fit checks melting temperatures only, not the worker liquid heat capacity, latent heat, or post-pulse temperature.',
            'The literature identifies gravitational radiation as a competing mode-energy loss; its dominance for this particular helium remnant is not calculated here.',
            'Historical source hash differences identify versions to recover; they do not establish that no copy exists elsewhere.'])
    for p,h in inputs.items():assert sha(p)==h,p
    with OUTPUT.open('x') as f:json.dump(output,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(report=str(OUTPUT),sha256=sha(OUTPUT),verified_files=sum(b['verified_files'] for b in batches),
        maximum_clean_control_change=max(abs(r['mean_deposit_relative_difference']) for r in controls),
        latent_gap_maximum_relative_residual=max(abs(r['relative_energy_residual']) for r in gap),
        quantum_melting=quantum,rates=rate_review,arithmetic=corrections)),flush=True)


if __name__=='__main__':main()
