#!/usr/bin/env python3
"""Bound missing OP opacity components in a finite-frequency plasma comparison.

Reuse corrected native OP spectra and a retained full TOPS control. Compare
continuous frequency integration with native OP means, propagate component
intervals through the signed thermodynamic interpolation, and keep frequency
tails and source-interpolation errors separate from the component bounds.
"""
from concurrent.futures import ProcessPoolExecutor
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from op_native_spectra_v2 import read_mesh
from op_plasma_transport import finite_means, log_lagrange_weights, interpolate_log_interval
from audit_tops_electron_dispersion import electron_moments, HBAR, ME, KB, KEV, RW
from audit_tops_spectral_means import source

ROOT=Path('/Users/greglaughlin/Projects/ember')
WORK=Path('/tmp/ember-op-plasma-transport-run-v1')
OUTPUT=ROOT/'docs/results/op_plasma_transport_v1.json'
PREVIOUS=ROOT/'docs/results/op_mixture_corrections_v1.json'
BASE=ROOT/'docs/results/op_native_mixture_comparison_v1.json'
MESH=Path('/tmp/ember-op-selected-spectra-v1/OP4STARS_1.3/mono/m01.mesh')
A0=5.29177210903e-9
MU=1.66053906660e-24
ELECTRON_CHARGE=4.80320471257e-10


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def analytic_controls():
    checks=[]
    def add(name,value,limit):
        checks.append(dict(name=name,value=float(value),limit=limit,passed=bool(abs(value)<=limit)))
    u=np.unique(np.r_[np.geomspace(.05,20.,97),1.])
    result=finite_means(u,np.full(len(u),6.5),np.full(len(u),.5),1.,0.,32,
                        actual_absorption=np.full(len(u),3.))
    def wr(x):return x**4*math.exp(-x)/(-math.expm1(-x))**2
    base=RW/quad(lambda x:wr(x)/7.,u[0],u[-1],epsabs=1e-13,epsrel=1e-12)[0]
    add('grey_finite_mean',result['uncut_finite_mean']/base-1,2e-11)
    for name,absorption in [('scattering_extreme_mean',0.),('absorption_extreme_mean',6.5),('actual_component_mean',3.)]:
        scattering=7.-absorption
        def f(x):
            n=math.sqrt(1-1/x**2)
            return wr(x)*n**3/(absorption+n*scattering)
        reference=RW/quad(f,1.,20.,epsabs=1e-13,epsrel=1e-12,limit=200)[0]
        add('grey_'+name,result[name]/reference-1,2e-8)
    for target in (1.1,1.8,3.7):
        weights=log_lagrange_weights([1.,1.5,2.5,4.],target)
        low=np.array([1.01,1.02,1.03,1.04]);high=low+.01
        interval=interpolate_log_interval(weights,low,high)
        corners=[math.exp(weights @ np.log(np.where(bits,high,low))) for bits in itertools.product((False,True),repeat=4)]
        add(f'interval_min_{target:g}',interval['lower']/min(corners)-1,2e-14)
        add(f'interval_max_{target:g}',interval['upper']/max(corners)-1,2e-14)
    return checks


def native_calculation(task):
    row,total,free_scattering=task
    _,u=read_mesh(MESH)
    remainder=total-free_scattering
    if np.any(remainder<=0):raise ValueError('known scattering exceeds total opacity')
    conversion=A0*A0/MU*(-np.expm1(-u))
    rest=remainder*conversion;scattering=free_scattering*conversion
    t=row['temperature_K'];ne=row['electron_density_cm3']
    electron=electron_moments(t,ne)
    cutoff=HBAR*math.sqrt(4*math.pi*ELECTRON_CHARGE**2*ne/ME)/(KB*t)*math.sqrt(electron['plasma_frequency_squared_ratio'])
    trials=[]
    for order in (8,16,32):
        calculated=finite_means(u,rest,scattering,cutoff,electron['vstar_squared'],order)
        trials.append(dict(order=order,**calculated))
        if len(trials)>1:
            keys=['uncut_finite_mean','scattering_extreme_mean','absorption_extreme_mean']
            error=max(abs(trials[-1][k]/trials[-2][k]-1) for k in keys)
            if error<1e-8:break
    if error>=1e-8:raise ValueError('finite spectral quadrature did not converge')
    best=trials[-1]
    return dict(label=row['label'],temperature_index=row['temperature_index'],electron_index=row['electron_index'],
        temperature_K=t,electron_density_cm3=ne,electron_response=electron,quadrature_trials=trials,
        quadrature_relative_change=error,final=best,
        continuous_vs_native_mean=best['uncut_finite_mean']/row['means'][-1]['rosseland_cm2_g']-1,
        unknown_component_interval_below_half_percent=best['component_interval_fraction']<.005)


def case_interpolations(case,lookup):
    results=[]
    for nt in (2,4):
        planes=[case['planes'][i] for i in ([1,2] if nt==2 else [0,1,2,3])]
        tw=log_lagrange_weights([p['temperature_K'] for p in planes],case['temperature_K'])
        for nr in (2,4):
            weights=[];samples=[]
            for plane,wt in zip(planes,tw):
                density=next(r for r in plane['density_interpolations'] if r['density_points']==nr)
                rw=log_lagrange_weights(density['native_densities'],case['density_g_cm3'])
                for jn,wr in zip(density['electron_indices'],rw):
                    samples.append(lookup[case['label'],plane['temperature_index'],jn])
                    weights.append(wt*wr)
            weights=np.array(weights)
            base=math.exp(weights @ np.log([r['final']['uncut_finite_mean'] for r in samples]))
            low=np.array([r['final']['ratio_lower'] for r in samples])
            high=np.array([r['final']['ratio_upper'] for r in samples])
            bound=interpolate_log_interval(weights,low,high)
            endmembers=dict(scattering=math.exp(weights @ np.log(low)),absorption=math.exp(weights @ np.log(high)))
            ne=math.exp(weights @ np.log([r['electron_density_cm3'] for r in samples]))
            results.append(dict(temperature_points=nt,density_points=nr,uncut_finite_mean=base,
                interpolated_electron_density_cm3=ne,ratio_bounds_for_this_interpolation=bound,
                bounds_interval_fraction=bound['upper']/bound['lower']-1,
                endmember_ratios=endmembers,endmember_means={k:base*v for k,v in endmembers.items()},
                component_mean_bounds=dict(lower=base*bound['lower'],upper=base*bound['upper']),
                interpolation_terms=[dict(temperature_index=r['temperature_index'],electron_index=r['electron_index'],weight=float(w)) for r,w in zip(samples,weights)]))
    return dict(label=case['label'],temperature_K=case['temperature_K'],density_g_cm3=case['density_g_cm3'],interpolations=results)


def tops_comparison(inputs):
    report=ROOT/'docs/results/tops_cool_full_spectra_v4.json'
    inputs[str(report)]=digest(report)
    data=json.loads(report.read_text())
    row=next(r for r in data['records'] if r['X']==.1 and r['temperature_keV']==.008 and r['density_atomic_g_cm3']==.053367)
    root=Path(row['full_spectrum_source'])
    for p in [root/'source.txt',root/'receipt.json',root/'request.json']:inputs[str(p)]=digest(p)
    spectral,_=source(root);s=spectral[(.008,.053367)]
    electron=electron_moments(.008*KEV/KB,row['free_electron_density_cm3'])
    _,op_u=read_mesh(MESH)
    trials=[]
    for order in (8,16):
        means=finite_means(s[:,0]/.008,s[:,2],s[:,3],row['cutoff_u'],electron['vstar_squared'],order,lower=op_u[0],upper=op_u[-1],actual_absorption=s[:,2])
        trials.append(dict(order=order,**means))
    best=trials[-1]
    error=max(abs(trials[-1][k]/trials[-2][k]-1) for k in ['uncut_finite_mean','actual_component_mean'])
    if error>=1e-8:raise ValueError('TOPS finite-interval quadrature did not converge')
    return dict(temperature_keV=.008,density_g_cm3=.053367,electron_density_cm3=row['free_electron_density_cm3'],
        electron_response=electron,quadrature_trials=trials,quadrature_relative_change=error,
        final=best,full_source_mean=row['reference_rosseland_atomic_cm2_g'],full_source_uncut_mean=row['reference_uncut_rosseland_atomic_cm2_g'],
        full_source_ratio=row['reference_refractive_ratio'],
        finite_interval_mean_change=best['actual_component_mean']/row['reference_rosseland_atomic_cm2_g']-1,
        finite_interval_ratio_change=best['ratio_upper']/row['reference_refractive_ratio']-1)


def main():
    if OUTPUT.exists():raise FileExistsError(OUTPUT)
    inputs={str(p):digest(p) for p in [Path(__file__).resolve(),ROOT/'scripts/op_plasma_transport.py',ROOT/'scripts/op_native_spectra_v2.py',
        ROOT/'scripts/audit_tops_electron_dispersion.py',ROOT/'scripts/audit_tops_spectral_means.py',PREVIOUS,BASE,MESH]}
    corrections=json.loads(PREVIOUS.read_text());base=json.loads(BASE.read_text())
    if corrections['outcome']!='completed_conditional_comparison' or not all(r['passed'] for r in corrections['controls']):raise ValueError('mixture correction checks failed')
    artifact=corrections['spectral_artifact'];p=Path(artifact['path'])
    if digest(p)!=artifact['sha256']:raise ValueError('corrected spectra changed')
    inputs[str(p)]=digest(p)
    checks=analytic_controls()
    tasks=[]
    with np.load(p) as spectra:
        for row in corrections['native_states']:
            if not row['used_in_mean_interpolation']:raise ValueError('nonpositive corrected source')
            prefix=f"{row['label']}_t{row['temperature_index']}_n{row['electron_index']}"
            tasks.append((row,spectra[prefix+'_total'].copy(),spectra[prefix+'_free_electron_scattering'].copy()))
    records=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for record in pool.map(native_calculation,tasks):records.append(record)
    lookup={(r['label'],r['temperature_index'],r['electron_index']):r for r in records}
    cases=[case_interpolations(case,lookup) for case in base['cases']]
    tops=tops_comparison(inputs)
    for row in cases[0]['interpolations']:
        row['relative_to_TOPS_finite_interval']={k:v/tops['final']['actual_component_mean']-1 for k,v in row['endmember_means'].items()}
        row['ratio_relative_to_TOPS_finite_interval']={k:v/tops['final']['ratio_upper']-1 for k,v in row['endmember_ratios'].items()}
        row['electron_density_relative_to_TOPS']=row['interpolated_electron_density_cm3']/tops['electron_density_cm3']-1
    for p,h in inputs.items():
        if digest(p)!=h:raise ValueError('input changed: '+p)
    passed=all(r['passed'] for r in checks)
    result=dict(scope=__doc__,outcome='completed_conditional_component_bounds' if passed else 'numerical_controls_failed',
        accepted_for_stellar_opacity=False,input_sha256=inputs,analytic_controls=checks,native_states=records,cases=cases,tops_control=tops,
        limitations=['Component bounds cover any absorption/scattering division of the supplied remainder only within the fixed dispersion convention.',
            'The finite frequency interval is the OP range. Missing tails and errors in source spectra or thermodynamic interpolation are not bounded by these component intervals.',
            'The interpolation bounds propagate negative Lagrange coefficients explicitly; they do not establish physical accuracy between native states.',
            'The free-electron scattering correction is held fixed as in the retained TOPS comparison. Its change in a refractive medium, bound-electron dispersion, and collision damping are not added.',
            'The four missing trace-element spectra and the limitations of the OP free-free correction remain. No stellar or atmosphere continuation is accepted.'])
    OUTPUT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    summary=dict(outcome=result['outcome'],controls=len(checks),failed=[r for r in checks if not r['passed']],native_states=len(records),
        component_interval_range=[min(r['final']['component_interval_fraction'] for r in records),max(r['final']['component_interval_fraction'] for r in records)],
        maximum_continuous_vs_native_mean=max(abs(r['continuous_vs_native_mean']) for r in records),
        cases=[dict(label=c['label'],four_point=next(r for r in c['interpolations'] if r['temperature_points']==r['density_points']==4)) for c in cases],
        TOPS_finite_ratio_change=tops['finite_interval_ratio_change'])
    print(json.dumps(summary),flush=True)
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
