#!/usr/bin/env python3
"""Check analytic electron-energy elimination against constrained energy integrals.

Reuse retained stellar collision moments. Independent momentum quadrature
evaluates the underlying sum of species collision integrals and its current
constraints, rather than only reassembling the reduced matrix.
"""
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit

from electron_elastic_mixture import elastic_mixture_form
from electron_ion_heat import current_basis, physical_prefactor
from electron_ion_born import coulomb_bracket
from diffusion_composition_forces import velocity_basis
from ion_collision_integrals import MU
from fetch_tops_composition import digest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('preserve completed mixture checks')
    files=[Path(__file__),Path('scripts/electron_elastic_mixture.py'),
           Path('scripts/electron_ion_heat.py'),Path('scripts/diffusion_composition_forces.py'),
           Path('docs/results/stellar_electron_heat_3890gyr_v1.json'),
           Path('docs/results/diffusion_material_regime_3890gyr_v1.json'),
           Path('include/ember/gs98_mixture.hpp')]
    inputs={str(p.resolve()):digest(p) for p in files}
    profile=json.loads(files[4].read_text());regime=json.loads(files[5].read_text())
    if profile['outcome']!='completed_model_comparison':raise ValueError('profile comparison did not complete')
    for k,h in profile['inputs_sha256'].items():
        if digest(Path(k))!=h:raise ValueError('retained source differs')
        inputs[k]=h
    metal=[tuple(map(float,s.split(','))) for s in re.findall(r'^\s*\{([^{}]+)\}, //',files[6].read_text(),re.M)]
    mass=np.array([1.,3.,4.]+[r[1] for r in metal]);charge=np.array([1.,2.,2.]+[r[0] for r in metal])
    checks=[];cases=[];guard_count=0
    def check(name,error,tolerance,**where):
        checks.append(dict(name=name,normalized_error=float(error),tolerance=tolerance,passed=bool(np.isfinite(error) and error<=tolerance),**where))
    gauss={n:leggauss(n) for n in [128,256]}
    for index in [0,100,234,300,396,421]:
        row=profile['records'][index];state=regime['records'][index]
        x=np.array([state['X'],state['Y3'],state['Y4']]+[(1-state['X']-state['Y3']-state['Y4'])*m[3] for m in metal])
        n=state['density_baryonic_g_cm3']/MU*x/mass
        b=velocity_basis(mass_fractions=x,mass_numbers=mass,charges=charge,
                         independent_ions=[0,1],reference_ion=2)
        basis=current_basis(row['eta_nonrelativistic']);eta=basis.eta
        for case in row['cases']:
            I=np.array(case['dimensionless_collision_moments']);L=np.array(case['dimensionless_full_energy_response'])
            kw=dict(velocity_map=b,ion_density=n,ion_charges=charge,collision_moments=I,full_energy_response=L)
            trial=elastic_mixture_form(**kw,relaxation='one_heat_variable')
            full=elastic_mixture_form(**kw)
            where=dict(zone=index,screening=case['screening'])
            scale=np.sqrt(np.diag(trial.collision_matrix));normal=np.outer(scale,scale)
            unit=physical_prefactor(1.,1.);weights=n*charge**2
            # Independent pair assembly for the restricted two-variable form.
            pair=np.zeros((3,3))
            for i,weight in enumerate(weights):
                A=np.array([list(b[-1]-b[i])+[0.],[0.,0.,1.]])
                pair+=unit*weight*(A.T@I@A)
            check('independent pair assembly',np.max(np.abs((pair-trial.collision_matrix)/normal)),2e-12,**where)
            check('relaxation cannot increase collision cost',max(0.,-np.linalg.eigvalsh((trial.collision_matrix-full.collision_matrix)/normal)[0]),2e-11,**where)
            check('positive constrained electron cost',max(0.,-np.linalg.eigvalsh(full.collision_matrix/normal)[0]),2e-12,**where)
            permutation=np.arange(len(n))[::-1]
            permuted=elastic_mixture_form(**{**kw,'ion_density':n[permutation],'ion_charges':charge[permutation],
                                            'velocity_map':np.vstack([b[:-1][permutation],b[-1]])})
            check('ion-order invariance',np.max(np.abs((permuted.collision_matrix-full.collision_matrix)/normal)),2e-12,**where)
            shifted=elastic_mixture_form(**{**kw,'velocity_map':b+np.array([.13,-.27])})
            check('common-velocity invariance',np.max(np.abs((shifted.collision_matrix-full.collision_matrix)/normal)),2e-12,**where)
            scaled=elastic_mixture_form(**{**kw,'ion_density':n*4})
            check('scatterer-density normalization',np.max(np.abs((scaled.collision_matrix/4-full.collision_matrix)/normal)),2e-12,**where)
            quadratures={}
            for order in [128,256]:
                # Integrate in momentum t=sqrt(epsilon/kT), with independent
                # nodes and weights on both sides of the Fermi surface.
                nodes,weight=gauss[order];edges=sorted(set([0.,math.sqrt(max(eta,0)+70)]+
                    [math.sqrt(eta+y) for y in [-20.,-5.,0.,5.,20.] if eta+y>0]))
                ts=[];ws=[]
                for lo,hi in zip(edges[:-1],edges[1:],strict=True):
                    ts.extend((lo+hi)/2+(hi-lo)/2*nodes);ws.extend((hi-lo)/2*weight)
                t=np.array(ts);w=np.array(ws);energy=t*t;s=energy-eta-basis.centered_mean
                occupation=expit(eta-energy)*expit(energy-eta)
                kernel=coulomb_bracket(case['b_thermal']*energy)
                current=2*t**4*occupation*w/basis.normalization
                collision=2*t*kernel*occupation*w
                quadratures[order]=(t,s,kernel,current,collision)
                check('independent current normalization',abs(np.sum(current)-1),2e-11,order=order,**where)
                check('independent centered heat moment',abs(np.sum(current*s))/math.sqrt(basis.variance),2e-11,order=order,**where)
            # Basis and mixed directions, scaled so trace-species columns are
            # tested on the same entropy scale as abundant species.
            directions=[*np.eye(3),[1.,1.,0.],[1.,0.,1.],[0.,1.,1.],[.2,-.7,.4]]
            worst=0.
            for number,direction in enumerate(directions):
                variable=np.asarray(direction)/scale
                ui=b[:-1]@variable[:-1];ue=float(b[-1]@variable[:-1]);heat=variable[-1]
                ubar=float(weights@ui/np.sum(weights))
                lagrange=np.linalg.solve(L,np.array([ue-ubar,heat]))
                expected=float(variable@full.collision_matrix@variable)
                costs=[]
                for order,(t,s,kernel,current,collision) in quadratures.items():
                    u=ubar+t**3/kernel/basis.normalization*(lagrange[0]+lagrange[1]*s)
                    measured_u=float(current@u);measured_heat=float((current*s)@u)
                    uscale=max(abs(ue),abs(ubar),abs(heat),np.max(np.abs(ui)),np.finfo(float).tiny)
                    check('energy solution satisfies prescribed electron current',abs(measured_u-ue)/uscale,3e-10,direction=number,order=order,**where)
                    check('energy solution satisfies prescribed electron heat',abs(measured_heat-heat)/uscale,3e-10,direction=number,order=order,**where)
                    cost=unit*sum(weight*float(collision@((u-ion_u)**2)) for weight,ion_u in zip(weights,ui,strict=True))
                    costs.append(cost);error=abs(cost/expected-1);worst=max(worst,error)
                    check('independent full-energy collision integral',error,3e-10,direction=number,order=order,**where)
                check('momentum quadrature resolution',abs(costs[1]/costs[0]-1),3e-10,direction=number,**where)
            cases.append(dict(**where,maximum_energy_cost_relative_difference=worst,
                              minimum_scaled_cost_eigenvalue=float(np.linalg.eigvalsh(full.collision_matrix/normal)[0])))
    # One scatterer population has no relative-ion variance contribution.
    single=elastic_mixture_form(velocity_map=[[0.],[1.]],ion_density=[1e20],ion_charges=[2.],
                                collision_moments=I,full_energy_response=L)
    expect=physical_prefactor(1e20,2.)*np.linalg.inv(L)
    check('single-ion exact elastic response',np.max(np.abs((single.collision_matrix-expect)/np.sqrt(np.outer(np.diag(expect),np.diag(expect))))),2e-12)
    check('single-ion variance vanishes',np.max(np.abs(single.scatterer_variance_matrix)),0.)
    for change in [dict(ion_density=np.zeros_like(n)),dict(ion_charges=-charge),
                   dict(velocity_map=np.full_like(b,np.nan)),dict(collision_moments=[[1.,2.],[2.,1.]]),
                   dict(full_energy_response=[[1.,0.],[1.,1.]]),dict(relaxation='unspecified')]:
        failed=False
        try:elastic_mixture_form(**{**kw,**change})
        except ValueError:failed=True
        guard_count+=1;check('invalid-input guard',0. if failed else 1.,0.,guard=guard_count)
    for name,h in inputs.items():
        if digest(Path(name))!=h:raise ValueError('input changed during the audit')
    result=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),
                outcome='passed' if all(c['passed'] for c in checks) else 'failed',
                accepted_for_stellar_evolution=False,checks=checks,cases=cases,guards=guard_count,
                maximum_energy_cost_relative_difference=max(c['maximum_energy_cost_relative_difference'] for c in cases),
                new_scattering_calculations=0,inputs_sha256=inputs,
                limitations=['Common-shape screened Born collisions only; fixed-ion elastic approximation.',
                             'No electron-electron collisions, ion heat coupling, recoil, correlations or stellar flux.'])
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['outcome','guards','maximum_energy_cost_relative_difference','new_scattering_calculations']}),flush=True)
    print('checks',len(checks),'failed',sum(not c['passed'] for c in checks),flush=True)
    if result['outcome']!='passed':raise SystemExit(1)


if __name__=='__main__':main()
