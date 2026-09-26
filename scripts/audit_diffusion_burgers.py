#!/usr/bin/env python3
"""Verify the offline transport equations; synthetic collisions are not a model.

Checks include analytic binary drift, charge/baryon conservation, equilibrium
under arbitrary electron compressibility, trace limits, heat-flow coupling,
species permutations, dimensional scaling and explicit invalid-input guards.
No stellar velocities, diffusion coefficients or settling times are inferred.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from diffusion_burgers import KB, ME, MU, solve
from fetch_tops_composition import digest


def fixture(ratio=.1, ion_a=(1.,4.), ion_z=(1.,2.)):
    ions = np.array([ratio,1.])*1e25
    n = np.append(ions, np.dot(ions,ion_z))
    a, q = np.append(ion_a, ME/MU), np.append(ion_z,-1.)
    # Deliberately synthetic coefficients, with the physical n_i*n_j trace
    # scaling. No plasma collision prescription is represented by these data.
    k = np.outer(n/1e25,n/1e25)*1e9
    zz = np.full((3,3),.6)
    return dict(number_density=n, mass_u=a, charge=q,
                ion_pressure_gradient=ions*np.array([-.8, -2.])*MU*2e5,
                gravity=2e5, temperature=1e7, temperature_gradient=-.03,
                resistance=k, z=zz, zprime=np.full((3,3),1.3),
                zdoubleprime=np.full((3,3),2.), heat_flow='classical')


def independent_balances(p, result):
    """Evaluate the original pair-force equations, not the assembled matrix."""
    n,a,q = [np.asarray(p[k]) for k in ('number_density','mass_u','charge')]
    w,h = result.velocity_cm_s, result.residual_heat_velocity_cm_s
    k,z,zp,zpp = [np.asarray(p[v]) for v in ('resistance','z','zprime','zdoubleprime')]
    extra = p.get('ion_extra_force',np.zeros(len(n)-1))
    residual, scales = [], []
    for s in range(len(n)-1):
        terms = [p['ion_pressure_gradient'][s], n[s]*a[s]*MU*p['gravity'],
                 -n[s]*q[s]*result.electric_force_dyn, -n[s]*extra[s]]
        for t in range(len(n)):
            if s == t:
                continue
            terms += [-k[s,t]*(w[t]-w[s]),
                      -k[s,t]*z[s,t]*(a[t]*h[s]-a[s]*h[t])/(a[s]+a[t])]
        residual.append(abs(sum(terms)))
        scales.append(sum(abs(v) for v in terms))
    if p['heat_flow'] == 'classical':
        for s in range(len(n)):
            terms = [2.5*n[s]*KB*p['temperature_gradient'], .4*k[s,s]*zpp[s,s]*h[s]]
            for t in range(len(n)):
                if s == t:
                    continue
                m_s,m_t = a[s],a[t]
                terms += [2.5*k[s,t]*z[s,t]*m_t/(m_s+m_t)*(w[t]-w[s]),
                          k[s,t]*((3*m_s*m_s+m_t*m_t*zp[s,t])
                                  +.8*m_s*m_t*zpp[s,t])*h[s]/(m_s+m_t)**2,
                          -k[s,t]*m_s*m_t*(3+zp[s,t]-.8*zpp[s,t])*h[t]/(m_s+m_t)**2]
            residual.append(abs(sum(terms)))
            scales.append(sum(abs(v) for v in terms))
    speed = max(np.max(np.abs(w)),1e-100)
    return dict(pair_equation_relative_residual=max(residual)/max(max(scales),1e-100),
                baryon_flux_relative_residual=abs(np.dot(n[:-1]*a[:-1],w[:-1]))/(np.dot(n[:-1],a[:-1])*speed),
                current_relative_residual=abs(np.dot(n*q,w))/(np.dot(n,np.abs(q))*speed))


def analytic_binary(p):
    n,a,q = [np.asarray(p[k]) for k in ('number_density','mass_u','charge')]
    k = p['resistance']
    alpha = n[0]*a[0]/(n[1]*a[1])
    beta = (n[0]*q[0]-n[1]*q[1]*alpha)/n[2]
    c0 = (-k[0,1]*(1+alpha)+k[0,2]*(beta-1))/n[0]
    c1 = (k[1,0]*(1+alpha)+k[1,2]*(beta+alpha))/n[1]
    b = p['ion_pressure_gradient']/n[:-1]+a[:-1]*MU*p['gravity']
    w0 = (b[0]*q[1]-b[1]*q[0])/(c0*q[1]-c1*q[0])
    electric = (b[0]-c0*w0)/q[0]
    return np.array([w0,-alpha*w0,beta*w0]),electric


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError('preserve completed transport reports')
    inputs={str(p.resolve()):digest(p) for p in
            (Path(__file__),Path('scripts/diffusion_burgers.py'))}
    checks=[]
    def check(name, error, tolerance=1e-10, **metadata):
        if not np.isfinite(error) or error > tolerance:
            raise ValueError(f'{name}: {error} exceeds {tolerance}')
        checks.append(dict(name=name,error=float(error),tolerance=tolerance,**metadata))
    def compare(name, actual, expected, scale=None, **kw):
        actual,expected=np.asarray(actual),np.asarray(expected)
        denominator=np.max(np.abs(expected)) if scale is None else scale
        check(name,np.max(np.abs(actual-expected))/max(denominator,1e-100),**kw)

    records=[]
    for ratio in (1.,.1,1e-4,1e-8,1e-12,1e-16):
        p=fixture(ratio)
        p['z']=np.zeros((3,3))
        expected,field=analytic_binary(p)
        for mode in ('classical','suppressed'):
            p['heat_flow']=mode
            result=solve(**p)
            compare(f'analytic binary velocities {ratio} {mode}',result.velocity_cm_s,expected)
            compare(f'analytic binary electric force {ratio} {mode}',result.electric_force_dyn,field)
            balance=independent_balances(p,result)
            for label,value in balance.items():
                check(f'binary {ratio} {mode} {label}',value)
            records.append(dict(test='binary',number_ratio=ratio,heat_flow=mode,
                                scaled_condition_number=result.scaled_condition_number,**balance))

    # Arbitrary-compressibility equilibria derived independently by combining
    # charge neutrality, ideal-ion force balance and dP_e/dln(n_e).
    # B = (dP_e/dln n_e)/(n_e kT); B=1 is classical, B>>1 is degenerate.
    equilibrium=[]
    for ratio in (1.,1e-4,1e-12):
        for b in (1.,2.,10.,1e4):
            p=fixture(ratio)
            p['temperature_gradient']=0.
            n,a,q=[p[k] for k in ('number_density','mass_u','charge')]
            field=MU*p['gravity']*np.dot(n[:-1],q[:-1]*a[:-1])/(np.dot(n[:-1],q[:-1]**2)+n[-1]/b)
            p['ion_pressure_gradient']=n[:-1]*(q[:-1]*field-a[:-1]*MU*p['gravity'])
            # Independent electron momentum plus derivative of neutrality.
            dne_dr=np.dot(q[:-1],p['ion_pressure_gradient'])/(KB*p['temperature'])
            pe_gradient=b*KB*p['temperature']*dne_dr
            check('equilibrium electron force',abs(pe_gradient+n[-1]*field)/(n[-1]*abs(field)),tolerance=1e-8)
            check('equilibrium total hydrostatic force',abs(np.sum(p['ion_pressure_gradient'])+pe_gradient
                  +np.dot(n[:-1],a[:-1])*MU*p['gravity'])/(np.dot(n[:-1],a[:-1])*MU*p['gravity']),tolerance=1e-8)
            for mode in ('classical','suppressed'):
                p['heat_flow']=mode
                result=solve(**p)
                # Drift scale from one particle's gravitational force/drag.
                reference=MU*p['gravity']/(np.max(p['resistance']/n[:,None]))
                compare('zero equilibrium drift',result.velocity_cm_s,np.zeros(3),scale=reference)
                compare('zero equilibrium residual heat flow',result.residual_heat_velocity_cm_s,np.zeros(3),scale=reference)
                compare('equilibrium electric field',result.electric_force_dyn,field)
            equilibrium.append(dict(hydrogen_helium_number_ratio=ratio,electron_stiffness_over_nkT=b,
                                    electric_force_over_mu_g=field/(MU*p['gravity'])))
    # In a pure helium isothermal background the analytic limits are 4/3 and 2.
    classical,degenerate=equilibrium[-4],equilibrium[-1]
    check('helium classical electric limit',abs(classical['electric_force_over_mu_g']/(4/3)-1))
    check('helium degenerate electric limit',abs(degenerate['electric_force_over_mu_g']/2-1),tolerance=5.1e-5)

    p=fixture(.2)
    base=solve(**p)
    for label,value in independent_balances(p,base).items():
        check('coupled heat '+label,value)
    # Nonzero heat-flow coupling must actually affect this numerical fixture.
    suppressed=solve(**{**p,'heat_flow':'suppressed'})
    relative_heat_effect=float(np.linalg.norm(base.velocity_cm_s-suppressed.velocity_cm_s)/np.linalg.norm(base.velocity_cm_s))
    if relative_heat_effect < 1e-3:
        raise ValueError('fixture does not exercise heat-flow coupling')
    for factor in (1e-18,1e-6,1e6,1e18):
        changed={**p,'resistance':p['resistance']*factor}
        result=solve(**changed)
        compare('inverse collision-rate velocity scaling',result.velocity_cm_s*factor,base.velocity_cm_s)
        compare('inverse collision-rate heat scaling',result.residual_heat_velocity_cm_s*factor,base.residual_heat_velocity_cm_s)
        compare('collision-independent field scaling',result.electric_force_dyn,base.electric_force_dyn)
        changed={**p,**{key:p[key]*factor for key in
                       ('number_density','ion_pressure_gradient','resistance')}}
        result=solve(**changed)
        compare('equation density-unit scaling',result.velocity_cm_s,base.velocity_cm_s)
        compare('equation density-unit electric scaling',result.electric_force_dyn,base.electric_force_dyn)
    permutation=np.array([1,0,2])
    changed={**p,**{key:p[key][permutation] for key in ('number_density','mass_u','charge')},
             'ion_pressure_gradient':p['ion_pressure_gradient'][[1,0]],
             **{key:p[key][np.ix_(permutation,permutation)] for key in
                ('resistance','z','zprime','zdoubleprime')}}
    result=solve(**changed)
    compare('ion permutation velocity',result.velocity_cm_s,base.velocity_cm_s[permutation])
    compare('ion permutation heat',result.residual_heat_velocity_cm_s,base.residual_heat_velocity_cm_s[permutation])
    compare('ion permutation electric field',result.electric_force_dyn,base.electric_force_dyn)
    reversed_forces={**p,'ion_pressure_gradient':-p['ion_pressure_gradient'],
                     'temperature_gradient':-p['temperature_gradient'],
                     'ion_extra_force':2*p['mass_u'][:-1]*MU*p['gravity']}
    reversed_result=solve(**reversed_forces)
    compare('force reversal velocities',reversed_result.velocity_cm_s,-base.velocity_cm_s)
    compare('force reversal heat',reversed_result.residual_heat_velocity_cm_s,-base.residual_heat_velocity_cm_s)
    compare('force reversal electric force',reversed_result.electric_force_dyn,-base.electric_force_dyn)

    identical=fixture(.2,ion_a=(4.,4.),ion_z=(2.,2.))
    identical['temperature_gradient']=0.
    identical['ion_pressure_gradient']=identical['number_density'][:-1]*(-2*MU*identical['gravity'])
    identical_result=solve(**identical)
    compare('equal-species no separation',identical_result.velocity_cm_s,np.zeros(3),scale=np.max(np.abs(base.velocity_cm_s)))

    guards=[]
    def reject(name,override):
        try:
            solve(**{**p,**override})
        except ValueError as error:
            guards.append(dict(name=name,reason=str(error)))
        else:
            raise ValueError(f'invalid input accepted: {name}')
    reject('negative temperature',{'temperature':-1.})
    reject('nonfinite force',{'gravity':float('nan')})
    reject('absent species',{'number_density':np.array([0.,1e25,2e25])})
    reject('non-neutral plasma',{'charge':np.array([1.,2.,-2.])})
    reject('missing heat approximation',{'heat_flow':'auto'})
    reject('disconnected populations',{'resistance':np.eye(3)})
    asymmetric=p['resistance'].copy();asymmetric[0,1]*=2
    reject('asymmetric resistance',{'resistance':asymmetric})
    reject('nonfinite coefficient',{'z':np.full((3,3),float('nan'))})
    for name,sha in inputs.items():
        if digest(Path(name))!=sha:
            raise ValueError('transport implementation changed during verification')
    report=dict(scope=__doc__,status='pass',accepted_for_stellar_evolution=False,
                physical_collision_coefficients_supplied=False,numpy_version=np.__version__,
                barycentric_frame='zero ion baryonic mass flux; electron mass retained only in collision kinematics',
                checks=checks,check_count=len(checks),invalid_input_controls=guards,
                binary_controls=records,equilibrium_controls=equilibrium,
                synthetic_fixture_relative_heat_flow_effect=relative_heat_effect,
                largest_scaled_linear_condition_number=max(r['scaled_condition_number'] for r in records),
                limitations=['Synthetic collision matrices test the equations, not physical diffusion rates.',
                             'No stellar velocities, time-dependent abundances or settling times were calculated.',
                             'The classical heat equations are not valid for degenerate electrons.',
                             'Suppressed residual heat flow is an explicit comparison, not an accepted stellar prescription.',
                             'Ionization, nonideal ion forces, screened collision integrals and energy coupling remain required.'],
                references={'equations_and_equilibrium_limits':'https://arxiv.org/html/1710.08424v2#S3',
                            'classical_multicomponent_transport':'https://arxiv.org/abs/astro-ph/9304005'},
                input_sha256=inputs)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({key:report[key] for key in ('status','check_count','accepted_for_stellar_evolution',
                      'synthetic_fixture_relative_heat_flow_effect','largest_scaled_linear_condition_number')}))


if __name__=='__main__':
    main()
