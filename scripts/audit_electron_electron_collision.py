#!/usr/bin/env python3
"""Bounded independent normalization and saved-state electron collision audit."""
import concurrent.futures
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.integrate import quad
from scipy.special import expit
from electron_electron_collision import collision_matrix, energy_basis, born_electron_ion_moments
from electron_ion_heat import physical_prefactor
from electron_ion_born import KB, ME, HBAR, E2


def classical_heat_reference(eta, bthermal, *, exchange, constant=False):
    """Independent dimensional Burgers Omega22 integral; no pair-energy rule."""
    temp=1e7;energy=KB*temp
    density=2/(2*math.pi*HBAR)**3*(2*ME*energy)**1.5*math.pi**1.5*math.exp(eta)
    amplitude2=(E2/(2*energy))**2
    def viscosity_cross(gamma):
        if constant:return 8*math.pi/3*amplitude2
        def integral(y):
            direct=1/(4/bthermal+2*gamma*gamma*y)
            other=1/(4/bthermal+2*gamma*gamma*(1-y))
            cross=direct*direct+other*other-(direct*other if exchange else 0)
            return 4*y*(1-y)*cross
        return amplitude2*4*math.pi*quad(integral,0,.5,epsabs=1e-100,epsrel=3e-11)[0]
    omega22=math.sqrt(energy/(math.pi*ME))*quad(
        lambda gamma:math.exp(-gamma*gamma)*gamma**7*viscosity_cross(gamma),
        0,12,epsabs=1e-100,epsrel=3e-10)[0]
    resistance=.16*(16/3)*density*density*(ME/2)*omega22
    return resistance/physical_prefactor(density,1)


def execute(task):
    name,args=task
    start=time.monotonic()
    result=collision_matrix(**args)
    result['elapsed_seconds']=time.monotonic()-start
    return name,result


def main():
    work=Path('/tmp/ember-electron-electron-collision-run-v1')
    report=Path('docs/results/electron_electron_collision_v1.json')
    assert work.is_dir() and not report.exists()
    checks=[];rows={}
    def check(name,error,tolerance,**meta):
        checks.append(dict(name=name,error=float(error),tolerance=tolerance,
                           passed=bool(np.isfinite(error) and abs(error)<=tolerance),**meta))
    rng=np.random.default_rng(1977)
    # Compare the reduced dot products with direct changes of all four vectors.
    for eta in [-16.,0.,6.563364251342554,20.]:
        basis=energy_basis(eta)
        for i in range(64):
            G=rng.normal(size=3)*math.sqrt(max(eta,1)/3)
            g=rng.normal(size=3)*2;gp=rng.normal(size=3);gp*=np.linalg.norm(g)/np.linalg.norm(gp)
            before=np.array([G+g/2,G-g/2]);after=np.array([G+gp/2,G-gp/2])
            x=np.sum(before*before,axis=1);xp=np.sum(after*after,axis=1)
            delta=basis.polynomials(x)@before-basis.polynomials(xp)@after
            b,bp=G@g,G@gp;d=b*g-bp*gp;h=b*b-bp*bp
            z=2*(G@G+g@g/4-eta-basis.mean)-basis.skew_ratio
            expected=np.array([np.zeros(3),d/basis.variance,
                               (z*d+2*h*G)/basis.quadratic_norm])
            scale=1+np.max(abs(delta))
            check('direct four-vector polynomial change',np.max(abs(delta-expected))/scale,2e-12,eta=eta,sample=i)
            f11=d@d/basis.variance**2
            f12=(z*(d@d)+2*h*h)/(basis.variance*basis.quadratic_norm)
            f22=(z*z*(d@d)+4*(z+G@G)*h*h)/basis.quadratic_norm**2
            algebra=np.array([[f11,f12],[f12,f22]])
            check('reduced versus direct collision dot products',
                  np.max(abs(algebra-delta[1:]@delta[1:].T))/(1+np.max(abs(algebra))),3e-12,eta=eta,sample=i)
            pre=np.prod(expit(eta-x))*np.prod(expit(xp-eta))
            post=np.prod(expit(eta-xp))*np.prod(expit(x-eta))
            check('Pauli detailed balance',(pre-post)/max(pre,post,1e-290),2e-12,eta=eta,sample=i)
            a=G@G+g@g/4-eta
            cosh=1/(4*(math.cosh(a)+math.cosh(b))*(math.cosh(a)+math.cosh(bp)))
            check('independent four-occupation weight',(pre-cosh)/max(pre,cosh,1e-290),2e-12,eta=eta,sample=i)
    source=Path('docs/results/stellar_electron_heat_3890gyr_v1.json')
    saved=json.loads(source.read_text())
    selected=[saved['records'][i] for i in [0,300,396]]
    low=(40,20,20,24,24);high=(64,32,32,40,40)
    tasks=[]
    for r in selected:
        eta=r['eta_nonrelativistic'];b=r['cases'][0]['b_thermal']
        for name,orders in [('low',low),('high',high)]:
            tasks.append((f"zone_{r['zone']}_{name}",dict(eta=eta,b_thermal=b,orders=orders)))
        got=born_electron_ion_moments(energy_basis(eta),b)
        old=np.asarray(r['cases'][0]['dimensionless_collision_moments'])
        check('retained two-mode electron-ion moments',np.max(abs(got[:2,:2]-old))/np.max(abs(old)),2e-10,zone=r['zone'])
    controls=[('maxwell_constant',dict(eta=-16.,b_thermal=3.,statistics='maxwell',exchange=False,scattering='constant')),
              ('maxwell_born',dict(eta=-16.,b_thermal=3.,statistics='maxwell',exchange=True)),
              ('fermi_dilute',dict(eta=-16.,b_thermal=3.,statistics='fermi',exchange=True))]
    for name,args in controls:tasks.append((name,{**args,'orders':(48,24,16,40,24)}))
    # Each independent state is computed once; completed outputs are retained
    # separately so a failed comparison need not repeat any collision integral.
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(execute,t):t[0] for t in tasks}
        for future in concurrent.futures.as_completed(futures):
            name,result=future.result();rows[name]=result
            (work/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n')
            print(name,'completed',flush=True)
    for name,result in rows.items():
        matrix=np.asarray(result['collision_matrix'])
        check('momentum null mode',np.max(abs(matrix[0])),0,case=name)
        check('reciprocity',np.max(abs(matrix-matrix.T)),0,case=name)
        eigen=np.linalg.eigvalsh(matrix[1:,1:])
        check('positive heat collision matrix',max(0.,-float(eigen[0])),0,case=name)
    for name,constant in [('maxwell_constant',True),('maxwell_born',False)]:
        reference=classical_heat_reference(-16.,3.,exchange=not constant,constant=constant)
        actual=rows[name]['collision_matrix'][1][1]
        check('independent dimensional classical heat integral',actual/reference-1,2e-7,case=name,reference=reference)
    exact_constant=8/(2.5**2*math.sqrt(2))*math.exp(-16.)
    check('analytic constant cross-section normalization',rows['maxwell_constant']['collision_matrix'][1][1]/exact_constant-1,2e-7)
    dilute=np.asarray(rows['fermi_dilute']['collision_matrix'])[1:,1:]
    maxwell=np.asarray(rows['maxwell_born']['collision_matrix'])[1:,1:]
    check('dilute Fermi to classical matrix',np.max(abs(dilute-maxwell))/np.max(abs(maxwell)),2e-6)
    physics=[]
    for r in selected:
        i=r['zone'];lo=np.asarray(rows[f'zone_{i}_low']['collision_matrix'])
        hi=np.asarray(rows[f'zone_{i}_high']['collision_matrix'])
        scale=np.sqrt(np.diag(hi)[1:]);difference=np.max(abs((lo-hi)[1:,1:]/scale[:,None]/scale[None,:]))
        check('five-dimensional quadrature refinement',difference,.005,zone=i)
        eta=r['eta_nonrelativistic'];b=r['cases'][0]['b_thermal']
        ei=born_electron_ion_moments(energy_basis(eta),b)
        # These first comparisons use fixed pure-H scatterers: n_i=n_e. They
        # evaluate the saved eta and screening, not the star's mixture fluxes.
        heat_noee=np.linalg.inv(ei[1:,1:])[0,0]
        heat_ee=np.linalg.inv((ei+hi)[1:,1:])[0,0]
        heat_two=1/(ei[1,1]+hi[1,1])
        physics.append(dict(zone=i,eta=eta,b_thermal=b,
            comparison='fixed pure-H scatterers at the saved eta and screening; not a stellar flux',
            three_mode_electron_ion_matrix=ei.tolist(),electron_electron_matrix=hi.tolist(),
            quadrature_relative_change=float(difference),
            ee_to_ei_heat_diagonal=float(hi[1,1]/ei[1,1]),
            heat_response_with_ee_over_without=float(heat_ee/heat_noee),
            two_mode_over_three_mode_with_ee=float(heat_two/heat_ee)))
    inputs=[Path(__file__),Path('scripts/electron_electron_collision.py'),
            Path('scripts/electron_ion_heat.py'),Path('scripts/electron_ion_born.py'),source]
    output=dict(created_utc=datetime.now(timezone.utc).isoformat(),
        scope='Born/Pauli electron collision brackets, three energy modes; independent classical normalization and saved-eta comparisons',
        outcome='passed' if all(c['passed'] for c in checks) else 'failed',
        accepted_for_stellar_evolution=False,checks=checks,calculations=rows,physical_comparisons=physics,
        inputs_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        limitations=['No higher-mode convergence beyond quadratic energy polynomials.',
                     'No full quantum phase-shift, dynamic screening, mixture correlations or relativity.',
                     'No strongly degenerate quadrature validation or stellar mixture coupling.',
                     'Pure-H fixed-scatterer heat comparisons do not calculate stellar diffusion.'])
    report.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(dict(outcome=output['outcome'],checks=len(checks),failed=[c for c in checks if not c['passed']],physical_comparisons=physics)),flush=True)
    if output['outcome']!='passed':raise SystemExit(1)


if __name__=='__main__':main()
