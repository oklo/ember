#!/usr/bin/env python3
"""Independent current, collision, classical and degenerate electron heat checks.

Momentum-space Gauss quadrature provides an independent integration path.
Saved displaced-distribution drag controls are reused. Analytic Lorentz and
constant-frequency cases check both the two-mode and exact elastic response.
Differences between those responses are measured model truncation, not failures
to be hidden by changing tolerances. No stellar transport is accepted here.
"""
import argparse,hashlib,json,math
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.special import expit,gamma
from electron_ion_born import KB,ME,HBAR,fermi_half,chemical_potential
from electron_ion_heat import current_basis,born_response,collision_response,physical_prefactor


def direct_momentum(eta,b,order=128):
    """Independent t=p/sqrt(2 m kT) quadrature, including a direct log formula."""
    nodes,weights=leggauss(order)
    low=max(0.,eta-60);high=max(eta,0)+60
    breaks=sorted(set([math.sqrt(low),math.sqrt(high)]+[math.sqrt(eta+q) for q in
        [-40,-20,-5,0,5,20,40] if low<eta+q<high]))
    ts=[];ws=[]
    for left,right in zip(breaks,breaks[1:]):
        ts.extend((left+right)/2+(right-left)/2*nodes);ws.extend(weights*(right-left)/2)
    t=np.array(ts);dx_weights=np.array(ws)*2*t;x=t*t;y=x-eta
    occupation=expit(y)*expit(-y)
    current=dx_weights*x**1.5*occupation;norm=sum(current)
    mean=current@y/norm;variance=current@(y-mean)**2/norm
    # longdouble avoids the small-b cancellation in this independent formula.
    bq=np.asarray(b*x,dtype=np.longdouble)
    bracket=np.asarray(np.log1p(bq)-bq/(1+bq),dtype=float)
    p=np.array([np.ones_like(y),(y-mean)/variance])
    moment=np.einsum('an,n,bn->ab',p,dx_weights*bracket*occupation,p)
    s=np.array([np.ones_like(y),y-mean])
    inverse=np.einsum('an,n,bn->ab',s,dx_weights*x**3/bracket*occupation,s)/norm**2
    return norm,mean,variance,moment,inverse


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('report',type=Path);a=parser.parse_args()
    checks=[];comparisons=[];cache={}
    def check(name,actual,reference,scale,tolerance=2e-8,**meta):
        error=float(np.max(abs(np.asarray(actual)-reference)/scale))
        checks.append(dict(name=name,normalized_error=error,tolerance=tolerance,
                           passed=bool(np.isfinite(error) and error<=tolerance),**meta))
    def norm(matrix):
        d=np.sqrt(np.diag(matrix));return d[:,None]*d[None,:]
    for eta in [-25.,-8.,-2.,0.,6.,20.,100.,1e4]:
        basis=current_basis(eta)
        check('number-current normalization',basis.normalization,1.5*fermi_half(eta),basis.normalization,2e-10,eta=eta)
        # Bulk kinetic enthalpy from occupation, not its derivative or current basis.
        upper=math.sqrt(max(eta,0)+60);points=[math.sqrt(eta)] if eta>0 else None
        f32=quad(lambda t:2*t**4*expit(eta-t*t),0,upper,points=points,epsabs=1e-100,epsrel=2e-12,limit=200)[0]
        enthalpy=(5/3)*f32/fermi_half(eta)
        check('independent ideal electron enthalpy',eta+basis.centered_mean,enthalpy,enthalpy,2e-10,eta=eta)
        for b in [.05,.1,1.,10.,100.]:
            result=born_response(eta,b,basis=basis);cache[eta,b]=result
            m=np.array(result['collision_moments']);exact=np.array(result['exact_elastic_mobility'])
            direct=direct_momentum(eta,b);fine=direct_momentum(eta,b,256)
            meta=dict(eta=eta,b_thermal=b)
            check('independent current normalization',basis.normalization,fine[0],fine[0],2e-9,**meta)
            check('independent centered current mean',basis.centered_mean,fine[1],1+abs(fine[1]),2e-9,**meta)
            check('independent current variance',basis.variance,fine[2],fine[2],2e-9,**meta)
            check('independent momentum-space collision matrix',m,fine[3],norm(fine[3]),2e-8,**meta)
            check('independent exact elastic response',exact,fine[4],norm(fine[4]),2e-8,**meta)
            check('momentum-space collision quadrature order',direct[3],fine[3],norm(fine[3]),2e-8,**meta)
            check('momentum-space exact-response quadrature order',direct[4],fine[4],norm(fine[4]),2e-8,**meta)
            trial=np.array(result['two_mode_mobility']);loss=(exact-trial)/norm(exact)
            check('variational response does not exceed exact response',max(0.,-min(np.linalg.eigvalsh(loss))),0.,1.,2e-8,**meta)
            check('two-mode heat bound',max(0.,result['two_mode_to_exact_zero_drift_heat']-1),0.,1.,2e-8,**meta)
            comparisons.append(dict(**meta,zero_drift_heat_ratio=result['two_mode_to_exact_zero_drift_heat'],
                isothermal_drift_ratio=float(trial[0,0]/exact[0,0]),
                collision_cross_over_drag=float(m[0,1]/m[0,0]),
                collision_heat_over_drag=float(m[1,1]/m[0,0])))
        constant=collision_response(eta,lambda x:x**1.5,basis=basis)
        expected=np.diag([basis.normalization,basis.normalization/basis.variance])
        check('constant-frequency collision oracle',constant['collision_moments'],expected,norm(expected),2e-10,eta=eta)
        inv=np.diag([1/basis.normalization,basis.variance/basis.normalization])
        check('constant-frequency exact response',constant['exact_elastic_mobility'],inv,norm(inv),2e-10,eta=eta)
        check('constant-frequency two-mode completeness',constant['two_mode_mobility'],inv,norm(inv),2e-10,eta=eta)
    # The Maxwellian constant-log Lorentz kernel has exact Gamma-function moments.
    eta=-25.;lorentz=collision_response(eta,lambda x:1.)
    expected=math.exp(eta)*np.array([[1.,-.6],[-.6,.52]])
    check('classical Coulomb collision moments',lorentz['collision_moments'],expected,norm(expected),2e-9)
    exact=math.exp(-eta)/gamma(2.5)**2*np.array([[6.,9.],[9.,37.5]])
    check('classical exact Lorentz heat and drift response',lorentz['exact_elastic_mobility'],exact,norm(exact),2e-9)
    classical=current_basis(eta)
    check('classical current mean',eta+classical.centered_mean,2.5,2.5,2e-9)
    check('classical current heat variance',classical.variance,2.5,2.5,2e-9)
    # Reuse physical finite-displacement drag references from the preceding audit.
    old_path=Path('docs/results/electron_ion_born_v1.json');old=json.loads(old_path.read_text())
    assert old['status']=='pass'
    for path,h in old['input_sha256'].items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==h
    for control in old['checks']:
        if control['name']=='displaced_distribution':
            eta,b=control['eta'],control['b_thermal'];result=cache[eta,b]
            value=physical_prefactor(1e25,2)*result['collision_moments'][0][0]
            check('retained finite displaced-distribution drag',value,control['reference'],control['reference'],2e-8,eta=eta,b_thermal=b)
        elif control['name']=='classical_drag_limit':
            eta,b=-25.,control['b_thermal'];result=cache.get((eta,b)) or born_response(eta,b)
            check('retained classical Born drag',result['collision_moments'][0][0]/math.exp(eta),control['reference'],control['reference'],2e-9,b_thermal=b)
    # Fixed-ne cold limit: the saved PBHY reference independently integrated q.
    cold_reference=next(c['reference'] for c in old['checks'] if c['name']=='degenerate_PBHY_limit')
    ne=1e27;ni=ne/2;length=1e-9;pf=HBAR*(3*math.pi**2*ne)**(1/3);cold=[]
    for theta in [1e-2,1e-3,1e-4]:
        T=theta*pf**2/(2*ME*KB);eta=chemical_potential(ne,T);b=8*ME*KB*T*length**2/HBAR**2
        result=born_response(eta,b);m=np.array(result['collision_moments'])*physical_prefactor(ni,2)
        kappa=ne**2*KB**2*T/m[1,1]
        reference=math.pi**2*ne**2*KB**2*T/(3*cold_reference)
        cold.append(dict(T_over_TF=theta,eta=eta,b_thermal=b,
            current_variance=result['current_variance'],cross_over_drag=float(m[0,1]/m[0,0]),
            heat_over_drag=float(m[1,1]/m[0,0]),kappa_relative_to_degenerate=kappa/reference-1,
            two_mode_to_exact_heat=result['two_mode_to_exact_zero_drift_heat']))
    check('degenerate current variance',cold[-1]['current_variance'],math.pi**2/3,math.pi**2/3,2e-7)
    check('degenerate heat-to-drag limit',cold[-1]['heat_over_drag'],3/math.pi**2,3/math.pi**2,2e-7)
    check('degenerate conductivity limit',cold[-1]['kappa_relative_to_degenerate'],0.,1.,2e-7)
    check('degenerate thermal coupling vanishes',abs(cold[-1]['cross_over_drag']),0.,1.,2e-4)
    check('degenerate heat truncation vanishes',cold[-1]['two_mode_to_exact_heat'],1.,1.,2e-7)
    # Tail and quadrature controls test the entire matrix, including small cross terms.
    for eta,b in [(-8.,.1),(6.,10.),(1e4,.05)]:
        fine=born_response(eta,b,tail=60.,epsrel=2e-12);base=cache[eta,b]
        for key in ['collision_moments','exact_elastic_mobility']:
            reference=np.array(fine[key]);check('tail and adaptive tolerance '+key,base[key],reference,norm(reference),2e-8,eta=eta,b_thermal=b)
    guards=0
    for fn in [lambda:current_basis(np.nan),lambda:current_basis(1,tail=10),
               lambda:born_response(1,0),lambda:born_response(1,1,basis=current_basis(2)),
               lambda:collision_response(1,lambda x:-1),lambda:physical_prefactor(0,1)]:
        try:fn()
        except ValueError:guards+=1
    check('invalid inputs',guards,6,1.,0.)
    paths=[Path(__file__),Path('scripts/electron_ion_heat.py'),Path('scripts/electron_ion_born.py'),old_path]
    passed=all(c['passed'] for c in checks)
    result=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed' if passed else 'failed',
                accepted_for_stellar_evolution=False,checks=checks,guards=guards,comparisons=comparisons,
                classical_constant_log_heat_ratio=lorentz['two_mode_to_exact_zero_drift_heat'],degenerate_convergence=cold,
                inputs_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                limitations=['Nonrelativistic fixed-ion elastic scattering','Born Yukawa potential and no correlations',
                             'Two-mode truncation explicitly measured against exact Lorentz response',
                             'No electron-electron collisions, moving-ion heat coupling or stellar diffusion'])
    with a.report.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(outcome=result['outcome'],checks=len(checks),comparisons=len(comparisons),
                         failures=[c for c in checks if not c['passed']])),flush=True)
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
