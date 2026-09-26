#!/usr/bin/env python3
"""Independent analytic-mixture and Burgers checks of composition forces.

Ideal classical ions and an analytic classical-plus-polytropic electron free
energy provide a thermodynamically consistent test EOS. Collision matrices
are synthetic. Finite differences of chemical potentials at fixed pressure
test the nonthermal force convention, independently of its analytic Hessian.
No stellar diffusion velocity or accepted thermal closure is claimed.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from diffusion_burgers import KB,MU,ME,solve as burgers
from diffusion_composition_forces import chemical_acceleration,solve_momentum,velocity_basis

R=KB/MU


class AnalyticMixture:
    def __init__(self,a,z,T,rho,x,electron_compressibility_excess):
        self.a=np.asarray(a);self.z=np.asarray(z);self.gamma=5/3
        ne=rho/MU*np.dot(x,self.z/self.a)
        self.K=electron_compressibility_excess*KB*T/(self.gamma*ne**(self.gamma-1))

    def particles(self,T,rho,x):
        n=rho/MU*x/self.a;ne=n@self.z
        C=1+self.gamma*self.K*ne**(self.gamma-1)/(KB*T)
        return n,ne,C

    def pressure(self,T,rho,x):
        n,ne,C=self.particles(T,rho,x)
        return (sum(n)+ne)*KB*T+self.K*ne**self.gamma

    def chemical_exchange(self,T,rho,x,indices,reference):
        n,ne,C=self.particles(T,rho,x)
        # Constants inside the logarithms set arbitrary chemical references.
        mu_i=KB*T*(np.log(n/1e25)-1.5*np.log(T/1e7))
        mu_e=KB*T*(np.log(ne/1e25)-1.5*np.log(T/1e7))
        mu_e+=self.gamma/(self.gamma-1)*self.K*ne**(self.gamma-1)
        neutral=(mu_i+self.z*mu_e)/(self.a*MU)
        return neutral[indices]-neutral[reference]

    def derivatives(self,T,rho,x,indices,reference):
        n,ne,C=self.particles(T,rho,x);a,z=self.a,self.z
        ye=np.dot(x,z/a);alpha=z[indices]/a[indices]-z[reference]/a[reference]
        density_gradient=R*(1/a[indices]-1/a[reference]+C*alpha)
        h=R*(np.diag(1/(a[indices]*x[indices]))+
             np.ones((len(indices),len(indices)))/(a[reference]*x[reference])+
             C/ye*np.outer(alpha,alpha))
        pr=(sum(n)+C*ne)*KB*T;pt=(sum(n)+ne)*KB*T
        return h,density_gradient,pt/pr,pr,pt

    def rho_at_pressure(self,T,x,P,rho):
        return rho*np.exp(brentq(lambda lr:self.pressure(T,rho*np.exp(lr),x)/P-1,
                                -.2,.2,xtol=1e-14,rtol=1e-14))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output',type=Path)
    a=parser.parse_args();checks=[];convention_differences=[];cases=[]
    def check(name,error,tolerance=1e-10,**metadata):
        error=float(error)
        checks.append(dict(name=name,error=error,tolerance=tolerance,passed=bool(np.isfinite(error) and error<tolerance),**metadata))
    def relative(actual,expected):
        return np.max(np.abs(np.asarray(actual)-expected))/max(np.max(np.abs(expected)),1e-100)
    for masses,charges in [([1.,4.],[1.,2.]),([1.,3.,4.],[1.,2.,2.]),
                            ([1.,3.,4.,12.,56.],[1.,2.,2.,6.,26.])]:
        masses=np.array(masses);charges=np.array(charges);ions=len(masses)
        for trace in [.15,1e-8,1e-16]:
            x=np.arange(1,ions+1,dtype=float);x[0]=0;x*=(1-trace)/sum(x);x[0]=trace
            ref=ions-1;ind=np.arange(ions-1);rho=1300.;T=1.1e7;g=2.1e5;length=1e9
            for excess in [0.,.5,5.,50.]:
                eos=AnalyticMixture(masses,charges,T,rho,x,excess)
                h,gr,delta,pr,pt=eos.derivatives(T,rho,x,ind,ref)
                n,ne,C=eos.particles(T,rho,x);pop=np.append(n,ne)
                mass=np.append(masses,ME/MU);charge=np.append(charges,-1.)
                reduced=mass[:,None]*mass[None,:]/(mass[:,None]+mass[None,:])
                k=np.outer(pop/1e25,pop/1e25)*np.sqrt(reduced)*1e9
                for thermal in [0.,-.2/length]:
                    cx=.03*x[ind]*np.sin(ind+1)/length
                    dx=np.zeros(ions);dx[ind]=cx;dx[ref]=-sum(cx)
                    lr=(-rho*g-rho*T*gr@cx-pt*thermal)/pr
                    force=chemical_acceleration(temperature=T,hessian_phi=h,gradient_phi_lnrho=gr,
                        material_delta=delta,composition_gradient=cx,logarithmic_density_gradient=lr,
                        logarithmic_temperature_gradient=thermal)
                    args=dict(density=rho,mass_fractions=x,mass_numbers=masses,charges=charges,
                              independent_ions=ind,reference_ion=ref,resistance=k)
                    solution=solve_momentum(**args,chemical_driving_acceleration=force)
                    meta=dict(ions=ions,trace_fraction=trace,electron_compressibility=C,
                              temperature_logarithmic_gradient=thermal)
                    v=solution.velocity_cm_s;speed=max(max(abs(v)),1e-100)
                    check('baryon flux',abs(x@v[:-1])/speed,**meta)
                    check('electric current',abs((pop*charge)@v)/(np.sum(pop*abs(charge))*speed),**meta)
                    check('positive dissipation balance',relative(solution.driving_power_erg_cm3_s,
                          solution.frictional_heating_erg_cm3_s),**meta)
                    check('nonnegative dissipation',max(-solution.frictional_heating_erg_cm3_s,0),**meta)
                    dp=n*KB*T*(lr+dx/x+thermal)
                    old=burgers(number_density=pop,mass_u=mass,charge=charge,ion_pressure_gradient=dp,
                         gravity=g,temperature=T,temperature_gradient=T*thermal,resistance=k,
                         z=np.zeros_like(k),zprime=np.ones_like(k),zdoubleprime=np.ones_like(k),heat_flow='suppressed')
                    # Independent projection of the explicit pressure-gradient
                    # forces in the original Burgers momentum equations.
                    electron_dp=ne*KB*T*(C*(lr+np.dot(dx,charges/masses)/np.dot(x,charges/masses))+thermal)
                    alpha=charges[ind]/masses[ind]-charges[ref]/masses[ref]
                    pressure_force=dp[ind]/(rho*x[ind])-dp[ref]/(rho*x[ref])+alpha*electron_dp/(ne*MU)
                    pressure_solution=solve_momentum(**args,chemical_driving_acceleration=pressure_force)
                    check('projected explicit Burgers forces',relative(pressure_solution.velocity_cm_s,old.velocity_cm_s),**meta)
                    if thermal==0 or excess==0:
                        check('isothermal or classical chemical-force equivalence',relative(v,old.velocity_cm_s),**meta)
                    else:
                        expected_difference=R*T*((delta-1)*(1/masses[ind]-1/masses[ref])+
                                               (delta*C-1)*alpha)*thermal
                        check('thermal-convention transformation',relative(force-pressure_force,expected_difference),1e-8,**meta)
                        convention_differences.append(dict(**meta,
                            maximum_force_difference_cm_s2=float(max(abs(force-pressure_force))),
                            relative_velocity_difference=float(relative(v,old.velocity_cm_s))))
                    # Change reference ion and rederive all exchanges; the
                    # physical solution cannot depend on that choice.
                    alt=0;others=np.arange(1,ions)
                    if trace==.15:
                        hh,gg,dd,_,_=eos.derivatives(T,rho,x,others,alt)
                        alt_force=chemical_acceleration(temperature=T,hessian_phi=hh,gradient_phi_lnrho=gg,
                            material_delta=dd,composition_gradient=dx[others],logarithmic_density_gradient=lr,
                            logarithmic_temperature_gradient=thermal)
                        other=solve_momentum(**{**args,'independent_ions':others,'reference_ion':alt},
                                              chemical_driving_acceleration=alt_force)
                        check('reference-ion invariance',relative(other.velocity_cm_s,v),**meta)
                        # Finite differences of full chemical potentials along
                        # r, subtracting a separate constant-P temperature path.
                        measured=[];P=eos.pressure(T,rho,x)
                        for step in [1e-4,5e-5]:
                            distance=length*step
                            radial=(eos.chemical_exchange(T*np.exp(thermal*distance),rho*np.exp(lr*distance),x+dx*distance,ind,ref)-
                                    eos.chemical_exchange(T*np.exp(-thermal*distance),rho*np.exp(-lr*distance),x-dx*distance,ind,ref))/(2*distance)
                            tp,tm=T*np.exp(step),T*np.exp(-step)
                            rp=eos.rho_at_pressure(tp,x,P,rho);rm=eos.rho_at_pressure(tm,x,P,rho)
                            temp=(eos.chemical_exchange(tp,rp,x,ind,ref)-eos.chemical_exchange(tm,rm,x,ind,ref))/(tp-tm)
                            measured.append(radial-temp*T*thermal)
                        check('independent chemical-potential force',relative(measured[-1],force),2e-6,**meta)
                        check('chemical-potential finite-difference spacing',relative(measured[0],measured[1]),2e-6,**meta)
                    cases.append(dict(**meta,scaled_condition_number=solution.scaled_condition_number))
    # A stationary metal inventory imposes an additional explicit constraint;
    # the two H/He exchange directions still conserve mass and charge.
    x=np.array([.1,.03,.85,.02]);mass=np.array([1.,3.,4.,56.]);z=np.array([1.,2.,2.,26.])
    b=velocity_basis(mass_fractions=x,mass_numbers=mass,charges=z,independent_ions=[0,1],reference_ion=2)
    check('stationary metal basis',np.max(abs(b[3])))
    check('stationary metal baryon constraints',np.max(abs(x@b[:-1])))
    check('stationary metal charge constraints',np.max(abs((x*z/mass)@b[:-1]-np.dot(x,z/mass)*b[-1])))
    base=dict(density=1.,mass_fractions=x,mass_numbers=mass,charges=z,
              independent_ions=[0,1],reference_ion=2,resistance=np.ones((5,5)),
              chemical_driving_acceleration=[1.,-2.])
    result=solve_momentum(**base);check('stationary metal velocity',abs(result.velocity_cm_s[3]))
    zero_e=base['resistance'].copy();zero_e[-1,:]=0;zero_e[:,-1]=0
    result=solve_momentum(**{**base,'resistance':zero_e})
    check('zero-electron-drag determined limit',result.scaled_backward_error)
    guard_count=0
    for patch in [dict(density=0),dict(mass_fractions=[0,.13,.85,.02]),
                  dict(charges=[1,2,2,26.1]),dict(independent_ions=[0,0]),
                  dict(reference_ion=0),dict(resistance=np.zeros((5,5))),
                  dict(chemical_driving_acceleration=[np.nan,1]),
                  dict(resistance=np.triu(np.ones((5,5))))]:
        try:solve_momentum(**{**base,**patch})
        except ValueError:guard_count+=1
    check('invalid input guards',0 if guard_count==8 else 1)
    paths=[Path(__file__),Path('scripts/diffusion_composition_forces.py'),Path('scripts/diffusion_burgers.py')]
    passed=all(c['passed'] for c in checks)
    report=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),
                outcome='passed' if passed else 'failed',accepted_for_stellar_evolution=False,
                cases=cases,checks=checks,guard_count=guard_count,check_count=len(checks),
                thermal_convention_differences=convention_differences,
                inputs_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                references=['https://doi.org/10.1093/mnras/stw751','https://arxiv.org/abs/1710.08424'],
                remaining_work=['Physical EOS family acceptance and actual chemical forces',
                                'Thermal transport in the matching force convention',
                                'Active-species EOS derivatives and conservative evolution coupling'])
    with a.output.open('x') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(outcome=report['outcome'],checks=len(checks),failures=[c for c in checks if not c['passed']],
                         cases=len(cases),thermal_convention_cases=len(convention_differences))))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
