#!/usr/bin/env python3
"""Classical projected transport versus independently assembled Burgers equations.

Synthetic nonzero momentum/heat collision coefficients test coupled velocities,
heat flux, conservation, reciprocity and positive entropy. Maxwellian ideal
mixtures supply hydrostatic forces; this does not accept degenerate transport.
"""
import argparse,hashlib,json
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from audit_diffusion_composition_forces import AnalyticMixture,R
from diffusion_burgers import solve as burgers,KB,MU,ME
from diffusion_composition_forces import chemical_acceleration
from diffusion_reciprocal_transport import maxwellian_operator
from diffusion_thermal_transport import HeatFluxConvention


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('report',type=Path);args=p.parse_args()
    checks=[];cases=[]
    def check(name,actual,expected,tolerance=2e-10,**meta):
        reference_scale=float(np.max(abs(np.asarray(expected))))
        scale=reference_scale if reference_scale>0 else 1.
        error=float(np.max(abs(np.asarray(actual)-expected))/scale)
        checks.append(dict(name=name,relative_error=error,tolerance=tolerance,passed=bool(np.isfinite(error) and error<=tolerance),**meta))
    for masses,charges in [([1.,4.],[1.,2.]),([1.,3.,4.],[1.,2.,2.]),
                           ([1.,3.,4.,12.,56.],[1.,2.,2.,6.,26.])]:
        a,z=np.array(masses),np.array(charges);ions=len(a);ref=ions-1;ind=np.arange(ions-1)
        for trace in [.15,1e-8,1e-16]:
            x=np.arange(1,ions+1,dtype=float);x[0]=0;x*=(1-trace)/sum(x);x[0]=trace
            rho,T,g,length=1300.,1.1e7,2.1e5,1e9
            eos=AnalyticMixture(a,z,T,rho,x,0.);h,gr,delta,pr,pt=eos.derivatives(T,rho,x,ind,ref)
            nion=rho/MU*x/a;pop=np.append(nion,nion@z);mass=np.append(a,ME/MU);charge=np.append(z,-1.)
            reduced=mass[:,None]*mass[None,:]/(mass[:,None]+mass[None,:])
            k=np.outer(pop/1e25,pop/1e25)*np.sqrt(reduced)*1e9
            cz=np.full_like(k,.35);zp=np.full_like(k,1.2);zpp=np.full_like(k,1.4)
            common=dict(density=rho,temperature=T,mass_fractions=x,mass_numbers=a,charges=z,
                        independent_ions=ind,reference_ion=ref,resistance=k,collision_z=cz,
                        collision_zprime=zp,collision_zdoubleprime=zpp,energy_scale=R*T)
            op=maxwellian_operator(**common)
            enthalpy=2.5*R*T*(1/a[ind]-1/a[ref]+z[ind]/a[ind]-z[ref]/a[ref])
            convention=HeatFluxConvention(enthalpy,R*T)
            for lt in [0.,-.2/length]:
                dc=.03*x[ind]*np.sin(ind+1)/length;dx=np.zeros(ions);dx[ind]=dc;dx[ref]=-sum(dc)
                lr=(-rho*g-rho*T*gr@dc-pt*lt)/pr
                acceleration=chemical_acceleration(temperature=T,hessian_phi=h,gradient_phi_lnrho=gr,
                    material_delta=delta,composition_gradient=dc,logarithmic_density_gradient=lr,
                    logarithmic_temperature_gradient=lt)
                force=np.append(-acceleration/T,-R*lt)
                result=op.solve(force)
                old=burgers(number_density=pop,mass_u=mass,charge=charge,ion_pressure_gradient=nion*KB*T*(lr+dx/x+lt),
                    gravity=g,temperature=T,temperature_gradient=T*lt,resistance=k,z=cz,zprime=zp,zdoubleprime=zpp,heat_flow='classical')
                meta=dict(ions=ions,hydrogen=trace,temperature_logarithmic_gradient=lt)
                check('independent Burgers velocities',result.velocity_cm_s,old.velocity_cm_s,**meta)
                check('independent Burgers residual heat velocities',result.residual_heat_velocity_cm_s,old.residual_heat_velocity_cm_s,**meta)
                heat=float((pop*KB*T)@old.residual_heat_velocity_cm_s)
                check('independent residual heat flux',result.reduced_heat_flux_erg_cm2_s,heat,**meta)
                measured=np.append(result.independent_mass_flux_g_cm2_s,result.reduced_heat_flux_erg_cm2_s/(R*T))
                check('mass and heat response',measured,op.mobility@force,**meta)
                check('entropy balance',result.entropy_production_erg_cm3_s_K,result.collision_entropy_erg_cm3_s_K,**meta)
                check('positive entropy',max(0.,-result.entropy_production_erg_cm3_s_K),0.,**meta)
                scale=np.sqrt(np.diag(op.mobility));normalized=op.mobility/scale[:,None]/scale[None,:]
                check('Onsager reciprocal response',normalized,normalized.T,**meta)
                check('positive transport matrix',max(0.,-min(np.linalg.eigvalsh(normalized))),0.,**meta)
                speed=max(abs(result.velocity_cm_s))
                check('baryonic mass constraint',abs(x@result.velocity_cm_s[:-1])/speed,0.,1e-13,**meta)
                check('zero current constraint',abs((pop*charge)@result.velocity_cm_s)/(sum(pop*abs(charge))*speed),0.,1e-13,**meta)
                full=convention.full_fluxes(measured)
                direct_energy=float((pop*KB*T)@(old.residual_heat_velocity_cm_s+2.5*old.velocity_cm_s))
                check('independent full material energy flux',full[-1]*R*T,direct_energy,**meta)
                full_l=convention.full_matrix(op.mobility);full_force=convention.full_forces(force)
                check('full-energy force convention',full_l@full_force,full,**meta)
                zero=op.solve(np.zeros(ions))
                check('zero-force equilibrium',np.append(zero.velocity_cm_s,zero.reduced_heat_flux_erg_cm2_s),0.,**meta)
                cases.append(dict(**meta,condition=op.scaled_condition_number,backward_error=op.scaled_backward_error))
            # One deliberately stationary metal still exchanges heat.
            if ions==5 and trace==.15:
                partial=maxwellian_operator(**{**common,'independent_ions':[0,1],'reference_ion':2})
                stationary=partial.solve([1e-4,-2e-4,3e-4])
                check('stationary unselected metals',stationary.velocity_cm_s[3:5],0.)
                check('stationary-metal entropy balance',stationary.entropy_production_erg_cm3_s_K,stationary.collision_entropy_erg_cm3_s_K)
    guards=0
    for patch in [dict(density=0),dict(energy_scale=0),dict(collision_z=np.triu(cz)),
                  dict(resistance=np.zeros_like(k)),dict(collision_z=np.full_like(k,1e6)),
                  dict(mass_fractions=[0.,.1,.2,.3,.4])]:
        try:maxwellian_operator(**{**common,**patch})
        except ValueError:guards+=1
    check('invalid operator guards',guards,6,0.)
    passed=all(c['passed'] for c in checks)
    paths=[Path(__file__),Path('scripts/diffusion_reciprocal_transport.py'),Path('scripts/diffusion_thermal_transport.py'),
           Path('scripts/diffusion_burgers.py'),Path('scripts/diffusion_composition_forces.py'),Path('scripts/audit_diffusion_composition_forces.py')]
    report=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed' if passed else 'failed',
                accepted_for_stellar_evolution=False,cases=cases,checks=checks,guards=guards,
                inputs_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                limitations=['Classical Maxwellian particles only','Synthetic collision coefficients','No stellar settling or evolved composition'])
    with args.report.open('x') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(outcome=report['outcome'],cases=len(cases),checks=len(checks),failures=[c for c in checks if not c['passed']])))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
