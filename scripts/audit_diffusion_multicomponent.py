#!/usr/bin/env python3
"""Check the frozen transport solver with several ion species and heat flow.

All collision coefficients are synthetic numerical controls. This report
accepts neither a collision prescription nor stellar diffusion velocities.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from audit_diffusion_burgers import independent_balances
from diffusion_burgers import KB, ME, MU, solve
from fetch_tops_composition import digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    if args.output.exists():
        raise FileExistsError('preserve completed multicomponent checks')
    paths=[Path(__file__),Path('scripts/diffusion_burgers.py'),
           Path('scripts/audit_diffusion_burgers.py'),Path('docs/results/diffusion_burgers_equations_v1.json')]
    inputs={str(p.resolve()):digest(p) for p in paths}
    controls=[]
    masses=np.array([1.,3.,4.,12.,16.,20.,24.,28.,32.,40.,56.])
    charges=np.array([1.,2.,2.,6.,8.,10.,12.,14.,16.,20.,26.])
    def check(name,error,**details):
        error=float(error)
        if not np.isfinite(error) or error>1e-10:
            raise ValueError(f'{name}: {error}')
        controls.append(dict(name=name,error=error,tolerance=1e-10,**details))
    def compare(name,actual,expected,**details):
        expected=np.asarray(expected)
        check(name,np.max(np.abs(np.asarray(actual)-expected))/max(np.max(np.abs(expected)),1e-100),**details)
    runs=[]
    for ions in (3,5,11):
        for trace in (1e-4,1e-16):
            abundance=np.logspace(-3,0,ions);abundance[0]=trace
            n=np.append(abundance,np.dot(abundance,charges[:ions]))*1e25
            a=np.append(masses[:ions],ME/MU);q=np.append(charges[:ions],-1.)
            reduced=a[:,None]*a[None,:]/(a[:,None]+a[None,:])
            k=np.outer(n/1e25,n/1e25)*np.sqrt(reduced)*1e9
            zz=.45+.1/(1+np.abs(a[:,None]-a[None,:]))
            case=dict(number_density=n,mass_u=a,charge=q,
                      ion_pressure_gradient=n[:-1]*KB*(-1e-3*(1+.05*np.arange(ions))-.02),
                      temperature=1e7,temperature_gradient=-.02,gravity=2e5,
                      resistance=k,z=zz,zprime=np.full(k.shape,1.3),
                      zdoubleprime=np.full(k.shape,2.))
            for mode in ('classical','suppressed'):
                case['heat_flow']=mode
                base=solve(**case)
                details=dict(ion_species=ions,trace_number_ratio=trace,heat_flow=mode)
                balance=independent_balances(case,base)
                for name,value in balance.items():
                    check(name,value,**details)
                runs.append({**details,**balance,'scaled_condition_number':base.scaled_condition_number})
                permutation=np.append(np.arange(ions)[::-1],ions)
                other={**case,**{key:case[key][permutation] for key in ('number_density','mass_u','charge')},
                       'ion_pressure_gradient':case['ion_pressure_gradient'][::-1],
                       **{key:case[key][np.ix_(permutation,permutation)] for key in
                          ('resistance','z','zprime','zdoubleprime')}}
                result=solve(**other)
                compare('species permutation velocity',result.velocity_cm_s,base.velocity_cm_s[permutation],**details)
                if mode=='classical':
                    compare('species permutation heat flow',result.residual_heat_velocity_cm_s,
                            base.residual_heat_velocity_cm_s[permutation],**details)
                compare('species permutation electric field',result.electric_force_dyn,base.electric_force_dyn,**details)
                pieces=[]
                for keep in ('gravity','ion_pressure_gradient','temperature_gradient'):
                    drive={**case,'gravity':0.,'ion_pressure_gradient':np.zeros(ions),'temperature_gradient':0.}
                    drive[keep]=case[keep]
                    pieces.append(solve(**drive))
                compare('separate physical-force contributions',sum(r.velocity_cm_s for r in pieces),base.velocity_cm_s,**details)
                if mode=='classical':
                    compare('separate heat-flow contributions',sum(r.residual_heat_velocity_cm_s for r in pieces),
                            base.residual_heat_velocity_cm_s,**details)
                compare('separate electric-field contributions',sum(r.electric_force_dyn for r in pieces),base.electric_force_dyn,**details)
    for p,h in inputs.items():
        if digest(Path(p))!=h:
            raise ValueError('transport source changed')
    report=dict(scope=__doc__,status='pass',accepted_for_stellar_evolution=False,
                checks=controls,check_count=len(controls),cases=runs,
                maximum_pair_equation_residual=max(r['pair_equation_relative_residual'] for r in runs),
                maximum_baryon_flux_residual=max(r['baryon_flux_relative_residual'] for r in runs),
                maximum_current_residual=max(r['current_relative_residual'] for r in runs),
                maximum_scaled_condition_number=max(r['scaled_condition_number'] for r in runs),
                input_sha256=inputs)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('scope','checks','cases','input_sha256')}))


if __name__=='__main__':
    main()
