#!/usr/bin/env python3
"""Compare electron and ion momentum drag on a saved hot stellar structure.

For each selected isotope, prescribe its motion relative to all other ions
moving together. Impose zero baryonic mass flux and electric current, then
compare the electron and ion contributions to dissipated power in this mode.
These are relative-motion comparisons, not solutions for stellar velocities.
Two specified electron screening lengths are sensitivities, not error bounds.
"""
import argparse
import json
import math
from pathlib import Path
import re

import numpy as np
from scipy.interpolate import CubicSpline

from electron_ion_born import KB, ME, HBAR, E2, CLIGHT, chemical_potential, pair_resistance
from ion_collision_integrals import MU
from fetch_tops_composition import digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('preserve stellar drag comparisons')
    paths=[Path(__file__),Path('scripts/electron_ion_born.py'),Path('scripts/ion_collision_integrals.py'),
           Path('docs/results/electron_ion_born_v1.json'),
           Path('docs/results/stellar_ion_collisions_3890gyr_v1.json'),
           Path('docs/results/diffusion_material_regime_3890gyr_v1.json'),
           Path('docs/results/yukawa_collision_table_v1.json'),Path('include/ember/gs98_mixture.hpp')]
    inputs={str(p.resolve()):digest(p) for p in paths}
    audit=json.loads(paths[3].read_text())
    if audit['status']!='pass':raise ValueError('electron collision controls have not passed')
    pair_profile=json.loads(paths[4].read_text());regime=json.loads(paths[5].read_text())
    tr=json.loads(paths[6].read_text());tp=Path(tr['table'])
    if tr['status']!='pass' or digest(tp)!=tr['table_sha256']:raise ValueError('ion table changed')
    inputs[str(tp.resolve())]=digest(tp)
    table=json.loads(tp.read_text())
    interpolation=CubicSpline(table['log10_strength'],np.log(table['dimensionless_integrals']),axis=0,extrapolate=False)
    metals=[tuple(map(float,s.split(','))) for s in re.findall(r'^\s*\{([^{}]+)\}, //',paths[7].read_text(),re.M)]
    if len(metals)!=19:raise ValueError('metal inventory changed')
    mass=np.array([1.,3.,4.]+[r[1] for r in metals])
    charge=np.array([1.,2.,2.]+[r[0] for r in metals])
    rows=[]
    for r,old in zip(regime['records'],pair_profile['records'],strict=True):
        if r['zone']!=old['zone']:raise ValueError('profiles are not the same zones')
        T,rho=r['temperature_K'],r['density_baryonic_g_cm3']
        fractions=np.array([r['X'],r['Y3'],r['Y4']]+[(1-r['X']-r['Y3']-r['Y4'])*m[3] for m in metals])
        n=rho/MU*fractions/mass;ne=float(n@charge)
        if np.any(n<=0) or abs(ne/r['electron_density_cm3']-1)>1e-12:raise ValueError('population convention differs')
        numerical=next(c for c in old['cases'] if c['screening']=='supplied_pressure_derivative')
        ion_length=numerical['screening_length_cm']
        electron_length=math.sqrt(old['electron_stiffness_over_kT']*KB*T/(4*math.pi*E2*ne))
        eta=chemical_potential(ne,T)
        electron_cases=[]
        for label,length in [('total_effective_static_screening',ion_length),
                             ('electron_only_static_screening',electron_length)]:
            # One integral per electron state and screening length. All charges
            # share this Born moment, with the exact Z^2 scaling of this model.
            e=pair_resistance(electron_density=ne,ion_density=1.,temperature=T,
                             charge=1.,screening_length=length,eta=eta)
            e['label']=label;e['screening_length_cm']=length
            electron_cases.append(e)
        modes=[]
        for i,name in enumerate(['H1','He3','He4']):
            selected=np.arange(len(n))!=i
            mu=MU*mass[i]*mass[selected]/(mass[i]+mass[selected])
            strength=charge[i]*charge[selected]*E2/(KB*T*ion_length)
            moments=np.exp(interpolation(np.log10(strength)))[:,0]
            if not np.all(np.isfinite(moments)):raise ValueError('ion collision table does not cover pair')
            omega=np.sqrt(2*math.pi/mu)*(charge[i]*charge[selected]*E2)**2/(KB*T)**1.5*moments
            k=16/3*n[i]*n[selected]*mu*omega
            ion_total=float(np.sum(k))
            c=charge[i]*n[i]/ne
            # Unit relative speed between the selected isotope and other ions.
            ion_velocity=np.full(len(n),-fractions[i]);ion_velocity[i]+=1
            electron_velocity=float((n*charge)@ion_velocity/ne)
            mass_residual=float(abs(fractions@ion_velocity))
            current_residual=float(abs((n*charge)@ion_velocity-ne*electron_velocity)/ne)
            if max(mass_residual,current_residual)>1e-12:raise ValueError('prescribed mode violates constraints')
            cases=[]
            for e in electron_cases:
                kei=e['force_per_ion_per_drift_g_s']*n*charge**2
                electron_total=float(kei[i]*(1-c)**2+np.sum(kei[selected])*c*c)
                direct=float(np.sum(kei*(ion_velocity-electron_velocity)**2))
                if abs(direct/electron_total-1)>1e-11:raise ValueError('electron dissipation reduction differs')
                cases.append(dict(screening=e['label'],screening_length_cm=e['screening_length_cm'],
                                  electron_resistance_g_cm3_s=electron_total,
                                  electron_over_ion_resistance=electron_total/ion_total,
                                  fraction_of_total_resistance=electron_total/(ion_total+electron_total)))
            modes.append(dict(isotope=name,ion_resistance_g_cm3_s=ion_total,
                              mass_constraint_residual=mass_residual,current_constraint_residual=current_residual,
                              cases=cases))
        rows.append(dict(zone=r['zone'],mass_fraction=r['mass_fraction'],convective=r['convective'],
                         temperature_K=T,density_g_cm3=rho,eta_nonrelativistic=eta,
                         fermi_momentum_over_mec=electron_cases[0]['fermi_momentum_over_mec'],
                         hydrogen_Born_parameter=electron_cases[0]['characteristic_Born_parameter'],
                         helium_Born_parameter=2*electron_cases[0]['characteristic_Born_parameter'],
                         largest_metal_Born_parameter=float(max(charge)*electron_cases[0]['characteristic_Born_parameter']),
                         electron_only_over_total_screening_length=electron_length/ion_length,modes=modes))
    summary={}
    for i,name in enumerate(['H1','He3','He4']):
        summary[name]={}
        for j,label in enumerate(['total_effective_static_screening','electron_only_static_screening']):
            all_ratios=[r['modes'][i]['cases'][j]['electron_over_ion_resistance'] for r in rows]
            radiative_ratios=[r['modes'][i]['cases'][j]['electron_over_ion_resistance'] for r in rows if not r['convective']]
            summary[name][label]=dict(all_hot_range=[min(all_ratios),max(all_ratios)],
                                     nonconvective_range=[min(radiative_ratios),max(radiative_ratios)])
    # Nested profile inputs include the actual numerical-electron EOS probe.
    inputs.update(pair_profile['input_sha256'])
    for path,h in inputs.items():
        if digest(Path(path))!=h:raise ValueError('stellar collision input changed')
    result=dict(scope=__doc__,status='completed_diagnostic',accepted_for_stellar_evolution=False,
                hot_zones=len(rows),hot_mass_fraction=regime['selected_hot_mass_fraction'],summary=summary,
                central_diagnostic=rows[0],records=rows,
                maximum_fermi_momentum_over_mec=max(r['fermi_momentum_over_mec'] for r in rows),
                helium_Born_parameter_range=[min(r['helium_Born_parameter'] for r in rows),max(r['helium_Born_parameter'] for r in rows)],
                maximum_metal_Born_parameter=max(r['largest_metal_Born_parameter'] for r in rows),
                limitations=['Nonrelativistic first Born scattering of displaced Fermi electrons from stationary ions.',
                             'Both screening choices omit an independently validated ion structure factor and dynamic screening.',
                             'Fully stripped ion charges, especially metals, are assumed; Born accuracy is not established by its parameter alone.',
                             'The prescribed relative-motion modes are not stellar diffusion velocities, timescales or rigorous error bounds.',
                             'Electron heat perturbations, ion chemical forces and abundance/energy coupling remain absent.'],
                input_sha256=inputs)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=result['status'],hot_zones=len(rows),summary=summary,
                         helium_Born_parameter_range=result['helium_Born_parameter_range'])),flush=True)


if __name__=='__main__':main()
