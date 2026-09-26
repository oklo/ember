#!/usr/bin/env python3
"""Evaluate a fully ionized H-He collision diagnostic on a saved structure.

Uses directly integrated classical repulsive-ion moments, comparing the
source's approximate electron screening with numerical Fermi-gas pressure
derivatives. The binary diffusion coefficient is a reference mixture value,
not a multicomponent stellar velocity or settling-time prediction.
"""
import argparse
import json
import math
from pathlib import Path
import re
import subprocess

import numpy as np
from scipy.interpolate import CubicSpline

from fetch_tops_composition import digest
from ion_collision_integrals import KB,MU,HBAR,E2,screening_length,reduced_integral


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('preserve profile collision diagnostics')
    regime_path=Path('docs/results/diffusion_material_regime_3890gyr_v1.json')
    regime=json.loads(regime_path.read_text())
    table_report=Path('docs/results/yukawa_collision_table_v1.json')
    tr=json.loads(table_report.read_text());table_path=Path(tr['table'])
    if tr['status']!='pass' or digest(table_path)!=tr['table_sha256']:
        raise ValueError('numerical ion table has not passed or changed')
    table=json.loads(table_path.read_text())
    interp=CubicSpline(table['log10_strength'],np.log(table['dimensionless_integrals']),axis=0,extrapolate=False)
    header=Path('include/ember/gs98_mixture.hpp')
    metals=[tuple(map(float,row.split(','))) for row in re.findall(r'^\s*\{([^{}]+)\}, //',header.read_text(),re.M)]
    if len(metals)!=19 or abs(sum(m[3] for m in metals)-1)>1e-12:
        raise ValueError('GS98 source changed')
    paths=[Path(__file__),Path('scripts/ion_collision_integrals.py'),a.probe,
           Path('scripts/diffusion_electron_probe.cpp'),regime_path,table_report,table_path,header]
    inputs={str(p.resolve()):digest(p) for p in paths}
    inputs.update(regime['input_sha256'])
    rows=regime['records']
    queries=[]
    def query(row,rho_factor=1.):
        return [row['density_baryonic_g_cm3']*rho_factor,row['temperature_K'],
                row['X'],row['Y3'],row['Y4'],1-row['X']-row['Y3']-row['Y4']]
    queries=[query(r) for r in rows]
    derivative_indices=list(range(0,len(rows),max(1,len(rows)//20)))
    for i in derivative_indices:
        queries.extend([query(rows[i],math.exp(-1e-4)),query(rows[i],math.exp(1e-4))])
    result=subprocess.run([str(a.probe)],input=''.join(' '.join(format(v,'.17g') for v in q)+'\n' for q in queries),
                          text=True,capture_output=True,check=True,timeout=60)
    values=np.array([list(map(float,line.split())) for line in result.stdout.splitlines()])
    if values.shape!=(len(queries),4) or not np.all(np.isfinite(values)) or np.any(values[:,[0,1,2]]<=0):
        raise ValueError('invalid numerical electron response')
    electron_derivative_controls=[]
    for j,i in enumerate(derivative_indices):
        lo,hi=values[len(rows)+2*j:len(rows)+2*j+2]
        numerical=(hi[1]-lo[1])/(2e-4)
        error=abs(numerical/values[i,2]-1)
        electron_derivative_controls.append(dict(zone=rows[i]['zone'],relative_error=error))
        if error>1e-6:raise ValueError('electron pressure derivative fails independent difference')
    output=[]
    for i,row in enumerate(rows):
        rho,T=row['density_baryonic_g_cm3'],row['temperature_K'];x,y3,y4=row['X'],row['Y3'],row['Y4']
        z=1-x-y3-y4
        n=np.array([rho/MU*x,rho/MU*y3/3,rho/MU*y4/4]+[rho/MU*z*f/mass for charge,mass,atomic,f in metals])
        q=np.array([1.,2.,2.]+[m[0] for m in metals])
        ne=float(n@q)
        if abs(ne/values[i,0]-1)>1e-8:raise ValueError('electron abundance convention differs')
        source=screening_length(n,q,T)
        numerical=screening_length(n,q,T,electron_stiffness_erg=values[i,2]/values[i,0])
        # Explicitly test the screening sensitivity to metals treated as fully
        # stripped; omitting them is a limiting diagnostic, not a new mixture.
        no_metal_ions=screening_length(n[:3],q[:3],T)
        cases=[]
        for screen in (source,numerical):
            strength=2*E2/(KB*T*screen['length_cm'])
            k=np.exp(interp(math.log10(strength)))
            if not np.all(np.isfinite(k)):raise ValueError('stellar pair exceeds the numerical collision table')
            mu=MU*4/5;factor=math.sqrt(2*math.pi/mu)*(2*E2)**2/(KB*T)**1.5
            omega=factor*k
            resistance=16/3*n[0]*n[2]*mu*omega[0]
            binary_d=3*KB*T/(16*(n[0]+n[2])*mu*omega[0])
            if min(resistance,binary_d)<=0:raise ValueError('invalid H-He pair coefficient')
            cases.append(dict(screening=screen['electron_screening'],screening_length_cm=screen['length_cm'],
                              collision_strength=strength,omega11_cm3_s=float(omega[0]),
                              resistance_H1_He4_g_cm3_s=resistance,binary_H1_He4_diffusivity_cm2_s=binary_d,
                              z=float(1-.4*k[1]/k[0]),zprime=float(2.5-2*k[1]/k[0]+.4*k[2]/k[0]),
                              zdoubleprime=float(k[3]/k[0]),
                              printed_fit_resistance_relative_change=reduced_integral(strength,1,1)/k[0]-1,
                              de_broglie_over_screening=HBAR/math.sqrt(2*mu*KB*T)/screen['length_cm']))
        output.append(dict(zone=row['zone'],mass_fraction=row['mass_fraction'],convective=row['convective'],
                           temperature_K=T,density_baryonic_g_cm3=rho,X=x,
                           helium_coupling=source['ion_sphere_coupling'][2],
                           electron_stiffness_over_kT=values[i,2]/(values[i,0]*KB*T),
                           approximate_screening_vs_FD_resistance_change=cases[0]['resistance_H1_He4_g_cm3_s']/cases[1]['resistance_H1_He4_g_cm3_s']-1,
                           fully_stripped_metals_screening_comparison=no_metal_ions['length_cm']/source['length_cm']-1,
                           cases=cases))
    for p,h in inputs.items():
        if digest(Path(p))!=h:raise ValueError('profile collision source changed')
    report=dict(scope=__doc__,status='completed_diagnostic',accepted_for_stellar_evolution=False,
                hot_zones=len(output),hot_mass_fraction=regime['selected_hot_mass_fraction'],
                central_diagnostic=output[0],records=output,electron_derivative_controls=electron_derivative_controls,
                maximum_electron_derivative_error=max(r['relative_error'] for r in electron_derivative_controls),
                maximum_screening_prescription_resistance_change=max(abs(r['approximate_screening_vs_FD_resistance_change']) for r in output),
                collision_strength_range=[min(c['collision_strength'] for r in output for c in r['cases']),
                                          max(c['collision_strength'] for r in output for c in r['cases'])],
                numerical_screening_binary_diffusivity_range_cm2_s=[min(r['cases'][1]['binary_H1_He4_diffusivity_cm2_s'] for r in output),
                                                                 max(r['cases'][1]['binary_H1_He4_diffusivity_cm2_s'] for r in output)],
                maximum_quantum_length_ratio=max(c['de_broglie_over_screening'] for r in output for c in r['cases']),
                maximum_printed_fit_resistance_difference=max(abs(c['printed_fit_resistance_relative_change']) for r in output for c in r['cases']),
                limitations=['All species are assumed fully ionized; no selected stellar ionization prescription is validated here.',
                             'The binary H-He reference coefficient is not a full mixture diffusion velocity or settling time.',
                             'Ion interaction forces and thermodynamic composition gradients are not evaluated.',
                             'Electron-ion collisions, thermal diffusion and composition/energy coupling remain incomplete.'],
                input_sha256=inputs)
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ('status','hot_zones','maximum_electron_derivative_error',
                      'maximum_screening_prescription_resistance_change','collision_strength_range',
                      'numerical_screening_binary_diffusivity_range_cm2_s','maximum_quantum_length_ratio',
                      'maximum_printed_fit_resistance_difference')}),flush=True)


if __name__=='__main__':main()
