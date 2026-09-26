#!/usr/bin/env python3
"""Identify the kinetic regime for microscopic diffusion on a saved structure.

These fully ionized diagnostics do not calculate diffusion coefficients,
velocities, abundance changes or settling times. They identify which plasma
approximations a future transport calculation must test. No stellar input or
evolution output is changed.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import re

from fetch_tops_composition import digest


MU = 1.66053906660e-24
ME = 9.1093837015e-28
KB = 1.380649e-16
HBAR = 1.054571817e-27
C = 2.99792458e10
E2 = 2.3070775523417355e-19


def diagnostics(rho, t, x, y3, y4, moments):
    z = 1-x-y3-y4
    if min(x, y3, y4, z) < 0 or not rho > 0 or not t > 0:
        raise ValueError('invalid baryonic material state')
    ni = rho/MU * (x+y3/3+y4/4+z*moments[0])
    ne = rho/MU * (x+2*y3/3+y4/2+z*moments[1])
    ae = (3/(4*math.pi*ne))**(1/3)
    q = HBAR*(3*math.pi**2*ne)**(1/3)/(ME*C)
    ef = ME*C*C*q*q/(math.sqrt(1+q*q)+1)
    nhe = rho*y4/(4*MU)
    tp = HBAR*math.sqrt(4*math.pi*4*E2*nhe/(4*MU))/KB
    return {'ion_density_cm3': ni, 'electron_density_cm3': ne,
            'fermi_momentum_over_mc': q, 'fermi_temperature_K': ef/KB,
            'temperature_over_fermi_temperature': t*KB/ef,
            'electron_sphere_radius_cm': ae,
            'helium_coulomb_coupling': 2**(5/3)*E2/(ae*KB*t),
            'helium4_plasma_temperature_K': tp,
            'temperature_over_helium4_plasma_temperature': t/tp}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--profile', type=Path, required=True)
    p.add_argument('--evolution-audit', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('preserve earlier regime reports')
    header = Path('include/ember/gs98_mixture.hpp')
    metals = [tuple(map(float, row.split(','))) for row in re.findall(r'^\s*\{([^{}]+)\}, //', header.read_text(), re.M)]
    if len(metals) != 19 or abs(sum(r[3] for r in metals)-1) > 1e-12:
        raise ValueError('unexpected GS98 source inventory')
    moments = [sum(f*z**power/m for z,m,atomic,f in metals) for power in (0,1)]
    with a.profile.open() as f:
        profile = [{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]
    audit = json.loads(a.evolution_audit.read_text())
    if len(profile) != 512:
        raise ValueError('unexpected saved structure size')
    mass = [row['mass_g'] for row in profile]
    if any(v <= 0 for v in mass) or any(v >= w for v,w in zip(mass,mass[1:])):
        raise ValueError('invalid enclosed mass grid')
    edges = [0.] + [(v+w)/2 for v,w in zip(mass,mass[1:])] + [mass[-1]]
    weights = [(v-u)/mass[-1] for u,v in zip(edges,edges[1:])]
    convective = [None]*len(profile)
    for region in audit['mixing_regions']:
        begin,end = region['begin'],region['end_exclusive']
        if abs(sum(weights[begin:end])-region['mass_fraction']) > 1e-12:
            raise ValueError('saved mixing region does not match the current profile mass grid')
        if any(v is not None for v in convective[begin:end]):
            raise ValueError('overlapping mixing regions')
        convective[begin:end] = [region['convective']]*(end-begin)
    if any(v is None for v in convective):
        raise ValueError('incomplete mixing classification')
    rows = []
    omitted = []
    for i,row in enumerate(profile):
        if row['temperature_K'] < 2e6:
            omitted.append({'zone':i, 'reason':'outside the selected hot, fully ionized diagnostic scope'})
            continue
        d = diagnostics(row['density_g_cm3'],row['temperature_K'],row['X'],row['Y3'],row['Y4'],moments)
        rows.append({'zone':i,'mass_fraction':weights[i], 'enclosed_mass_fraction':mass[i]/mass[-1],
                     'convective':convective[i], 'X':row['X'], 'Y3':row['Y3'], 'Y4':row['Y4'],
                     'temperature_K':row['temperature_K'], 'density_baryonic_g_cm3':row['density_g_cm3'], **d})
    # Independent SI evaluation catches electrostatic, density and mass-unit
    # factors. This matters especially for the ion plasma temperature.
    e_si = 1.602176634e-19
    eps0 = 8.8541878128e-12
    controls = []
    for row in rows[::max(1,len(rows)//12)]:
        ne_si = row['electron_density_cm3']*1e6
        rho_si = row['density_baryonic_g_cm3']*1000
        mu_si = MU*.001
        p_si = HBAR*1e-7*(3*math.pi**2*ne_si)**(1/3)
        me_si,c_si,k_si = ME*.001,C*.01,KB*1e-7
        ef_si = p_si*p_si*c_si*c_si/(math.hypot(me_si*c_si*c_si,p_si*c_si)+me_si*c_si*c_si)
        nhe_si = rho_si*row['Y4']/(4*mu_si)
        tp_si = HBAR*1e-7*math.sqrt(nhe_si*(2*e_si)**2/(eps0*4*mu_si))/k_si
        ae_si = (3/(4*math.pi*ne_si))**(1/3)
        gamma_si = 2**(5/3)*e_si*e_si/(4*math.pi*eps0*ae_si*k_si*row['temperature_K'])
        controls.append({'zone':row['zone'],'fermi_temperature_relative_error':ef_si/k_si/row['fermi_temperature_K']-1,
                         'plasma_temperature_relative_error':tp_si/row['helium4_plasma_temperature_K']-1,
                         'coupling_relative_error':gamma_si/row['helium_coulomb_coupling']-1})
    max_unit_error = max(abs(v) for row in controls for k,v in row.items() if k!='zone')
    if max_unit_error > 1e-8:
        raise ValueError('SI/cgs diagnostic mismatch')
    pure = diagnostics(1e4,1e7,0.,0.,1.,moments)
    ai = (3/(4*math.pi*(1e4/(4*MU))))**(1/3)
    gamma_direct = 4*E2/(ai*KB*1e7)
    pure_error = abs(gamma_direct/pure['helium_coulomb_coupling']-1)
    if pure_error > 1e-12:
        raise ValueError('pure helium ion-sphere limit differs')
    selected = [r for r in rows if not r['convective']]
    result = {'scope':__doc__, 'accepted_diffusion_prescription':False,
        'profile':str(a.profile), 'evolution_audit':str(a.evolution_audit),
        'fully_ionized_diagnostic_minimum_temperature_K':2e6,
        'selected_hot_zones':len(rows), 'selected_hot_mass_fraction':sum(r['mass_fraction'] for r in rows),
        'selected_nonconvective_mass_fraction':sum(r['mass_fraction'] for r in selected),
        'central_diagnostics':rows[0],
        'nonconvective_mass_fraction_with_T_over_TF_below_0p3':sum(r['mass_fraction'] for r in selected if r['temperature_over_fermi_temperature']<.3),
        'nonconvective_mass_fraction_with_helium_Gamma_above_0p1':sum(r['mass_fraction'] for r in selected if r['helium_coulomb_coupling']>.1),
        'nonconvective_mass_fraction_with_helium_Gamma_above_1':sum(r['mass_fraction'] for r in selected if r['helium_coulomb_coupling']>1),
        'maximum_hot_helium_Gamma':max(r['helium_coulomb_coupling'] for r in rows),
        'minimum_hot_temperature_over_helium_plasma_temperature':min(r['temperature_over_helium4_plasma_temperature'] for r in rows),
        'unit_controls':controls, 'maximum_SI_cgs_relative_difference':max_unit_error,
        'pure_helium_ion_sphere_relative_difference':pure_error, 'records':rows, 'omissions':omitted,
        'definitions':{'TF':'Zero-temperature relativistic kinetic Fermi energy divided by kB; rest mass excluded.',
            'Gamma_He':'Z^(5/3) e^2/(a_e kB T), Z=2, a_e=(3/(4 pi n_e))^(1/3).',
            'Tp_He4':'hbar sqrt(4 pi Z^2 e^2 n_He4/(4 m_u))/kB; helium-component reference, not a mixed-plasma mode spectrum.'},
        'limitations':['Fully ionized particle counts diagnose the hot material; they do not validate an ionization treatment for diffusion.',
            'These are local dimensionless regimes, not velocities, settling times, or a prediction of the surface composition.',
            'Small T/TF requires electron degeneracy in the electric-field and chemical-potential balance.',
            'The report does not establish that diffusion can be omitted until a particular age.'],
        'references':{'burgers_transport':'https://arxiv.org/abs/astro-ph/9304005',
                      'degenerate_electron_diffusion':'https://arxiv.org/abs/1710.08424',
                      'ion_collision_transport':'https://doi.org/10.1103/PhysRevE.93.043203'},
        'input_sha256':{str(p.resolve()):digest(p) for p in [a.profile,a.evolution_audit,header,Path(__file__)]}}
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('scope','records','omissions','input_sha256','definitions','limitations','references','unit_controls')}))


if __name__ == '__main__':
    main()
