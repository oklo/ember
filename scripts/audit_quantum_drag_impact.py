#!/usr/bin/env python3
"""Measure how direct quantum scattering changes prescribed stellar motion modes.

Combines accepted numerical H/He and metal scattering comparisons with the
unchanged ion resistances and electron population model. This is a comparison
of resistance for prescribed motions, not evolved diffusion velocities.
"""
import argparse
import json
import math
from pathlib import Path
import re

import numpy as np

from electron_ion_born import pair_resistance
from ion_collision_integrals import MU
from fetch_tops_composition import digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('preserve completed drag impact comparisons')
    paths=[Path(__file__),Path('scripts/electron_ion_born.py'),Path('scripts/ion_collision_integrals.py'),
           Path('include/ember/gs98_mixture.hpp'),Path('docs/results/diffusion_material_regime_3890gyr_v1.json'),
           Path('docs/results/stellar_electron_drag_3890gyr_v1.json'),
           Path('docs/results/stellar_quantum_drag_3890gyr_v1.json'),
           Path('docs/results/stellar_quantum_metal_drag_3890gyr_v1.json')]
    inputs={str(p.resolve()):digest(p) for p in paths}
    metals=[tuple(map(float,row.split(','))) for row in re.findall(r'^\s*\{([^{}]+)\}, //',paths[3].read_text(),re.M)]
    if len(metals)!=19:raise ValueError('metal inventory differs')
    material={r['zone']:r for r in json.loads(paths[4].read_text())['records']}
    old={r['zone']:r for r in json.loads(paths[5].read_text())['records']}
    corrections={}
    for path in paths[6:]:
        report=json.loads(path.read_text())
        if report['status']!='pass':raise ValueError('quantum thermal comparisons have not passed')
        for case in report['cases']:
            if not case['passed']:raise ValueError('failed quantum thermal case')
            if 'samples_file' in case:
                cp=Path(case['samples_file'])
                if digest(cp)!=case['samples_sha256']:raise ValueError('metal samples changed')
                inputs[str(cp.resolve())]=digest(cp)
            key=(case['zone'],case['screening'],case['charge'])
            if key in corrections:raise ValueError('duplicate charge correction')
            corrections[key]=case
    charge=np.array([1.,2.,2.]+[m[0] for m in metals])
    mass=np.array([1.,3.,4.]+[m[1] for m in metals])
    rows=[]
    for zone in sorted({k[0] for k in corrections}):
        r,prior=material[zone],old[zone]
        rho,T=r['density_baryonic_g_cm3'],r['temperature_K']
        fractions=np.array([r['X'],r['Y3'],r['Y4']]+[(1-r['X']-r['Y3']-r['Y4'])*m[3] for m in metals])
        n=rho/MU*fractions/mass;ne=float(n@charge)
        for screen in prior['modes'][0]['cases']:
            label=screen['screening'];length=screen['screening_length_cm']
            factors=np.array([corrections[(zone,label,int(q))]['quantum_over_Born_drag'] for q in charge])
            base=pair_resistance(electron_density=ne,ion_density=1.,temperature=T,charge=1,
                                 screening_length=length,eta=prior['eta_nonrelativistic'])['force_per_ion_per_drift_g_s']
            kb=base*n*charge**2;kq=kb*factors
            for i,name in enumerate(['H1','He3','He4']):
                c=charge[i]*n[i]/ne;rest=np.arange(len(n))!=i
                old_e=float(kb[i]*(1-c)**2+np.sum(kb[rest])*c*c)
                new_e=float(kq[i]*(1-c)**2+np.sum(kq[rest])*c*c)
                ion=prior['modes'][i]['ion_resistance_g_cm3_s']
                stored=next(v for v in prior['modes'][i]['cases'] if v['screening']==label)
                if abs(old_e/stored['electron_resistance_g_cm3_s']-1)>1e-11:
                    raise ValueError('unchanged Born comparison no longer matches')
                rows.append(dict(zone=zone,layer=corrections[(zone,label,1)]['layer'],isotope=name,screening=label,
                                 electron_resistance_Born_g_cm3_s=old_e,electron_resistance_quantum_g_cm3_s=new_e,
                                 ion_resistance_g_cm3_s=ion,electron_quantum_over_ion=new_e/ion,
                                 relative_change_in_total_resistance=(new_e-old_e)/(ion+old_e)))
    for path,h in inputs.items():
        if digest(Path(path))!=h:raise ValueError('drag comparison input changed')
    result=dict(scope=__doc__,status='completed_diagnostic',accepted_for_stellar_evolution=False,
                records=rows,quantum_corrections=len(corrections),
                maximum_absolute_total_resistance_change=max(abs(r['relative_change_in_total_resistance']) for r in rows),
                isotope_maximum_total_resistance_changes={name:max(abs(r['relative_change_in_total_resistance']) for r in rows if r['isotope']==name)
                                                          for name in ['H1','He3','He4']},
                limitations=['Three layers and two specified static screening prescriptions, with all ions assumed fully stripped.',
                             'Quantum scattering changes electron drag only; ion interactions and thermodynamic driving forces are unchanged.',
                             'Prescribed isotope-against-background motions are not a full diffusion solution or a general error bound.',
                             'Ionic correlations, dynamic screening, ionization, thermal diffusion and abundance/energy evolution remain unaccepted.'],
                input_sha256=inputs)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['status','quantum_corrections','maximum_absolute_total_resistance_change',
                                         'isotope_maximum_total_resistance_changes']}),flush=True)


if __name__=='__main__':main()
