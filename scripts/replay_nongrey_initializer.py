#!/usr/bin/env python3
"""Replay a stopped starting model with the canonical atmosphere equations.

The stopped model is only a guess. Its original files and stop receipt are
retained, and the new canonical solve must pass all ordinary source checks.
This avoids requiring convergence of an auxiliary iteration scheme before
the physical atmosphere solver can even begin.
"""
import argparse
import itertools
import json
from pathlib import Path

from generate_nongrey_grid import (archive, atmosphere_inputs, composition,
    execute, input_fingerprint, resample_initial_structure, sequence, temperatures)
from import_nongrey_grid import source_inputs, source_state
from nongrey_opacity import validate_table
from prepare_nongrey_sources import digest


def replay(work, plane, temperature, gravity, initializer):
    work=Path(work).resolve(); initializer=Path(initializer).resolve()
    if not initializer.is_relative_to(work):
        raise ValueError('initializer must be preserved inside the grid workspace')
    stopped=json.loads((initializer/'stopped.json').read_text())
    if stopped.get('status')!='stopped initializer; not an accepted atmosphere':
        raise ValueError('an explicit initializer stop receipt is required')
    for name, expected in stopped['files_sha256'].items():
        if digest(initializer/name)!=expected:
            raise ValueError('stopped initializer changed')
    if input_fingerprint(stopped['executable_sha256'],initializer)!=stopped['input_sha256']:
        raise ValueError('stopped initializer input fingerprint mismatch')
    spec=json.loads((work/'specification.json').read_text())
    prepared=json.loads((work/'provenance.json').read_text())
    if digest(prepared['tlusty'])!=prepared['executables']['tlusty']:
        raise ValueError('canonical executable changed')
    x,y=list(itertools.product(spec['hydrogen'],spec['helium3']))[plane]
    teff=spec['teff_K'][temperature]; logg=spec['log_g'][gravity]
    if list(map(float,(initializer/'fort.5').read_text().splitlines()[0].split()))!=[teff,logg]:
        raise ValueError('initializer coordinates differ')
    abundance,masses=composition(x,y,spec['metals'])
    table=work/f'plane-{plane:03d}/opacity/fort.63'
    opacity={'temperature_K':temperatures(spec),'density_g_cm3':sequence(spec['log_density'])}
    validate_table(table,abundance,opacity['temperature_K'],opacity['density_g_cm3'])
    if digest(table)!=(initializer/'opacity.sha256').read_text().strip():
        raise ValueError('initializer opacity differs')
    initial=resample_initial_structure((initializer/'fort.7').read_text(),spec['depths'])
    target=work/f'plane-{plane:03d}/model-{temperature:03d}-{gravity:03d}'
    if (target/'validated.json').exists():
        raise FileExistsError('accepted model already exists')
    atmosphere_inputs(target,prepared,spec,table,abundance,masses,teff,logg,initial)
    execute(prepared['tlusty'],target,['fort.7','fort.9'])
    names={'log':'run.log','convergence':'fort.9','atmosphere_input':'fort.5',
           'element_masses':'ember-masses.dat','parameters':'tas','initial_structure':'fort.8'}
    log=(target/'run.log').read_text()
    source_inputs({k:(target/v).read_text() for k,v in names.items()
                   if k not in ['log','convergence']},spec,x,y,teff,logg,log)
    state=source_state(log,(target/'fort.9').read_text(),teff,logg,opacity,spec['tau'])
    if state['depths']!=spec['depths']:
        raise ValueError('canonical depth count mismatch')
    files={}
    for name in [*stopped['files_sha256'],'stopped.json']:
        p=archive(initializer/name,initializer)
        files[str(p.relative_to(work))]=digest(p)
    record={'XH':x,'X3':y,'teff_K':teff,'log_g':logg,'diagnostics':state,
            'initialization':{'method':'stopped native CONREF iterate used only as a starting guess; independent canonical solve',
                              'executable_sha256':stopped['executable_sha256'],
                              'files_sha256':files}}
    for kind,name in names.items():
        p=archive(target/name,target)
        record[kind]=str(p.relative_to(work)); record[kind+'_sha256']=digest(p)
    (target/'validated.json').write_text(json.dumps(record,indent=2)+'\n')
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('work',type=Path);p.add_argument('plane',type=int)
    p.add_argument('temperature',type=int);p.add_argument('gravity',type=int)
    p.add_argument('initializer',type=Path)
    a=p.parse_args()
    print(json.dumps(replay(a.work,a.plane,a.temperature,a.gravity,a.initializer),indent=2))
