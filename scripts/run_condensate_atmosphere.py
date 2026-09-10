#!/usr/bin/env python3
"""Solve and validate one gas-control or settled-grain atmosphere experiment."""
import argparse
import json
import re
from pathlib import Path
from generate_nongrey_grid import (atmosphere_inputs,composition,execute,
    completed_initial_structure,checkpoint_initial_structure,resample_initial_structure,
    temperatures,sequence,continuation_structure)
from generate_condensate_opacity import depletion_inputs
from import_nongrey_grid import source_inputs,source_state,read_text
from nongrey_opacity import validate_table
from prepare_nongrey_sources import digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('prepared',type=Path);p.add_argument('specification',type=Path)
    p.add_argument('table',type=Path);p.add_argument('initial',type=Path);p.add_argument('work',type=Path)
    p.add_argument('--hydrogen',type=float,required=True);p.add_argument('--helium3',type=float,default=0)
    p.add_argument('--teff',type=float,required=True);p.add_argument('--logg',type=float,required=True)
    p.add_argument('--mode',choices=['gas','equilibrium'],required=True)
    p.add_argument('--step-limit',type=float,help='Newton temperature-step limit; final tolerances unchanged')
    p.add_argument('--checkpoint-initial',action='store_true',
                   help='use an attested unfinished profile only as an initial guess; final source checks remain mandatory')
    p.add_argument('--convective-iterations',type=int,
                   help='native CONREF initializer only, 0..200; zero disables CONREF while retaining the indexed chemistry initializer')
    a=p.parse_args();prepared=json.loads(a.prepared.read_text());spec=json.loads(a.specification.read_text())
    if a.step_limit is not None:
        if not 1<a.step_limit<=1.25:raise ValueError('invalid temperature-step limit')
        spec['temperature_step_limit']=a.step_limit
    if a.convective_iterations is not None:
        if not 0<=a.convective_iterations<=200 or prepared.get('initialization_only',{}).get('method')!='native CONREF with explicit bounds repairs':
            raise ValueError('convective initialization requires the verified CONREF source')
    if digest(prepared['tlusty'])!=prepared['executables']['tlusty']:raise ValueError('source executable changed')
    basics=(Path(prepared['tlusty']).parent/'BASICS.FOR').read_text()
    capacity=re.search(r'\bMTABT\s*=\s*(\d+)',basics)
    if not capacity or len(temperatures(spec))>int(capacity[1]):
        raise ValueError('opacity temperature grid exceeds compiled TLUSTY capacity; build an isolated capacity variant')
    abundance,masses=composition(a.hydrogen,a.helium3,spec['metals'])
    validate_table(a.table,abundance,temperatures(spec),sequence(spec['log_density']))
    if a.mode=='equilibrium':
        op=json.loads((a.table.parent/'provenance.json').read_text())
        if any(op[k]!=v for k,v in [('XH',a.hydrogen),('X3',a.helium3),('mode',a.mode)]):
            raise ValueError('depletion opacity physical inputs disagree')
        for key in ['metals','log_temperature','temperature_K','log_density','opacity_frequencies','opacity_method',
                    'wavelength_A','line_threshold','synthesis_spacing_A','microturbulence_km_s']:
            if op['specification'].get(key)!=spec.get(key):raise ValueError('opacity material or spectral inputs differ')
        if op['prepared']['depletion']!=prepared['depletion']:raise ValueError('opacity and atmosphere depletion physics differ')
    initial_deck=a.initial/'fort.5'
    if not initial_deck.exists():initial_deck=initial_deck.with_name('fort.5.gz')
    source_teff,source_logg=map(float,read_text(initial_deck).splitlines()[0].split())
    checkpoint=None
    if a.checkpoint_initial:
        checkpoint=json.loads((a.initial/'checkpoint.json').read_text())
        required={'provenance.json','ember-condensates.cfg','condensate-abundances.dat','condensates.sha256'}
        if not required.issubset(checkpoint['files']):
            raise ValueError('depleted checkpoint lacks attested physical inputs')
        initial=checkpoint_initial_structure(a.initial,a.hydrogen,a.helium3,source_teff,source_logg)
        original=json.loads((a.initial/'provenance.json').read_text())
        if original['mode']!=a.mode or original['prepared']['depletion']!=prepared['depletion']:
            raise ValueError('checkpoint depletion physics differs')
        initial=resample_initial_structure(initial,spec['depths'])
    else:
        initial=completed_initial_structure(a.initial,spec,source_teff,source_logg,spec['depths'])
    if initial is None:raise ValueError('initial model resolution exceeds target')
    initial=continuation_structure(initial,source_teff,source_logg,a.teff,a.logg)
    a.work.mkdir(parents=True,exist_ok=True)
    provenance={'prepared':prepared,'specification':spec,'mode':a.mode,'XH':a.hydrogen,'X3':a.helium3,
                'teff_K':a.teff,'log_g':a.logg,'opacity_sha256':digest(a.table),'initial':str(a.initial.resolve()),
                'initial_coordinates':{'teff_K':source_teff,'log_g':source_logg},
                'initialization':'pressure-preserving temperature/gravity continuation; final source solved anew'}
    if a.convective_iterations is not None:provenance['convective_iterations']=a.convective_iterations
    if checkpoint is not None:
        provenance['initialization_checkpoint']={'kind':'unconverged_initial_guess',
            'receipt_sha256':digest(a.initial/'checkpoint.json'),
            'source_executable_sha256':checkpoint['source_executable_sha256'],
            'scope':'Attested unfinished state used only as a guess; no convergence or flux acceptance inherited.'}
    initial_provenance=a.initial/'provenance.json'
    if initial_provenance.exists():
        previous=json.loads(initial_provenance.read_text())
        if checkpoint is None and 'initialization_only' in previous['prepared']:
            provenance['initialization_provenance']={'prepared':previous['prepared'],
                'source_receipt_sha256':digest(a.initial/'completed.json'),
                'provenance_sha256':digest(initial_provenance)}
    receipt=a.work/'provenance.json'
    if receipt.exists() and json.loads(receipt.read_text())!=provenance:raise ValueError('work inputs changed')
    receipt.write_text(json.dumps(provenance,indent=2)+'\n')
    atmosphere_inputs(a.work,prepared,spec,a.table,abundance,masses,a.teff,a.logg,initial)
    if a.convective_iterations is not None:
        tas=a.work/'tas'
        tas.write_text(tas.read_text()+f'ICONRE={a.convective_iterations},ICONRS=1,IMUCON=200,IDEEPC=3,CRFLIM=-1\n')
    depletion_inputs(a.work,prepared,spec,a.hydrogen,a.helium3,a.mode)
    execute(prepared['tlusty'],a.work,['fort.7','fort.9'])
    log=(a.work/'run.log').read_text()
    source_inputs({k:(a.work/v).read_text() for k,v in [('atmosphere_input','fort.5'),('parameters','tas'),
        ('element_masses','ember-masses.dat'),('initial_structure','fort.8')]},spec,a.hydrogen,a.helium3,a.teff,a.logg,log)
    state=source_state(log,(a.work/'fort.9').read_text(),a.teff,a.logg,
        {'temperature_K':temperatures(spec),'density_g_cm3':sequence(spec['log_density'])},spec['tau'])
    # A modified numerical initializer supplies a guess, never an accepted
    # atmosphere. Its original-source replay receives the final receipt.
    filename='initializer.json' if 'initialization_only' in prepared else 'validated.json'
    (a.work/filename).write_text(json.dumps({'provenance':provenance,'diagnostics':state},indent=2)+'\n')
    print(json.dumps(state,indent=2))


if __name__=='__main__':main()
