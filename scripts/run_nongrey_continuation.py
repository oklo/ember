#!/usr/bin/env python3
"""Continue a verified gas atmosphere, then replay the canonical source.

Native convective refinement is an optional initializer only. Both source
runs retain their input decks, complete outputs and independent checks.
"""
import argparse
import json
from pathlib import Path
from generate_nongrey_grid import (atmosphere_inputs,composition,execute,temperatures,
    sequence,completed_initial_structure,continuation_structure,archive)
from import_nongrey_grid import source_inputs,source_state
from nongrey_opacity import validate_table
from prepare_nongrey_sources import digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['prepared','specification','opacity','initial','work']:
        p.add_argument(name,type=Path)
    p.add_argument('--hydrogen',type=float,required=True)
    p.add_argument('--helium3',type=float,default=0)
    p.add_argument('--teff',type=float,required=True)
    p.add_argument('--logg',type=float,required=True)
    p.add_argument('--initializer',type=Path,required=True)
    p.add_argument('--convective-iterations',type=int,default=200)
    p.add_argument('--composition-continuation',action='store_true',
                   help='allow a verified starting guess at another H/He composition; target chemistry and opacity are solved anew')
    a=p.parse_args();prepared=json.loads(a.prepared.read_text());spec=json.loads(a.specification.read_text())
    init=json.loads(a.initializer.read_text())
    if 'depletion' in prepared or 'initialization_only' in prepared:
        raise ValueError('canonical gas source required')
    if digest(prepared['tlusty'])!=prepared['executables']['tlusty']:
        raise ValueError('canonical executable changed')
    if digest(init['executable'])!=init['executable_sha256'] or digest(
            Path(prepared['tlusty']).parent/'tlusty208.f')!=init['base_source_sha256']:
        raise ValueError('convective initializer identity mismatch')
    if not 1<=a.convective_iterations<=200:raise ValueError('invalid refinement count')
    initial_record=json.loads((a.initial/'validated.json').read_text())
    changed_composition=initial_record['XH']!=a.hydrogen or initial_record['X3']!=a.helium3
    if changed_composition and not a.composition_continuation:
        raise ValueError('continuation requires the same composition')
    original=completed_initial_structure(a.initial,spec,initial_record['teff_K'],
                                         initial_record['log_g'],spec['depths'])
    if original is None:raise ValueError('initial model exceeds requested depth count')
    text=continuation_structure(original,initial_record['teff_K'],initial_record['log_g'],a.teff,a.logg)
    abundance,masses=composition(a.hydrogen,a.helium3,spec['metals'])
    validate_table(a.opacity,abundance,temperatures(spec),sequence(spec['log_density']))
    provenance={'prepared':prepared,'specification':spec,'initializer':init,
        'convective_iterations':a.convective_iterations,'XH':a.hydrogen,'X3':a.helium3,
        'teff_K':a.teff,'log_g':a.logg,'initial':str(a.initial.resolve()),
        'initial_receipt_sha256':digest(a.initial/'completed.json'),
        'initial_validation_sha256':digest(a.initial/'validated.json'),'opacity_sha256':digest(a.opacity)}
    if changed_composition:
        provenance['initial_composition']={k:initial_record[k] for k in ['XH','X3']}
        provenance['composition_continuation']='starting guess only; target composition independently validated'
    a.work.mkdir(parents=True,exist_ok=True)
    receipt=a.work/'provenance.json'
    if receipt.exists() and json.loads(receipt.read_text())!=provenance:
        raise ValueError('continuation inputs changed; use a fresh directory')
    receipt.write_text(json.dumps(provenance,indent=2)+'\n')
    support={'temperature_K':temperatures(spec),'density_g_cm3':sequence(spec['log_density'])}
    files=[('log','run.log'),('convergence','fort.9'),('atmosphere_input','fort.5'),
           ('element_masses','ember-masses.dat'),('parameters','tas'),('initial_structure','fort.8')]
    for phase,executable in [('initializer',init['executable']),('final',prepared['tlusty'])]:
        d=a.work/phase
        atmosphere_inputs(d,prepared,spec,a.opacity,abundance,masses,a.teff,a.logg,text)
        if phase=='initializer':
            (d/'tas').write_text((d/'tas').read_text()+
                f'ICONRE={a.convective_iterations},ICONRS=1,IMUCON=200,IDEEPC=3,CRFLIM=-1\n')
        execute(executable,d,['fort.7','fort.9'])
        log=(d/'run.log').read_text()
        state=source_state(log,(d/'fort.9').read_text(),a.teff,a.logg,support,spec['tau'])
        if state['depths']!=spec['depths']:raise ValueError('atmosphere depth count differs')
        name='initializer.json' if phase=='initializer' else 'validated.json'
        record={k:provenance[k] for k in ['XH','X3','teff_K','log_g']}
        record['diagnostics']=state
        if phase=='final':
            source_inputs({k:(d/f).read_text() for k,f in files if k not in ['log','convergence']},
                          spec,a.hydrogen,a.helium3,a.teff,a.logg,log)
            for kind,file in files:
                target=archive(d/file,d)
                record[kind]=str(target.relative_to(a.work));record[kind+'_sha256']=digest(target)
            hashes={}
            for file in ['fort.5','fort.7','fort.8','fort.9','fort.15','tas','ember-masses.dat',
                         'opacity.sha256','physics.json','completed.json','run.log']:
                target=archive(a.work/'initializer'/file,a.work/'initializer')
                hashes[str(target.relative_to(a.work))]=digest(target)
            record['initialization']={'method':'pressure-preserving continuation and native CONREF, followed by canonical replay',
                                      'source':init,'files_sha256':hashes,'convective_iterations':a.convective_iterations}
            record['continuation_provenance']=provenance
        (d/name).write_text(json.dumps(record,indent=2)+'\n')
        text=(d/'fort.7').read_text()
        print(phase,json.dumps(state),flush=True)


if __name__=='__main__':main()
