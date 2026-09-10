#!/usr/bin/env python3
"""Continue a verified gas atmosphere, then replay the canonical source.

Native convective refinement is an optional initializer only. Both source
runs retain their input decks, complete outputs and independent checks.
"""
import argparse
import json
from pathlib import Path
from generate_nongrey_grid import (atmosphere_inputs,composition,execute,temperatures,
    sequence,completed_initial_structure,continuation_structure,resample_initial_structure,archive,
    truncate_initial_structure)
from import_nongrey_grid import source_inputs,source_state,read_text
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
    p.add_argument('--convective-iterations',type=int,default=200,
                   help='native refinement through this global iteration; zero solves the canonical source directly from the verified seed')
    p.add_argument('--initializer-depths',type=int,help='coarser trial mesh; final source always uses the specification resolution')
    p.add_argument('--initializer-frequencies',type=int,help='coarser trial frequency sampling; final canonical replay is unchanged')
    p.add_argument('--initial-bottom-tau',type=float,
                   help='truncate the verified seed at this measured optical depth; requires independent lower-boundary convergence checks')
    p.add_argument('--composition-continuation',action='store_true',
                   help='allow a verified starting guess at another H/He composition; target chemistry and opacity are solved anew')
    a=p.parse_args()
    cancellation=a.work.parent/'cancellation.json'
    if cancellation.exists():
        cancelled=json.loads(cancellation.read_text())
        if cancelled.get('work')!=str(a.work.parent.resolve()) or cancelled.get('plan_sha256')!=digest(a.work.parent/'plan.json'):
            raise ValueError('invalid enclosing plan cancellation record')
        raise RuntimeError('enclosing source plan cancelled before launch: '+cancelled['reason'])
    prepared=json.loads(a.prepared.read_text());spec=json.loads(a.specification.read_text())
    init=json.loads(a.initializer.read_text())
    if 'depletion' in prepared or 'initialization_only' in prepared:
        raise ValueError('canonical gas source required')
    if digest(prepared['tlusty'])!=prepared['executables']['tlusty']:
        raise ValueError('canonical executable changed')
    if digest(init['executable'])!=init['executable_sha256'] or digest(
            Path(prepared['tlusty']).parent/'tlusty208.f')!=init['base_source_sha256']:
        raise ValueError('convective initializer identity mismatch')
    if not 0<=a.convective_iterations<=200:raise ValueError('invalid refinement count')
    trial_spec={**spec,'depths':spec['depths'] if a.initializer_depths is None else a.initializer_depths,
                'atmosphere_frequencies':spec['atmosphere_frequencies'] if a.initializer_frequencies is None else a.initializer_frequencies}
    if not 20<=trial_spec['depths']<=spec['depths'] or not 2<=trial_spec['atmosphere_frequencies']<=spec['atmosphere_frequencies']:
        raise ValueError('invalid initializer resolution')
    coarse=a.convective_iterations>0 and trial_spec!=spec
    initial_record=json.loads((a.initial/'validated.json').read_text())
    changed_composition=initial_record['XH']!=a.hydrogen or initial_record['X3']!=a.helium3
    if changed_composition and not a.composition_continuation:
        raise ValueError('continuation requires the same composition')
    original=completed_initial_structure(a.initial,spec,initial_record['teff_K'],
                                         initial_record['log_g'],spec['depths'])
    if original is None:raise ValueError('initial model exceeds requested depth count')
    if a.initial_bottom_tau is not None:
        if not a.initial_bottom_tau > spec['tau']:
            raise ValueError('seed bottom must lie deeper than the matching depth')
        original=truncate_initial_structure(original,read_text(a.initial/'run.log.gz')
                    if (a.initial/'run.log.gz').exists() else (a.initial/'run.log').read_text(),a.initial_bottom_tau)
        original=resample_initial_structure(original,spec['depths'])
    text=continuation_structure(original,initial_record['teff_K'],initial_record['log_g'],a.teff,a.logg)
    if coarse:text=resample_initial_structure(text,trial_spec['depths'],allow_coarsen=True)
    abundance,masses=composition(a.hydrogen,a.helium3,spec['metals'])
    validate_table(a.opacity,abundance,temperatures(spec),sequence(spec['log_density']))
    provenance={'prepared':prepared,'specification':spec,'initializer':init,
        'convective_iterations':a.convective_iterations,'XH':a.hydrogen,'X3':a.helium3,
        'teff_K':a.teff,'log_g':a.logg,'initial':str(a.initial.resolve()),
        'initial_receipt_sha256':digest(a.initial/'completed.json'),
        'initial_validation_sha256':digest(a.initial/'validated.json'),'opacity_sha256':digest(a.opacity)}
    if coarse:provenance['initializer_resolution']={k:trial_spec[k] for k in ['depths','atmosphere_frequencies']}
    if coarse:provenance['initializer_sequence']='coarse native convection, full-resolution native convection, canonical replay'
    if a.initial_bottom_tau is not None:
        provenance['initial_bottom_tau']=a.initial_bottom_tau
        provenance['lower_boundary_note']='seed truncated at measured Rosseland depth; final optical depth changes during convergence and must be recorded'
    if changed_composition:
        provenance['initial_composition']={k:initial_record[k] for k in ['XH','X3']}
        provenance['composition_continuation']='starting guess only; target composition independently validated'
    a.work.mkdir(parents=True,exist_ok=True)
    receipt=a.work/'provenance.json'
    if receipt.exists() and json.loads(receipt.read_text())!=provenance:
        raise ValueError('continuation inputs changed; use a fresh directory')
    if (a.work/'final/validated.json').exists():
        # Completed runs archive their input decks. Rewriting those decks and
        # invoking execute would discard the expensive accepted calculation.
        from assemble_nongrey_grid import load_continuation
        _, state, _, _, _ = load_continuation(a.work)
        print('reused canonical', json.dumps(state), flush=True)
        return
    receipt.write_text(json.dumps(provenance,indent=2)+'\n')
    support={'temperature_K':temperatures(spec),'density_g_cm3':sequence(spec['log_density'])}
    files=[('log','run.log'),('convergence','fort.9'),('atmosphere_input','fort.5'),
           ('element_masses','ember-masses.dat'),('parameters','tas'),('initial_structure','fort.8')]
    phases=[('initializer',init['executable'])] if a.convective_iterations else []
    if coarse:phases.append(('refined_initializer',init['executable']))
    phases.append(('final',prepared['tlusty']))
    for phase,executable in phases:
        d=a.work/phase
        phase_spec=trial_spec if phase=='initializer' else spec
        atmosphere_inputs(d,prepared,phase_spec,a.opacity,abundance,masses,a.teff,a.logg,text)
        if phase!='final':
            (d/'tas').write_text((d/'tas').read_text()+
                f'ICONRE={a.convective_iterations},ICONRS=1,IMUCON=200,IDEEPC=3,CRFLIM=-1\n')
        execute(executable,d,['fort.7','fort.9'])
        log=(d/'run.log').read_text()
        state=source_state(log,(d/'fort.9').read_text(),a.teff,a.logg,support,spec['tau'],
                           max_flux_error=.05 if phase=='initializer' and coarse else .002)
        if state['depths']!=phase_spec['depths']:raise ValueError('atmosphere depth count differs')
        name='validated.json' if phase=='final' else 'initializer.json'
        record={k:provenance[k] for k in ['XH','X3','teff_K','log_g']}
        record['diagnostics']=state
        if phase=='final':
            source_inputs({k:(d/f).read_text() for k,f in files if k not in ['log','convergence']},
                          spec,a.hydrogen,a.helium3,a.teff,a.logg,log)
            for kind,file in files:
                target=archive(d/file,d)
                record[kind]=str(target.relative_to(a.work));record[kind+'_sha256']=digest(target)
            hashes={}
            for initial_phase,_ in phases[:-1]:
                for file in ['fort.5','fort.7','fort.8','fort.9','fort.15','tas','ember-masses.dat',
                             'opacity.sha256','physics.json','completed.json','run.log']:
                    target=archive(a.work/initial_phase/file,a.work/initial_phase)
                    hashes[str(target.relative_to(a.work))]=digest(target)
            record['initialization']={'method':('pressure-preserving continuation and native CONREF, followed by canonical replay'
                                              if a.convective_iterations else 'pressure-preserving continuation followed directly by canonical solve'),
                                      'source':init,'files_sha256':hashes,'convective_iterations':a.convective_iterations}
            record['continuation_provenance']=provenance
        (d/name).write_text(json.dumps(record,indent=2)+'\n')
        text=(d/'fort.7').read_text()
        print(phase,json.dumps(state),flush=True)


if __name__=='__main__':main()
