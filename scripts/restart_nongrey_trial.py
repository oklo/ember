#!/usr/bin/env python3
"""Reuse an unfinished atmosphere as a guess; accept only a new canonical solve."""
import argparse
import json
from pathlib import Path

from generate_nongrey_grid import (atmosphere_inputs, composition, execute,
    temperatures, sequence, archive, resample_initial_structure)
from import_nongrey_grid import source_inputs, source_state
from nongrey_opacity import validate_table
from prepare_nongrey_sources import digest
from assemble_nongrey_grid import load_continuation


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trial',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--temperature-step-limit',type=float)
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError('retain previous trials; choose a fresh output directory')
    if (args.trial/'final/validated.json').exists():
        raise ValueError('use the validated continuation path for a completed source')
    parent=json.loads((args.trial/'provenance.json').read_text())
    prepared=parent['prepared'];spec=dict(parent['specification'])
    if 'depletion' in prepared or 'initialization_only' in prepared:
        raise ValueError('canonical gas source required')
    executable=Path(prepared['tlusty'])
    if digest(executable)!=prepared['executables']['tlusty']:
        raise ValueError('canonical executable changed')
    if args.temperature_step_limit is not None:
        if not 1<args.temperature_step_limit<=1.25:
            raise ValueError('invalid temperature correction limit')
        spec['temperature_step_limit']=args.temperature_step_limit
    trial_checkpoint=args.trial/'final/fort.7'
    guess=trial_checkpoint.read_text()
    # Checks the full finite, positive, ordered structure. This does not
    # assert radiative/convective equilibrium or accept it as a table node.
    guess=resample_initial_structure(guess,spec['depths'])
    opacity=(args.trial/'final/opacity.bin').resolve(strict=True)
    if digest(opacity)!=parent['opacity_sha256']:
        raise ValueError('trial opacity changed')
    h,y,t,g=[parent[k] for k in ['XH','X3','teff_K','log_g']]
    abundance,masses=composition(h,y,spec['metals'])
    validate_table(opacity,abundance,temperatures(spec),sequence(spec['log_density']))
    args.output.mkdir()
    seed=args.output/'trial_seed';seed.mkdir()
    for name,path in [('fort.7',trial_checkpoint),('parent_provenance.json',args.trial/'provenance.json')]:
        (seed/name).write_bytes(path.read_bytes())
    initialization_files={str(p.relative_to(args.output)):digest(p) for p in seed.iterdir()}
    provenance=dict(prepared=prepared,specification=spec,XH=h,X3=y,teff_K=t,log_g=g,
        opacity_sha256=digest(opacity),trial_parent=str(args.trial.resolve()),
        trial_checkpoint_sha256=digest(trial_checkpoint),
        parent_provenance_sha256=digest(args.trial/'provenance.json'),
        initialization='Unfinished iterate is an initial guess only. New canonical solve and all independent source checks are required.',
        runner_sha256=digest(Path(__file__)))
    (args.output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    final=args.output/'final'
    atmosphere_inputs(final,prepared,spec,opacity,abundance,masses,t,g,guess)
    execute(executable,final,['fort.7','fort.9'])
    log=(final/'run.log').read_text()
    support=dict(temperature_K=temperatures(spec),density_g_cm3=sequence(spec['log_density']))
    diagnostics=source_state(log,(final/'fort.9').read_text(),t,g,support,spec['tau'])
    if diagnostics['depths']!=spec['depths']:
        raise ValueError('source depth count changed')
    files=[('log','run.log'),('convergence','fort.9'),('atmosphere_input','fort.5'),
           ('element_masses','ember-masses.dat'),('parameters','tas'),('initial_structure','fort.8')]
    source_inputs({kind:(final/name).read_text() for kind,name in files if kind not in ['log','convergence']},spec,h,y,t,g,log)
    record=dict(XH=h,X3=y,teff_K=t,log_g=g,diagnostics=diagnostics,
                initialization=dict(method=provenance['initialization'],files_sha256=initialization_files),
                continuation_provenance=provenance)
    for kind,name in files:
        path=archive(final/name,final)
        record[kind]=str(path.relative_to(args.output));record[kind+'_sha256']=digest(path)
    (final/'validated.json').write_text(json.dumps(record,indent=2)+'\n')
    key,state,_,_,_=load_continuation(args.output)
    print(json.dumps(dict(coordinates=key,state=state)),flush=True)


if __name__=='__main__':
    main()
