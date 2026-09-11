#!/usr/bin/env python3
"""Check an actual checkpoint continuation and its join to the preceding track.

The joined history retains all accepted rows without interpolation. The duplicate
restart row is removed only after its physical fields match the preceding endpoint
exactly; step diagnostics are absent at a restart. Source acceptance and executable
compatibility remain separate checks identified by the caller.
"""
import argparse
import datetime
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess

from summarize_remnant_endpoints import summarize


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['previous','continuation','output']:
        p.add_argument(name,type=Path)
    p.add_argument('--mixing-probe',type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--evidence',type=Path,nargs='+',required=True)
    p.add_argument('--joined-output',type=Path)
    a=p.parse_args()
    tracks=[];receipts=[];files={str(Path(__file__).resolve()):digest(__file__)}
    for path in [a.previous,a.continuation]:
        track=json.loads(path.read_text())
        receipt_path=path.with_suffix('.receipt.json');receipt=json.loads(receipt_path.read_text())
        if digest(path)!=receipt['output_sha256']:
            raise ValueError('stellar output checksum differs from receipt')
        if receipt['data_changed_during_run'] or receipt['restart_changed_during_run']:
            raise ValueError('physical input changed during evolution')
        exe=Path(receipt['command'][0])
        if digest(exe)!=receipt['executable_sha256']:
            raise ValueError('stellar executable checksum differs')
        for f in [path,receipt_path,path.with_suffix('.log'),exe]:files[str(f.resolve())]=digest(f)
        for f,expected in receipt.get('checkpoint_output_sha256',{}).items():
            if digest(f)!=expected:raise ValueError('output checkpoint changed')
            files[str(Path(f).resolve())]=expected
        if track['mass_Msun']!=.1 or track['points']!=512:
            raise ValueError('requires the 512-point 0.1 solar-mass calculation')
        for row in track['history']:
            for key,value in zip(track['columns'],row,strict=True):
                if value is None:
                    if key!='last_halfstep_gravothermal_Lsun':raise ValueError('missing history value')
                elif not math.isfinite(value):raise ValueError('nonfinite history')
        summarize(track)  # Also checks increasing ages and physical endpoint values.
        if len(track['profile'])!=512 or any(not all(math.isfinite(v) for v in row)
                or abs(sum(row[5:8])+.02-1)>1e-12 for row in track['profile']):
            raise ValueError('invalid final structure/composition')
        tracks.append(track);receipts.append(receipt)
    old,new=tracks;ro,rn=receipts
    if old['history'][0][0]!=0 or new['columns']!=old['columns']:
        raise ValueError('requires a complete original history with matching columns')
    start=dict(zip(new['columns'],new['history'][0],strict=True))
    endpoint=dict(zip(old['columns'],old['history'][-1],strict=True))
    step_only={'step_yr','step_error_norm','discrete_luminosity_balance','nuclear_rest_mass_balance',
               'last_halfstep_coupling_iterations','last_halfstep_gravothermal_Lsun'}
    for key,value in start.items():
        if key in step_only:
            if value!=(None if key=='last_halfstep_gravothermal_Lsun' else 0):
                raise ValueError('restart has unexpected step diagnostics')
        elif value!=endpoint[key]:raise ValueError('physical state changed at checkpoint join: '+key)
    checkpoint=Path(new['restart']['source'])
    if not checkpoint.is_absolute():checkpoint=Path(rn['working_directory'])/checkpoint
    checkpoint=checkpoint.resolve()
    old_checkpoints={str(Path(k).resolve()):v for k,v in ro['checkpoint_output_sha256'].items()}
    restart_inputs={str(Path(k).resolve()):v for k,v in rn['restart_input_sha256'].items()}
    if (str(checkpoint) not in old_checkpoints or restart_inputs.get(str(checkpoint))!=digest(checkpoint)
            or old_checkpoints[str(checkpoint)]!=digest(checkpoint)):
        raise ValueError('continuation is not linked to the preceding saved checkpoint')
    if (new['restart']['history_start_age_yr']!=endpoint['age_yr']
            or new['restart']['accepted_steps_before_restart']!=len(old['history'])-1):
        raise ValueError('restart age or accepted-step count does not match')
    # This first implementation permits only the checked atmosphere extension.
    if 'atmosphere_extension' not in new or 'eos_temperature_extension' in new or 'opacity_extension' in new:
        raise ValueError('expected a separately checked atmosphere-only extension')
    if digest(new['atmosphere_extension']['source_executable'])!=ro['executable_sha256']:
        raise ValueError('extension source executable differs from the original calculation')
    if 'nongrey:'+new['atmosphere_extension']['source_atmosphere']!=old['atmosphere_model']:
        raise ValueError('extension source atmosphere differs from the original selection')
    for key in ['mass_basis','nuclear_model','thermal_neutrino_model','transport_model','eos_model',
                'convection_criterion','step_error_tolerances','opacity_directory','atmosphere_tau_match']:
        if old[key]!=new[key]:raise ValueError('unexpected physical selection change: '+key)
    final=dict(zip(new['columns'],new['history'][-1],strict=True))
    if any(final[k]!=new['profile'][i][j] for k,i,j in
           [('central_X',0,5),('central_Y3',0,6),('central_T_K',0,3),('central_rho',0,2),('surface_X',-1,5)]):
        raise ValueError('final profile disagrees with history')
    log=a.continuation.with_suffix('.log').read_text()
    errors=[line for line in log.splitlines() if 'CompositionAtmosphereGrid:' in line and 'outside source grid' in line]
    if new['converged']:
        if rn['returncode']!=0:raise ValueError('inconsistent completion status')
    elif rn['returncode']!=1 or not errors or abs(final['Teff_K']-new['atmosphere_grid_support']['teff_K'][1])>1e-4:
        raise ValueError('unclassified stopping condition')
    profile=''.join(' '.join(format(v,'.17g') for v in row)+'\n' for row in new['profile'])
    response=subprocess.run([str(a.mixing_probe),new['eos_model'].split(':',1)[1],new['opacity_directory']],
                            input=profile,text=True,capture_output=True,check=True)
    regions=[json.loads(line) for line in response.stdout.splitlines()]
    if abs(sum(r['mass_fraction'] for r in regions if r['convective'])-final['convective_mass_fraction'])>1e-12:
        raise ValueError('independent mixing classification disagrees')
    files[str(a.mixing_probe.resolve())]=digest(a.mixing_probe)
    for p in a.evidence:files[str(p.resolve())]=digest(p)
    joined={**new,'history':old['history']+new['history'][1:],
            'history_assembly':{'scope':__doc__,'sources':[str(a.previous),str(a.continuation)],
                'source_sha256':[digest(a.previous),digest(a.continuation)],'physical_join_exact':True}}
    if a.joined_output:
        if a.joined_output.exists():raise FileExistsError('use a new joined output')
        a.joined_output.write_text(json.dumps(joined,separators=(',',':'),allow_nan=False)+'\n')
        files[str(a.joined_output.resolve())]=digest(a.joined_output)
    a.archive.mkdir(parents=True,exist_ok=False)
    archives=[]
    for i,(path,checksum) in enumerate(files.items()):
        raw=Path(path).read_bytes();packed=gzip.compress(raw,mtime=0)
        if hashlib.sha256(raw).hexdigest()!=checksum:raise ValueError('input changed during audit')
        target=a.archive/(f'{i:03d}-'+Path(path).name+'.gz');target.write_bytes(packed)
        archives.append({'source':path,'archive':str(target),'sha256':checksum,'gzip_sha256':digest(target)})
    result={'scope':__doc__,'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'endpoint':final,'continued_states':len(new['history']),'joined_states':len(joined['history']),
            'physical_join_exact':True,'duplicate_initial_diagnostics_excluded':sorted(step_only),
            'checkpoint_sha256':digest(checkpoint),'mixing_regions':regions,'milestones':summarize(joined),
            'requested_age_reached':new['converged'],'first_atmosphere_domain_rejection':errors[0] if errors else None,
            'last_atmosphere_domain_rejection':errors[-1] if errors else None,
            'segment_cpu_seconds':rn['child_user_seconds']+rn['child_system_seconds'],
            'segment_awake_seconds':rn['awake_elapsed_seconds'],'segment_utc_seconds':rn['utc_elapsed_seconds'],
            'segment_rejected_attempts':new['rejected_steps']-new['restart']['rejected_steps_before_restart'],
            'all_attempts_including_preceding_terminal_rejections':old['rejected_steps']+new['rejected_steps']-new['restart']['rejected_steps_before_restart'],
            'input_sha256':files,'recovery_archives':archives}
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(f"Checked {len(joined['history'])} joined states through {final['age_yr']/1e12:.4g} trillion years.")


if __name__=='__main__':main()
