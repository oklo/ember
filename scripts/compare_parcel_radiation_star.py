#!/usr/bin/env python3
"""Compare optional parcel radiation heat from a retained conditional checkpoint."""
import argparse
import json
import os
import selectors
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from audit_full_star_diffusion import digest
from audit_initial_diffusion_step import time_error, YEAR
from audit_full_star_timestep import compare


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--scratch',type=Path,required=True)
    args=parser.parse_args()
    assert not args.output.exists() and not args.scratch.exists()
    args.scratch.mkdir()
    root=Path(__file__).resolve().parents[1]
    prior=root/'docs/results/full_star_adaptive_diffusion_v1.json'
    prior_report=json.loads(prior.read_text())
    assert not prior_report['failures']
    identities={str(prior):digest(prior)}
    for p,h in prior_report['artifacts_sha256'].items():
        assert digest(p)==h,p
        identities[p]=h
    checkpoint=Path('/tmp/ember-adaptive-diffusion-v1/checkpoint.json')
    initial=json.loads(checkpoint.read_text())['model_record']
    assert initial['converged'] and initial['maximum_species_face']==396
    old_inputs=json.loads(Path('/tmp/ember-adaptive-diffusion-v1/inputs.json').read_text())
    for p,h in old_inputs['input_sha256'].items():
        if Path(p).suffix=='.dat':
            assert digest(p)==h,p
            identities[p]=h
    source_manifest=Path('/tmp/ember-full-star-timestep-source-v2/manifest.json')
    changed=[]
    allowed={'include/ember/eos_smooth_mixture.hpp','src/eos_smooth_mixture.cpp',
             'include/ember/screened_microscopic_transport.hpp','src/screened_microscopic_transport.cpp',
             'scripts/conditional_envelope_heat.hpp','scripts/full_star_diffusion_step_probe.cpp'}
    for original,r in json.loads(source_manifest.read_text()).items():
        assert digest(r['retained'])==r['sha256'],original
        if digest(original)!=r['sha256']:
            relative=str(Path(original).relative_to(root));assert relative in allowed,relative;changed.append(relative)
        identities[original]=digest(original)
    command=['/tmp/ember-parcel-star-probe-v2',*old_inputs['command'][1:],'select-radiation']
    audit=root/'docs/results/native_parcel_heat_v2.json'
    assert not json.loads(audit.read_text())['failures']
    for p in [Path(__file__).resolve(),root/'scripts/audit_initial_diffusion_step.py',root/'scripts/audit_full_star_timestep.py',
              Path(command[0]),Path('/tmp/ember-parcel-radiation-build-v1/src/libember.a'),source_manifest,audit]:
        identities[str(p)]=digest(p)
    inputs=dict(command=command,input_sha256=identities,changed_shared_sources=changed,
                initial_checkpoint=str(checkpoint),interval_years=100000,
                scope='One common retained state; material heat alone versus optional LTE parcel radiation heat.')
    (args.scratch/'inputs.json').write_text(json.dumps(inputs,indent=2)+'\n')
    stamps={p:(Path(p).stat().st_size,Path(p).stat().st_mtime_ns) for p in identities}
    mass=np.asarray(initial['model'])[:,0];weights=np.zeros(512);weights[0]=mass[0]
    weights[:-1]+=.5*np.diff(mass);weights[1:]+=.5*np.diff(mass)
    assert abs(weights.sum()/mass[-1]-1)<1e-14
    logs=[(args.scratch/name).open('x') for name in ('queries.txt','responses.jsonl','stderr.txt')]
    queries,raw,err=logs;start=time.monotonic()
    proc=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=err,text=True,bufsize=1)
    (args.scratch/'native_process.json').write_text(json.dumps(dict(pid=proc.pid,controller_pid=os.getpid(),command=command,
           start_utc=datetime.now(timezone.utc).isoformat(),maximum_elapsed_seconds=480),indent=2)+'\n')
    selector=selectors.DefaultSelector();selector.register(proc.stdout,selectors.EVENT_READ)
    failures=[];steps=[];endpoints={};accuracy={};coarse_endpoints={}
    def step(state,dt,radiation,label):
        model=np.asarray(state['model_log']);age=state['age_seconds']
        line=f'{radiation} 1 1 {dt:.17g} 1e-14 1e-10 {age:.17g} 1 512 '+' '.join(format(x,'.17g') for x in model.flat)+'\n'
        queries.write(line);queries.flush();proc.stdin.write(line);proc.stdin.flush();record=None
        while record is None:
            if not selector.select(120):raise TimeoutError('No native response within 120 seconds.')
            text=proc.stdout.readline()
            if not text:raise RuntimeError('Native process ended early.')
            raw.write(text);raw.flush();r=json.loads(text)
            if 'error' in r:raise RuntimeError(r['error'])
            if r.get('kind')=='preflight':assert r['supported_nodes']==512 and not r['errors'],r
            elif r.get('kind')=='step':record=r
            else:raise RuntimeError('Unknown native result kind.')
        compact={k:v for k,v in record.items() if k not in ('model','model_log','total_species_rates','nuclear_dXdt','mixing_regions','extra_radiation_luminosity')}
        compact['case']=label;steps.append(compact)
        assert record['converged'],record['message']
        assert record['input_preserved'] and record['parcel_radiation']==radiation
        assert record['initial_age_seconds']==age and record['age_seconds']==age+dt
        new=np.asarray(record['model']);logs=np.asarray(record['model_log'])
        assert new.shape==logs.shape==(512,7) and np.isfinite(new).all() and np.isfinite(logs).all()
        assert np.array_equal(new[:,0],mass) and np.min(new[:,1:4])>0
        assert np.min(new[:,5:7])>=0 and np.max(new[:,5:7].sum(axis=1))<=.98
        rates=np.zeros((513,2));rates[1:-1]=record['total_species_rates'];source=np.asarray(record['nuclear_dXdt'])
        storage=weights[:,None]*(new[:,5:7]-np.asarray(state['model'])[:,5:7]-dt*source)/mass[-1]
        local=storage+dt/mass[-1]*np.diff(rates,axis=0)
        compact.update(independent_local_species_balance=float(abs(local).max()),
                       independent_global_species_balance=float(abs(storage.sum(axis=0)).max()),
                       mixed_regions=[r for r in record['mixing_regions'] if r[1]>r[0]+1])
        assert abs(record['luminosity_balance'])<=2e-8 and abs(record['nuclear_mass_balance'])<=2e-6
        assert record['abundance_residual']<=1e-14 and record['material_heat_residual']<=1e-10
        assert max(compact['independent_local_species_balance'],compact['independent_global_species_balance'])<=1e-13
        for lo,hi in compact['mixed_regions']:assert np.max(abs(new[lo:hi,5:7]-new[lo,5:7]))==0
        extra=np.asarray(record['extra_radiation_luminosity']);assert extra.shape==(511,) and np.isfinite(extra).all()
        i=int(np.argmax(abs(extra)));boundary=np.zeros(513);boundary[1:-1]=extra
        compact.update(maximum_extra_luminosity_erg_s=float(abs(extra[i])),extra_luminosity_face=i,
                       maximum_extra_over_surface_L=float(abs(extra[i])/new[-1,4]),
                       maximum_extra_divergence_over_surface_L=float(abs(np.diff(boundary)).max()/new[-1,4]),
                       summed_extra_divergence_over_surface_L=float(np.diff(boundary).sum()/new[-1,4]))
        assert record['conductivity_difference']==0
        assert abs(compact['summed_extra_divergence_over_surface_L'])<1e-14
        print(json.dumps(compact),flush=True)
        return record
    try:
        for radiation in (0,1):
            coarse=step(initial,100000*YEAR,radiation,'coarse')
            first=step(initial,50000*YEAR,radiation,'half-1')
            fine=step(first,50000*YEAR,radiation,'half-2')
            accuracy[radiation]=time_error(coarse,fine);endpoints[radiation]=fine;coarse_endpoints[radiation]=coarse
        proc.stdin.close();proc.wait(timeout=10);assert proc.returncode==0
    except Exception as exc:failures.append(dict(exception=type(exc).__name__,message=str(exc)))
    finally:
        if proc.poll() is None:
            proc.terminate()
            try:proc.wait(timeout=3)
            except subprocess.TimeoutExpired:proc.kill();proc.wait(timeout=3)
        selector.close()
        for stream in logs:stream.close()
    for p,stamp in stamps.items():
        if (Path(p).stat().st_size,Path(p).stat().st_mtime_ns)!=stamp and digest(p)!=identities[p]:failures.append(dict(changed_input=p))
    effect={}
    if len(endpoints)==2:
        effect=dict(two_half_steps=compare(endpoints[0],endpoints[1]),one_full_step=compare(coarse_endpoints[0],coarse_endpoints[1]))
        for key,value in endpoints.items():
            (args.scratch/f'checkpoint-radiation-{key}.json').write_text(json.dumps(dict(accepted_for_stellar_evolution=False,
                initial_checkpoint=str(checkpoint),additional_elapsed_seconds=100000*YEAR,model_record=value))+'\n')
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='completed_conditional_parcel_radiation_comparison' if not failures else 'failed_conditional_parcel_radiation_comparison',
                accepted_for_stellar_evolution=False,selected_star_unchanged=True,initial_checkpoint=str(checkpoint),
                additional_interval_years=100000,new_steps=len(steps),steps=steps,time_accuracy=accuracy,physical_response=effect,
                elapsed_seconds=time.monotonic()-start,native_exit=proc.returncode,failures=failures,
                artifacts_sha256={str(p):digest(p) for p in args.scratch.iterdir()},
                limitations=['LTE, isobaric parcel enthalpy uses linear composition contrasts; its correlation assumptions still need assessment.',
                             'Cool microscopic drift and boundary motion remain conditional; opacity sampling remains unaccepted.',
                             'This finite-interval response does not bound long-term effects or extend the selected trajectory.'])
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('outcome','new_steps','physical_response','failures')},indent=2),flush=True)
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
