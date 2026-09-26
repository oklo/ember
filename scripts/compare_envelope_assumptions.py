#!/usr/bin/env python3
"""Measure envelope sensitivities at the retained 100 Myr conditional state."""
import argparse
import hashlib
import json
import os
import selectors
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from audit_full_star_timestep import compare

YEAR = 365.25*86400


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('output', type=Path)
    ap.add_argument('--scratch', type=Path, required=True)
    ap.add_argument('--experiment', choices=('envelope','cn'), default='envelope')
    args = ap.parse_args()
    assert not args.output.exists() and not args.scratch.exists()
    args.scratch.mkdir()
    identities = {}

    def pin(p):
        p = Path(p); data = p.read_bytes()
        identities[str(p)] = hashlib.sha256(data).hexdigest()
        return data

    prior = json.loads(pin('docs/results/full_star_adaptive_diffusion_v7.json'))
    assert prior['target_reached'] and not prior['failures']
    for p,h in prior['artifacts_sha256'].items():
        assert hashlib.sha256(pin(p)).hexdigest() == h, p
    parent = Path('/tmp/ember-adaptive-diffusion-v7')
    initial = json.loads(pin(parent/'checkpoint.json'))['model_record']
    inputs = json.loads(pin(parent/'inputs.json'))
    for p,h in inputs['input_sha256'].items():
        if Path(p).suffix == '.dat': assert hashlib.sha256(pin(p)).hexdigest() == h,p
    build = Path('/tmp/ember-cn-stellar-build-v1' if args.experiment=='cn' else '/tmp/ember-envelope-sensitivity-build-v1')
    for p,h in json.loads(pin(build/'manifest.json')).items():
        assert hashlib.sha256(pin(p)).hexdigest() == h,p
    command = [str(build/'star_probe'), *inputs['command'][1:-1], 'select-cn' if args.experiment=='cn' else 'select-envelope']
    pin(command[0]);pin(__file__)
    # factor: bridge opacity; micro: hot reduced microscopic carried heat only.
    # The prescribed cool heat is a conservative stress test, not an error bound.
    cases = [dict(name='baseline',factor=1.,lower=2e6,upper=3e6,micro=1.,stress=0.),
             dict(name='opacity_minus',factor=.97,lower=2e6,upper=3e6,micro=1.,stress=0.),
             dict(name='opacity_plus',factor=1.03,lower=2e6,upper=3e6,micro=1.,stress=0.),
             dict(name='microscopic_heat_zero',factor=1.,lower=2e6,upper=3e6,micro=0.,stress=0.),
             dict(name='microscopic_heat_double',factor=1.,lower=2e6,upper=3e6,micro=2.,stress=0.),
             dict(name='hotter_join',factor=1.,lower=2.5e6,upper=3e6,micro=1.,stress=0.),
             dict(name='cool_heat_minus',factor=1.,lower=2e6,upper=3e6,micro=1.,stress=-.001),
             dict(name='cool_heat_plus',factor=1.,lower=2e6,upper=3e6,micro=1.,stress=.001)]
    if args.experiment=='cn':
        cases=[dict(name=name,factor=1.,carbon=carbon) for name,carbon in (('baseline',-1),('nitrogen_only',0),('CN',1))]
    (args.scratch/'inputs.json').write_text(json.dumps(dict(command=command,cases=cases,
        initial_checkpoint=str(parent/'checkpoint.json'),interval_years=5e6,input_sha256=identities),indent=2)+'\n')
    a = np.asarray(initial['model']);m = a[:,0]
    w = np.zeros(512);w[0] = m[0];w[:-1] += .5*np.diff(m);w[1:] += .5*np.diff(m)
    queries=(args.scratch/'queries.txt').open('x');raw=(args.scratch/'responses.jsonl').open('x')
    stderr=(args.scratch/'stderr.txt').open('x')
    start=time.monotonic()
    proc=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=stderr,text=True,bufsize=1)
    (args.scratch/'native_process.json').write_text(json.dumps(dict(pid=proc.pid,controller_pid=os.getpid(),
        command=command,start_utc=datetime.now(timezone.utc).isoformat(),maximum_elapsed_seconds=720),indent=2)+'\n')
    sel=selectors.DefaultSelector();sel.register(proc.stdout,selectors.EVENT_READ)
    steps=[];results={};failures=[];cache={}

    def step(old,dt,c,part):
        prefix=[c['carbon']] if args.experiment=='cn' else [c['lower'],c['upper'],c['micro'],c['stress']]
        vals=[*prefix,1,c['factor'],dt,1e-10,1e-10,old['age_seconds'],1,512]
        line=' '.join(format(x,'.17g') for x in vals)+' '+' '.join(format(x,'.17g') for x in np.asarray(old['model_log']).flat)+'\n'
        if line in cache:return cache[line]
        assert time.monotonic()-start < 720
        queries.write(line);queries.flush();proc.stdin.write(line);proc.stdin.flush()
        while True:
            if not sel.select(120):raise TimeoutError('No native response within 120 seconds')
            text=proc.stdout.readline()
            if not text:raise RuntimeError('Native process ended')
            raw.write(text);raw.flush();r=json.loads(text)
            if 'error' in r:raise RuntimeError(r['error'])
            if r.get('kind')=='preflight':assert r['supported_nodes']==512 and not r['errors'],r
            elif r.get('kind')=='step':break
        compact={k:v for k,v in r.items() if k not in ('model','model_log','nuclear_dXdt','total_species_rates','mixing_regions')}
        compact.update(case=c['name'],part=part);steps.append(compact)
        if not r['converged']:raise RuntimeError(c['name']+': '+r['message'])
        b=np.asarray(r['model']);flux=np.zeros((513,2));flux[1:-1]=r['total_species_rates']
        storage=w[:,None]*(b[:,5:7]-np.asarray(old['model'])[:,5:7]-dt*np.asarray(r['nuclear_dXdt']))/m[-1]
        compact.update(local_species_balance=float(abs(storage+dt/m[-1]*np.diff(flux,axis=0)).max()),
                       global_species_balance=float(abs(storage.sum(axis=0)).max()))
        assert max(compact['local_species_balance'],compact['global_species_balance'])<=1e-13
        assert abs(r['luminosity_balance'])<=2e-8 and abs(r['nuclear_mass_balance'])<=2e-6
        assert r['input_preserved'] and r['initial_age_seconds']==old['age_seconds'] and r['age_seconds']==old['age_seconds']+dt
        assert np.isfinite(b).all() and np.array_equal(b[:,0],m) and np.min(b[:,5:7])>=0 and np.max(b[:,5:7].sum(axis=1))<=.98
        cache[line]=r
        print(json.dumps(compact),flush=True)
        return r

    try:
        for c in cases:
            try:
                full=step(initial,5e6*YEAR,c,'full')
                first=step(initial,2.5e6*YEAR,c,'first_half')
                fine=step(first,2.5e6*YEAR,c,'second_half')
                x,y=np.asarray(full['model_log']),np.asarray(fine['model_log'])
                delta=x[:,5:7]-y[:,5:7]
                species=np.column_stack((abs(delta),abs(delta.sum(axis=1))))
                fractions=np.maximum(abs(np.column_stack((x[:,5:7],.98-x[:,5:7].sum(axis=1)))),
                                     abs(np.column_stack((y[:,5:7],.98-y[:,5:7].sum(axis=1)))))
                errors=dict(log_structure=float(abs(x[:,1:4]-y[:,1:4]).max()/1e-4),
                            scaled_abundance=float((species/(1e-8+1e-4*fractions)).max()),
                            surface_luminosity=float(abs(x[-1,4]/y[-1,4]-1)/1e-3),
                            hydrogen_mass=float(w@abs(delta[:,0])/(w@y[:,5])/1e-6))
                regions=lambda s:[v for v in s['mixing_regions'] if v[1]>v[0]+1]
                assert max(errors.values())<=1 and regions(full)==regions(fine),errors
                results[c['name']]=dict(time_error_terms=errors,coarse=full,fine=fine)
                (args.scratch/(c['name']+'-endpoint.json')).write_text(json.dumps(fine)+'\n')
            except Exception as exc:
                failures.append(dict(case=c['name'],message=str(exc)))
                if c['name']=='baseline':break
        proc.stdin.close();proc.wait(timeout=10);assert proc.returncode==0
    finally:
        if proc.poll() is None:
            proc.terminate()
            try:proc.wait(timeout=3)
            except subprocess.TimeoutExpired:proc.kill();proc.wait(timeout=3)
        sel.close();queries.close();raw.close();stderr.close()
    comparisons={}
    if 'baseline' in results:
        base=results['baseline']
        for name,v in results.items():
            comparisons[name]=dict(time_error_terms=v['time_error_terms'],
                full_step=compare(base['coarse'],v['coarse']),two_half_steps=compare(base['fine'],v['fine']),
                relative_total_hydrogen=float(w@np.asarray(v['fine']['model'])[:,5]/(w@np.asarray(base['fine']['model'])[:,5])-1),
                same_mixed_regions=base['fine']['mixing_regions']==v['fine']['mixing_regions'])
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed' if not failures else 'incomplete',
                experiment=args.experiment,
                initial_conditional_elapsed_years=1e8,additional_interval_years=5e6,comparisons=comparisons,
                native_solves=len(steps),solver_seconds=sum(s['seconds'] for s in steps),wall_seconds=time.monotonic()-start,
                steps=steps,failures=failures,native_exit=proc.returncode,
                artifacts_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in args.scratch.iterdir()},
                accepted_for_stellar_evolution=False,
                limitations=['The imposed cool flux is a prescribed conservative stress test, not a bound or a microscopic transport law.',
                             'Changing the heat join changes conductivity and reduced microscopic heat together.',
                             'This measures the finite 5 Myr response at one evolved state, not a full-track sensitivity.'])
    if args.experiment=='cn':
        report['limitations']=['Closed CN cycle with fixed GS98 catalyst number; no historical carbon-conversion fuel is restored.',
                              'The nitrogen-only case is a catalyst sensitivity; it is not the inferred catalyst history.',
                              'This measures the finite 5 Myr response at one evolved state, not a full-track sensitivity.']
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('outcome','native_solves','solver_seconds','failures')},indent=2))
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
