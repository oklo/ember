#!/usr/bin/env python3
"""Refine OP hydrogen sampling using retained independent comparisons.

The first family, its physical evaluator and all completed mixtures remain
unchanged. Only missing composition/isotherm calculations are performed.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.interpolate import PchipInterpolator

from generate_op_warm_family import plane_task, connected_states, atomic_json

ROOT=Path(__file__).resolve().parents[1]
WORK=Path('/tmp/ember-op-warm-family-v2')
BASE=ROOT/'docs/results/op_warm_family_v1.json'
COMPARISON=ROOT/'docs/results/op_warm_runtime_v1.json'
OUTPUT=ROOT/'docs/results/op_warm_family_v2.json'


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def key(row):return row['temperature_index'],row['X'],row['Z']


def filename(task):
    body=json.dumps(task,separators=(',',':')).encode()
    return WORK/'planes'/('mixture_'+hashlib.sha256(body).hexdigest()+'.json')


def table_files(rows,config,step=.025):
    folder=WORK/'density0025';folder.mkdir()
    xs,zs,temps=config['X'],config['Z'],config['temperatures']
    lo,hi=config['logR_min'],config['logR_max']
    grid=lo+step*np.arange(round((hi-lo)/step)+1)
    lookup={key(r):r for r in rows};files=[]
    for component in ['ordinary','absorption','scattering']:
        for z in zs:
            lines=[f'EMBER_OPACITY_TABLE 2 {len(xs)} {len(temps)} {len(grid)} OP_GS98_17_elements_finite_frequency_{component}',
                   ' '.join(f'{v:.17g}' for v in grid),
                   ' '.join(f'{.025*it:.17g}' for it in temps)]
            for x in xs:
                lines.append(f'{x:.17g} {z:.17g}')
                for it in temps:
                    states=connected_states(lookup[it,x,z])
                    rr=np.array([s['logR'] for s in states])
                    count=int(np.searchsorted(grid,rr[-1],side='right'))
                    if count<4:raise ValueError('insufficient source support')
                    values=PchipInterpolator(rr,np.log10([s[component] for s in states]),extrapolate=False)(grid[:count])
                    if not np.isfinite(values).all():raise ValueError('native density extrapolation')
                    lines.append(str(count)+' '+' '.join(f'{v:.17g}' for v in values))
            path=folder/f'op_gs98_z{round(z*1000):03d}_{component}.dat'
            with path.open('x') as f:f.write('\n'.join(lines)+'\n')
            files.append(dict(path=str(path),sha256=digest(path),bytes=path.stat().st_size,Z=z,component=component))
        path=folder/f'op_gs98_{component}.dat'
        lines=[f'EMBER_OPACITY_MIXTURE 1 {len(zs)} logR OP_GS98_17_elements_finite_frequency_{component}']
        lines += [f'{z:.17g} "op_gs98_z{round(z*1000):03d}_{component}.dat"' for z in zs]
        path.write_text('\n'.join(lines)+'\n')
        files.append(dict(path=str(path),sha256=digest(path),bytes=path.stat().st_size,component=component))
    return files


def main():
    if OUTPUT.exists():raise FileExistsError(OUTPUT)
    WORK.mkdir(exist_ok=True);(WORK/'planes').mkdir(exist_ok=True)
    base=json.loads(BASE.read_text());comparison=json.loads(COMPARISON.read_text())
    inputs={str(p):digest(p) for p in [Path(__file__),BASE,COMPARISON,
        ROOT/'scripts/generate_op_warm_family.py']}
    for p,h in base['input_sha256'].items():
        if digest(p)!=h:raise ValueError('source calculation changed: '+p)
        inputs[p]=h
    config_path=Path(base['configuration'])
    assert digest(config_path)==base['configuration_sha256']
    oldconfig=json.loads(config_path.read_text());inputs[str(config_path)]=digest(config_path)
    direct=Path('/tmp/ember-op-warm-runtime-v1/withheld_composition_means.json')
    direct_hash=digest(direct)
    assert direct_hash==comparison['input_sha256'][str(direct)]
    inputs[str(direct)]=direct_hash
    oldx=oldconfig['X'];xs=set(oldx);rules=[]
    for a,b in zip(oldx,oldx[1:]):
        selected=[r for r in comparison['heldout_supported'] if r['Z']==.02 and r['X']==.5*(a+b)]
        if not selected:raise ValueError('missing interval comparison')
        error=max(abs(r['relative_opacity_difference']) for r in selected)
        depth=3 if error>.04 else 2 if error>.02 else 1 if error>.005 else 0
        segment=[a,b]
        for _ in range(depth):
            segment=sorted(set(segment+[.5*(l+h) for l,h in zip(segment,segment[1:])]))
        xs.update(segment)
        rules.append(dict(X_lower=a,X_upper=b,previous_maximum_error=error,subintervals=2**depth))
    config=dict(X=sorted(xs),Z=oldconfig['Z'],temperatures=oldconfig['temperatures'],
                logR_min=oldconfig['logR_min'],logR_max=oldconfig['logR_max'],
                physical_configuration=str(config_path),physical_configuration_sha256=digest(config_path),
                refinement_rules=rules,input_sha256=inputs)
    new_config=WORK/'configuration.json'
    if new_config.exists():assert json.loads(new_config.read_text())==config,'resumption inputs changed'
    else:atomic_json(new_config,config)
    known={};origins={}
    for f in base['plane_files']:
        assert digest(f['path'])==f['sha256'];row=json.loads(Path(f['path']).read_text())
        known[key(row)]=row;origins[key(row)]={**f,'kind':'original_family'}
    for index,row in enumerate(json.loads(direct.read_text())):
        if key(row) in known:assert known[key(row)]==row
        known[key(row)]=row;origins[key(row)]=dict(path=str(direct),sha256=direct_hash,index=index,kind='retained_independent_mixture')
    tasks=[(it,x,z) for it in config['temperatures'] for z in config['Z'] for x in config['X']]
    assert len(set(filename(t) for t in tasks))==len(tasks)
    rows=[];records=[];pending=[]
    for task in tasks:
        if task in known:
            rows.append(known[task]);records.append(dict(task=task,origin=origins[task]));continue
        p=filename(task)
        if p.exists():
            row=json.loads(p.read_text());assert key(row)==task
            rows.append(row);records.append(dict(task=task,origin=dict(path=str(p),sha256=digest(p),kind='resumed_completed_plane')))
        else:pending.append(task)
    reused=len(rows);start=time.monotonic();done=0
    print(json.dumps(dict(X_planes=len(config['X']),total=len(tasks),reused=reused,new=len(pending))),flush=True)
    with ProcessPoolExecutor(max_workers=8) as pool:
        futures={pool.submit(plane_task,t):t for t in pending}
        for future in as_completed(futures):
            task=futures[future];row=future.result();assert key(row)==task
            p=filename(task);atomic_json(p,row);rows.append(row)
            records.append(dict(task=task,origin=dict(path=str(p),sha256=digest(p),kind='new_plane')));done+=1
            if done%36==0 or done==len(pending):
                progress=dict(completed=len(rows),total=len(tasks),elapsed_seconds=time.monotonic()-start,
                              failed_native_states=sum(len(r['failures']) for r in rows))
                atomic_json(WORK/'progress.json',progress);print(json.dumps(progress),flush=True)
    rows.sort(key=lambda r:(r['temperature_index'],r['Z'],r['X']));records.sort(key=lambda r:r['task'])
    unsupported=[]
    for r in rows:
        try:connected_states(r)
        except ValueError as e:unsupported.append(dict(task=key(r),reason=str(e)))
    files=table_files(rows,config) if not unsupported else []
    for p,h in inputs.items():assert digest(p)==h,p
    report=dict(outcome='assembled_refined_conditional_family' if not unsupported else 'incomplete_native_support',
                accepted_for_stellar_opacity=False,input_sha256=inputs,configuration=str(new_config),
                configuration_sha256=digest(new_config),X_planes=len(config['X']),planes=len(rows),
                reused_planes=reused,new_planes=len(pending),source_records=records,table_files=files,
                unsupported_planes=unsupported,native_states=sum(len(r['native_states']) for r in rows),
                failed_native_states=sum(len(r['failures']) for r in rows),
                limitations=['Refinement is driven by previous independent errors; new independent comparisons remain required.',
                    'Metallicity coordinates and thermodynamic grids are unchanged; atomic physics and missing-element assumptions are unchanged.',
                    'Former independent mixtures incorporated as table nodes are no longer withheld validation data.',
                    'No new stellar model, atmosphere or selected opacity input is created.'])
    atomic_json(OUTPUT,report)
    print(json.dumps(dict(report=str(OUTPUT),sha256=digest(OUTPUT),outcome=report['outcome'],reused=reused,new=len(pending))),flush=True)


if __name__=='__main__':main()
