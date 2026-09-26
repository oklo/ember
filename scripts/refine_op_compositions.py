#!/usr/bin/env python3
"""Refine an existing OP composition family without repeating completed mixtures."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.interpolate import PchipInterpolator
from generate_op_warm_family import plane_task, connected_states, atomic_json


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def key(row):return row['temperature_index'],row['X'],row['Z']
def plane_path(work,task):
    return work/'planes'/('mixture_'+hashlib.sha256(json.dumps(task,separators=(',',':')).encode()).hexdigest()+'.json')


def table_files(work,rows,config,step=.025):
    folder=work/'density0025';folder.mkdir()
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
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['family','comparison','work','output']:parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args();work=args.work
    if args.output.exists():raise FileExistsError(args.output)
    work.mkdir(exist_ok=True);(work/'planes').mkdir(exist_ok=True)
    inputs={str(p):digest(p) for p in [Path(__file__),args.family,args.comparison]}
    base=json.loads(args.family.read_text());comparison=json.loads(args.comparison.read_text())
    assert base['outcome']=='assembled_refined_conditional_family'
    for p,h in base['input_sha256'].items():assert digest(p)==h;inputs[p]=h
    config_path=Path(base['configuration']);assert digest(config_path)==base['configuration_sha256']
    oldconfig=json.loads(config_path.read_text());inputs[str(config_path)]=digest(config_path)
    for p,h in comparison['input_sha256'].items():assert digest(p)==h
    comparison_config=Path(comparison['configuration'])
    assert digest(comparison_config)==comparison['configuration_sha256']
    direct_path=comparison_config.parent/'withheld_composition_means.json'
    direct_hash=digest(direct_path);assert direct_hash==comparison['input_sha256'][str(direct_path)]
    inputs[str(direct_path)]=direct_hash
    xs=set(oldconfig['X']);rules=[]
    for a,b in zip(oldconfig['X'],oldconfig['X'][1:]):
        selected=[r for r in comparison['interval_errors'] if r['X_lower']==a and r['X_upper']==b]
        assert len(selected)>=len(oldconfig['Z'])
        error=max(r['maximum_relative_difference'] for r in selected)
        depth=3 if error>.02 else 1 if error>.005 else 0
        segment=[a,b]
        for _ in range(depth):segment=sorted(set(segment+[.5*(l+h) for l,h in zip(segment,segment[1:])]))
        xs.update(segment);rules.append(dict(X_lower=a,X_upper=b,previous_maximum_error=error,subintervals=2**depth))
    config=dict(X=sorted(xs),Z=oldconfig['Z'],temperatures=oldconfig['temperatures'],
        logR_min=oldconfig['logR_min'],logR_max=oldconfig['logR_max'],
        physical_configuration=oldconfig['physical_configuration'],
        physical_configuration_sha256=oldconfig['physical_configuration_sha256'],
        refinement_rules=rules,input_sha256=inputs.copy())
    new_config=work/'configuration.json'
    if new_config.exists():assert json.loads(new_config.read_text())==config
    else:atomic_json(new_config,config)
    @lru_cache(maxsize=None)
    def retained_file(path,expected):
        assert digest(path)==expected,path
        return json.loads(Path(path).read_text())
    known={};origins={}
    for record in base['source_records']:
        origin=record['origin'];stored=retained_file(origin['path'],origin['sha256'])
        row=stored[origin['index']] if 'index' in origin else stored
        task=tuple(record['task']);assert key(row)==task
        known[task]=row;origins[task]=origin
    for index,row in enumerate(retained_file(str(direct_path),direct_hash)):
        task=key(row)
        if task in known:assert row==known[task]
        else:
            known[task]=row;origins[task]=dict(path=str(direct_path),sha256=direct_hash,index=index,kind='retained_independent_mixture')
    tasks=[(it,x,z) for it in config['temperatures'] for z in config['Z'] for x in config['X']]
    assert len({plane_path(work,t) for t in tasks})==len(tasks)
    rows=[];records=[];pending=[]
    for task in tasks:
        p=plane_path(work,task)
        if task in known:
            rows.append(known[task]);records.append(dict(task=task,origin=origins[task]))
        elif p.exists():
            row=json.loads(p.read_text());assert key(row)==task
            rows.append(row);records.append(dict(task=task,origin=dict(path=str(p),sha256=digest(p),kind='resumed_completed_plane')))
        else:pending.append(task)
    reused=len(rows);start=time.monotonic()
    print(json.dumps(dict(X_planes=len(config['X']),total=len(tasks),reused=reused,new=len(pending))),flush=True)
    with ProcessPoolExecutor(max_workers=8) as pool:
        futures={pool.submit(plane_task,t):t for t in pending}
        for future in as_completed(futures):
            task=futures[future];row=future.result();assert key(row)==task
            p=plane_path(work,task);atomic_json(p,row);rows.append(row)
            records.append(dict(task=task,origin=dict(path=str(p),sha256=digest(p),kind='new_plane')))
            if (len(rows)-reused)%36==0 or len(rows)==len(tasks):
                progress=dict(completed=len(rows),total=len(tasks),elapsed_seconds=time.monotonic()-start,
                    failed_native_states=sum(len(r['failures']) for r in rows))
                atomic_json(work/'progress.json',progress);print(json.dumps(progress),flush=True)
    rows.sort(key=lambda r:(r['temperature_index'],r['Z'],r['X']));records.sort(key=lambda r:r['task'])
    unsupported=[]
    for row in rows:
        try:connected_states(row)
        except ValueError as e:unsupported.append(dict(task=key(row),reason=str(e)))
    files=table_files(work,rows,config) if not unsupported else []
    for p,h in inputs.items():assert digest(p)==h,p
    report=dict(outcome='assembled_refined_conditional_family' if not unsupported else 'incomplete_native_support',
        accepted_for_stellar_opacity=False,input_sha256=inputs,configuration=str(new_config),
        configuration_sha256=digest(new_config),X_planes=len(config['X']),planes=len(rows),
        reused_planes=reused,new_planes=len(pending),source_records=records,table_files=files,
        unsupported_planes=unsupported,native_states=sum(len(r['native_states']) for r in rows),
        failed_native_states=sum(len(r['failures']) for r in rows),
        limitations=['Additional composition nodes require new independent comparisons before acceptance.',
            'Atomic physics, missing elements, metallicity and thermodynamic source grids are unchanged.',
            'Independent mixtures incorporated as table nodes are no longer independent validation data.',
            'No selected stellar opacity, atmosphere, history or checkpoint is modified.'])
    atomic_json(args.output,report)
    print(json.dumps(dict(report=str(args.output),sha256=digest(args.output),outcome=report['outcome'],
        reused=reused,new=len(pending))),flush=True)


if __name__=='__main__':main()
