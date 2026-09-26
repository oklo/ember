#!/usr/bin/env python3
"""Build a bounded reusable grid of the checked Born/Pauli pair matrices.

Each completed source is written once and retained independently. Generation
does not accept an interpolator or select stellar transport coefficients.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import time

from electron_electron_energy_modes_v2 import energy_mode_matrix


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def calculate(task):
    key,args,path=task
    start=time.process_time();wall=time.monotonic()
    result=energy_mode_matrix(**args)
    result.update(source_key=key,cpu_seconds=time.process_time()-start,
                  elapsed_seconds=time.monotonic()-wall,completed_utc=datetime.now(timezone.utc).isoformat())
    with Path(path).open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    return dict(key=key,path=path,sha256=digest(path),cpu_seconds=result['cpu_seconds'],
                elapsed_seconds=result['elapsed_seconds'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('specification',type=Path)
    args=parser.parse_args();spec=json.loads(args.specification.read_text())
    work=Path(spec['work']);assert work.is_dir()
    target=Path(spec['report']);assert not target.exists()
    for path,sha in spec['input_sha256'].items():assert digest(path)==sha,path
    source_dir=work/'sources';source_dir.mkdir()
    tasks=[]
    for i,eta in enumerate(spec['eta_grid']):
        for j,b in enumerate(spec['b_thermal_grid']):
            key=f'e{i:02d}-b{j:02d}'
            params=dict(eta=eta,b_thermal=b,degree=9,orders=(80,40,40,48,48),tail=50.,normalization_tail=80.)
            tasks.append((key,params,str(source_dir/(key+'.json'))))
    completed=[];failed=[]
    with ProcessPoolExecutor(max_workers=spec['workers']) as pool:
        futures={pool.submit(calculate,t):t for t in tasks}
        for future in as_completed(futures):
            key=futures[future][0]
            try:
                entry=future.result();completed.append(entry)
                print(json.dumps(dict(completed=len(completed),total=len(tasks),**entry)),flush=True)
            except Exception as error:
                failed.append(dict(key=key,error=repr(error)))
                print(json.dumps(failed[-1]),flush=True)
            progress=dict(completed=completed,failed=failed,total=len(tasks))
            temp=work/'progress.json.tmp';temp.write_text(json.dumps(progress,indent=2)+'\n');temp.replace(work/'progress.json')
    result=dict(outcome='completed_unvalidated_table' if not failed else 'incomplete',
                accepted_for_stellar_evolution=False,eta_grid=spec['eta_grid'],b_thermal_grid=spec['b_thermal_grid'],
                mode_count=10,completed_sources=sorted(completed,key=lambda x:x['key']),failed_sources=failed,
                total_source_cpu_seconds=sum(s['cpu_seconds'] for s in completed),
                specification=str(args.specification),specification_sha256=digest(args.specification),
                input_sha256=spec['input_sha256'],
                limitations=['Prescribed static Born/Pauli collision physics; no relativistic or correlated-mixture corrections.',
                             'Interpolation has not been independently tested or selected.',
                             'The domain covers the saved hot profile parameters, not the full future track or cold remnant.'])
    with target.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(dict(outcome=result['outcome'],sources=len(completed),source_cpu_seconds=result['total_source_cpu_seconds'],
                          report=str(target),sha256=digest(target))),flush=True)
    if failed:raise SystemExit(1)


if __name__=='__main__':main()
