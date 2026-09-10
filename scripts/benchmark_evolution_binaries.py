#!/usr/bin/env python3
"""Compare immutable evolution executables with identical inputs and ABBA timing.

Every JSON output must match byte for byte, including all accepted timesteps
and final composition/structure. Child CPU time excludes concurrent source
jobs, although scheduling and processor frequency can still affect timing.
"""
import argparse
import json
from pathlib import Path
import platform
import resource
import shutil
import subprocess
import time

from run_evolution_snapshot import ROOT, sha, input_data, data_label


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--after-step-workers',type=int,choices=[1,2],
                   help='change only the candidate executable CPU-worker setting')
    for name in ['before','after','work','output']:p.add_argument(name,type=Path)
    p.add_argument('arguments',nargs=argparse.REMAINDER);a=p.parse_args()
    args=a.arguments[1:] if a.arguments[:1]==['--'] else a.arguments
    if any(arg.startswith('--checkpoint') or arg=='--restart' for arg in args):
        raise ValueError('benchmark uses independent fresh starts without checkpoints')
    if a.work.exists():raise FileExistsError('use a fresh benchmark directory')
    a.work.mkdir(parents=True);a.work=a.work.resolve()
    executables={}
    for name,path in [('before',a.before),('after',a.after)]:
        target=a.work/('ember-evolve-'+name);shutil.copy2(path,target);executables[name]=target
    hashes={data_label(f):sha(f) for f in input_data(args)}
    binary_hashes={name:sha(path) for name,path in executables.items()}
    records=[];reference=None
    for i,name in enumerate(['before','after','after','before']):
        start=time.perf_counter();cpu=resource.getrusage(resource.RUSAGE_CHILDREN)
        command=[str(executables[name]),*args]
        if name=='after' and a.after_step_workers is not None:
            command.extend(['--step-workers',str(a.after_step_workers)])
        result=subprocess.run(command,cwd=ROOT,capture_output=True)
        end=resource.getrusage(resource.RUSAGE_CHILDREN)
        (a.work/f'run-{i}.json').write_bytes(result.stdout);(a.work/f'run-{i}.log').write_bytes(result.stderr)
        data=json.loads(result.stdout)
        if result.returncode or not data.get('converged'):raise RuntimeError('stellar benchmark failed')
        if reference is None:reference=result.stdout
        if result.stdout!=reference:raise ValueError('executable change alters the accepted numerical trajectory')
        record={'binary':name,'command':command,'wall_seconds':time.perf_counter()-start,
                'child_cpu_seconds':end.ru_utime+end.ru_stime-cpu.ru_utime-cpu.ru_stime,
                'output_sha256':sha(a.work/f'run-{i}.json'),'accepted_steps':len(data['history'])-1,
                'age_yr':data['history'][-1][0]}
        records.append(record);print(json.dumps(record),flush=True)
    if any(sha(ROOT/name)!=value for name,value in hashes.items()):raise ValueError('benchmark data changed')
    if any(sha(executables[name])!=value for name,value in binary_hashes.items()):raise ValueError('benchmark binary changed')
    means={name:sum(r['child_cpu_seconds'] for r in records if r['binary']==name)/2 for name in executables}
    wall_means={name:sum(r['wall_seconds'] for r in records if r['binary']==name)/2 for name in executables}
    report={'scope':__doc__,'architecture':platform.machine(),'arguments':args,
            'data_sha256':hashes,'executables_sha256':binary_hashes,'runs':records,
            'mean_child_cpu_seconds':means,'cpu_speedup_before_over_after':means['before']/means['after'],
            'mean_wall_seconds':wall_means,'wall_speedup_before_over_after':wall_means['before']/wall_means['after'],
            'all_outputs_byte_identical':True,'script_sha256':sha(Path(__file__))}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
