#!/usr/bin/env python3
"""Run ember-evolve with a copied executable, input hashes and explicit clocks.

The copied executable is immutable during this run. Data are hashed before
and after execution. On macOS, perf_counter may exclude system suspension;
UTC elapsed time and child user/system CPU seconds are recorded separately.
"""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import shlex
import shutil
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def input_data(arguments):
    """Include selected families and their references, even outside data/."""
    data={f.resolve() for directory in ['eos','opacity','conduction','atmosphere']
          for f in (ROOT/'data'/directory).glob('*.dat')}
    for arg in arguments:
        if arg.startswith(('nongrey:', 'metal:')):
            path=Path(arg.split(':',1)[1])
            if not path.is_absolute():path=ROOT/path
            path=path.resolve(strict=True);data.add(path)
            if arg.startswith('metal:'):data.update(f.resolve() for f in path.parent.glob('*.dat'))
    if '--opacity-directory' in arguments:
        directory=Path(arguments[arguments.index('--opacity-directory')+1])
        if not directory.is_absolute():directory=ROOT/directory
        for name in ['aesopus21_gs98_mixture.dat','tops_gs98_mixture_low.dat','tops_gs98_mixture_high.dat']:
            manifest=(directory/name).resolve(strict=True);data.add(manifest)
            rows=manifest.read_text().splitlines();header=rows[0].split()
            if len(header)!=5 or header[:2]!=['EMBER_OPACITY_MIXTURE','1'] or len(rows)-1!=int(header[2]):
                raise ValueError('invalid selected opacity manifest')
            for row in rows[1:]:
                fields=shlex.split(row)
                if len(fields)!=2:raise ValueError('invalid opacity source path')
                data.add((manifest.parent/fields[1]).resolve(strict=True))
    return sorted(data)


def data_label(path):
    try:return str(path.relative_to(ROOT))
    except ValueError:return str(path)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable',type=Path);p.add_argument('work',type=Path);p.add_argument('output',type=Path)
    p.add_argument('arguments',nargs=argparse.REMAINDER)
    a=p.parse_args();arguments=a.arguments[1:] if a.arguments[:1]==['--'] else a.arguments
    if len(arguments)<2:raise ValueError('supply ember-evolve mesh, duration and physical selections')
    duration=float(arguments[1])
    if not duration>0:raise ValueError('invalid duration')
    if a.work.exists():raise FileExistsError('use a fresh snapshot directory')
    if a.output.exists() or a.output.with_suffix('.log').exists():raise FileExistsError('output already exists')
    a.work.mkdir(parents=True);a.output.parent.mkdir(parents=True,exist_ok=True)
    executable=a.work.resolve()/'ember-evolve';shutil.copy2(a.executable,executable)
    data=input_data(arguments)
    sources=[f for directory in ['src','include','apps','examples'] for f in (ROOT/directory).rglob('*')
             if f.is_file() and f.suffix in ['.cpp','.hpp','.txt']]
    data_hashes={data_label(f):sha(f) for f in data}
    source_hashes={str(f.relative_to(ROOT)):sha(f) for f in sorted(sources)}
    restart_inputs={}
    if '--restart' in arguments:
        path=Path(arguments[arguments.index('--restart')+1]).resolve(strict=True)
        restart_inputs[str(path)]=sha(path)
    command=[str(executable),*arguments]
    receipt={'command':command,'working_directory':str(ROOT),'executable_sha256':sha(executable),
             'data_sha256':data_hashes,'source_sha256':source_hashes,
             'timing_note':'UTC elapsed includes suspension; perf_counter may exclude suspension on macOS. CPU sums this process\'s child user and system times. Concurrent source jobs are excluded from CPU seconds.'}
    if restart_inputs:receipt['restart_input_sha256']=restart_inputs
    record=a.work/'receipt.json'
    receipt['started_unix']=time.time();awake_start=time.perf_counter();cpu_start=resource.getrusage(resource.RUSAGE_CHILDREN)
    record.write_text(json.dumps(receipt,indent=2)+'\n')
    with a.output.open('w') as output,a.output.with_suffix('.log').open('w') as log:
        process=subprocess.Popen(command,cwd=ROOT,stdout=output,stderr=log)
        (a.work/'running.json').write_text(json.dumps({'pid':process.pid,'started_unix':receipt['started_unix'],'command':command},indent=2)+'\n')
        code=process.wait()
    receipt['finished_unix']=time.time();receipt['awake_elapsed_seconds']=time.perf_counter()-awake_start
    cpu_end=resource.getrusage(resource.RUSAGE_CHILDREN)
    receipt.update({'utc_elapsed_seconds':receipt['finished_unix']-receipt['started_unix'],
                    'child_user_seconds':cpu_end.ru_utime-cpu_start.ru_utime,
                    'child_system_seconds':cpu_end.ru_stime-cpu_start.ru_stime,
                    'returncode':code,'output_sha256':sha(a.output)})
    receipt['data_changed_during_run']=[name for name,value in data_hashes.items() if sha(ROOT/name)!=value]
    receipt['source_changed_during_run']=[name for name,value in source_hashes.items() if sha(ROOT/name)!=value]
    receipt['restart_changed_during_run']=[name for name,value in restart_inputs.items() if sha(Path(name))!=value]
    if '--checkpoint' in arguments:
        path=Path(arguments[arguments.index('--checkpoint')+1]).resolve()
        receipt['checkpoint_output_sha256']={str(path):sha(path)} if path.exists() else {}
    record.write_text(json.dumps(receipt,indent=2)+'\n');(a.work/'running.json').unlink()
    a.output.with_suffix('.receipt.json').write_text(record.read_text())
    if code or receipt['data_changed_during_run'] or receipt['restart_changed_during_run']:
        raise RuntimeError('run failed or input data changed; retained receipt and outputs')
    result=json.loads(a.output.read_text())
    if not result['converged'] or result['history'][-1][0]!=duration:raise RuntimeError('requested age was not reached')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['data_sha256','source_sha256']},indent=2))


if __name__=='__main__':main()
