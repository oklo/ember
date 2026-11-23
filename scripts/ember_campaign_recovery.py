"""Run finite, evidence-producing local jobs when a campaign lane has no runnable mass."""
import fcntl
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
from ember_column_restart import prepare, retry_job, seed

_workers = {}

def load(path):
    return json.loads(path.read_text())

def save(path,value):
    tmp=path.with_suffix('.new')
    tmp.write_text(json.dumps(value,indent=2)+'\n')
    tmp.replace(path)

def worker(directory,plan,lock_fd):
    with (directory/'log').open('a') as log:
        if 'seed' in plan:prepare(plan)
        proc=subprocess.Popen(plan['command'],stdout=log,stderr=log,start_new_session=True,
                              cwd=plan.get('cwd'),pass_fds=(lock_fd,),
                              env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
        try:
            code=proc.wait(timeout=plan['wall_seconds'])
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid,signal.SIGKILL)
            proc.wait()
            code=124
    passed=code==0 and all(Path(p).is_file() for p in plan['outputs'])
    save(directory/'result.json',dict(status='passed' if passed else 'failed',exit_code=code))

def external_live(job):
    pid=str(job['pid'])
    cmd=subprocess.run(['ps','-p',pid,'-o','command='],capture_output=True,text=True)
    if cmd.returncode==1 and not cmd.stdout.strip():return False
    cmd.check_returncode()
    if job['command_marker'] not in cmd.stdout:return False
    cwd=subprocess.run(['/usr/sbin/lsof','-a','-p',pid,'-d','cwd','-Fn'],text=True,capture_output=True)
    if cwd.returncode:raise RuntimeError('could not verify external worker directory')
    return 'n'+job['cwd'] in cwd.stdout.splitlines()

def attempt(root,job):
    d=root/'recovery'/job['name']
    if 'pid' in job:
        if external_live(job):return 'running'
        return 'passed' if all(Path(p).is_file() for p in job['outputs']) else 'failed'
    if not d.exists():return 'new'
    with (d/'lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return 'running'
    if (d/'result.json').is_file():return load(d/'result.json')['status']
    return 'interrupted'

def column_action(root,job,status):
    policy=job['resume_column']
    parent=Path(policy['directory'])
    record=root/'recovery'/job['name']/'selected-column.json'
    seen=set()
    current=job
    for number in range(1,policy['max_restarts']+2):
        if status in ('passed','interrupted'):
            spec=load(parent/'spec.json')
            receipt=parent/('hhe_arm_'+spec['anchor']+spec['tag']+'.json')
            if status=='interrupted' and receipt.is_file():status='passed'
        if status=='passed':
            r=load(receipt)
            if r['spec']!=spec or not r['quality']['usable']:return 'failed',None
            if not record.exists():
                record.parent.mkdir(parents=True,exist_ok=True)
                save(record,dict(directory=str(parent),receipt=str(receipt)))
            return 'passed',None
        if status not in ('failed','interrupted') or number>policy['max_restarts']:break
        previous=root/'recovery'/current['name']/'result.json'
        if previous.exists() and load(previous).get('exit_code') not in (-9,-15,-24,124,137,143,152):
            return status,None
        name=job['name']+f'-resume-{number:02}'
        receipt=root/'recovery'/name/'input.json'
        if receipt.exists():
            # Reconcile the existing child before examining its predecessor again.
            child=load(receipt)
        else:
            ck=seed(parent)
            if not ck or ck['profile_sha256'] in seen:return status,None
            child=retry_job(root,job,parent,ck,number)
            return 'ready',child
        seen.add(child['seed']['profile_sha256'])
        current=child
        parent=Path(child['column_directory']);status=attempt(root,child)
    return status,None

def tick(root, states, slots, allow):
    for pid in list(_workers):
        if _workers[pid].poll() is not None:
            del _workers[pid]
    specs=load(root/'recovery-plan.json')
    summary={};pending={}
    for job in specs:
        key=job['name']
        status=attempt(root,job)
        action=None
        if 'resume_column' in job and status!='new':
            try:
                status,action=column_action(root,job,status)
            except (ValueError,KeyError,TypeError,IndexError,OSError) as error:
                status='invalid column record'
                print(f'{key}: {error}',file=sys.stderr,flush=True)
        if status=='new':
            deps=[summary[k] for k in job.get('after',[])]
            if any(d in ('failed','interrupted','dependency failed','invalid column record') for d in deps):status='dependency failed'
            elif any(d!='passed' for d in deps):status='waiting'
            elif not any(s==job['trigger'] for s,_ in states.values()):status='unneeded'
            else:status,action='ready',job
        if action:pending[key]=action
        summary[key]=status
    active=sum(s=='running' for s in summary.values())
    for key,job in pending.items():
        if not allow or active>=slots:break
        directory=root/'recovery'/job['name']
        directory.mkdir(parents=True,exist_ok=False)
        save(directory/'input.json',job)
        with (directory/'lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            child=subprocess.Popen([sys.executable,__file__,str(directory),str(lock.fileno())],
                             pass_fds=(lock.fileno(),),start_new_session=True,stdin=subprocess.DEVNULL)
            _workers[child.pid]=child
        summary[key]='running';active+=1
    return active,summary

if __name__=='__main__':
    d=Path(sys.argv[1])
    try:
        worker(d,load(d/'input.json'),int(sys.argv[2]))
    except (ValueError,KeyError,TypeError,OSError) as error:
        save(d/'result.json',dict(status='failed',exit_code=1,error=str(error)))
