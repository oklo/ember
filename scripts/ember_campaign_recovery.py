"""Run finite, evidence-producing local jobs when a campaign lane has no runnable mass."""
import fcntl
from graphlib import TopologicalSorter
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

def retryable_column_failure(result_path, policy):
    if not result_path.is_file():
        return False
    code = load(result_path).get('exit_code')
    if code in (-9, -15, -24, 124, 137, 143, 152):
        return True
    if code != 1 or not policy.get('retry_detached_convection'):
        return False
    log = result_path.with_name('log')
    if not log.is_file():
        return False
    with log.open('rb') as stream:
        stream.seek(0, os.SEEK_END)
        stream.seek(max(0, stream.tell()-4096))
        tail = stream.read().rstrip().splitlines()
    return bool(tail) and tail[-1] in (
        b'ValueError: The current column requires detached regions or a different discretization',
        b'ValueError: Single-region search cycles; no column accepted',
    )

def column_action(root,job,status):
    policy=job['resume_column']
    parent=Path(policy['directory'])
    record=root/'recovery'/job['name']/'selected-column.json'
    seen=set()
    current=job
    for number in range(1,policy['max_restarts']+2):
        previous=root/'recovery'/current['name']/'result.json'
        if status=='failed':
            if not retryable_column_failure(previous,policy):
                return status,None
        if status in ('passed','failed','interrupted'):
            spec=load(parent/'spec.json')
            receipt=parent/('hhe_arm_'+spec['anchor']+spec['tag']+'.json')
            if receipt.is_file():status='passed'
        if status=='passed':
            r=load(receipt)
            if r.get('spec')!=spec or r.get('quality',{}).get('usable') is not True:return 'failed',None
            if not record.exists():
                record.parent.mkdir(parents=True,exist_ok=True)
                save(record,dict(directory=str(parent),receipt=str(receipt)))
            return 'passed',None
        if status not in ('failed','interrupted') or number>policy['max_restarts']:break
        name=job['name']+f'-resume-{number:02}'
        receipt=root/'recovery'/name/'input.json'
        if receipt.exists():
            # Reconcile the existing child before examining its predecessor again.
            child=load(receipt)
        else:
            ck=seed(parent)
            if (not ck or ck['profile_sha256'] in seen) and policy.get('alternative_seed'):
                ck=seed(parent,policy['alternative_seed'])
            if not ck or ck['profile_sha256'] in seen:return status,None
            child=retry_job(root,job,parent,ck,number)
            return 'ready',child
        seen.add(child['seed']['profile_sha256'])
        current=child
        parent=Path(child['column_directory']);status=attempt(root,child)
    return status,None

def column_job(job):
    return ('resume_column' in job or job.get('resource') == 'atmosphere'
            or any(Path(str(arg)).name == 'm03_midpoint_if_needed.py' for arg in job.get('command', [])))


def gate_allows(root, gate):
    if gate is None:
        return False
    path, field = gate
    source = root / path
    return source.is_file() and load(source).get(field) is True


def tick(root, states, slots, allow, release_routing=None, expand_plan=None):
    for pid in list(_workers):
        if _workers[pid].poll() is not None:
            del _workers[pid]
    specs=load(root/'recovery-plan.json')
    jobs={job['name']:job for job in specs}
    if len(jobs)!=len(specs):
        raise ValueError('duplicate job name in recovery plan')
    if expand_plan:
        specs = expand_plan(root, specs)
        jobs = {job['name']:job for job in specs}
    dependencies={key:job.get('after',[]) for key,job in jobs.items()}
    unknown={dep for deps in dependencies.values() for dep in deps if dep not in jobs}
    if unknown:
        raise ValueError(f'unknown recovery prerequisites: {sorted(unknown)}')
    priorities, gates = release_routing(root, jobs) if release_routing else ({}, {})
    summary={};pending={}
    for key in TopologicalSorter(dependencies).static_order():
        job=jobs[key]
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
            if any(d in ('failed','interrupted','dependency failed','invalid column record','await scientific assessment') for d in deps):status='dependency failed'
            elif any(d!='passed' for d in deps):status='waiting'
            elif key in gates and not gate_allows(root, gates[key]):status='await scientific assessment'
            elif not any(s==job['trigger'] for s,_ in states.values()):status='unneeded'
            else:status,action='ready',job
        if action:pending[key]=action
        summary[key]=status
    # Dispatch in the operator's priority order after resolving dependencies.
    pending={key:pending[key] for key in sorted(jobs, key=lambda name: priorities.get(name, len(jobs))) if key in pending}
    active=sum(s=='running' for s in summary.values())
    columns={job['name'] for job in specs if column_job(job)}
    column_running=any(summary[key]=='running' for key in columns)
    while allow and pending and active<slots:
        key=next(iter(pending))
        if column_running and key in columns:
            # Give a ready evolution or assessment a lane while a long atmosphere
            # solve is active. A second column can still use an otherwise idle lane.
            key=next((name for name in pending if name not in columns),key)
        job=pending.pop(key)
        directory=root/'recovery'/job['name']
        directory.mkdir(parents=True,exist_ok=False)
        save(directory/'input.json',job)
        with (directory/'lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            child=subprocess.Popen([sys.executable,__file__,str(directory),str(lock.fileno())],
                             pass_fds=(lock.fileno(),),start_new_session=True,stdin=subprocess.DEVNULL)
            _workers[child.pid]=child
        summary[key]='running';active+=1
        column_running |= key in columns
    return active,summary

if __name__=='__main__':
    d=Path(sys.argv[1])
    try:
        worker(d,load(d/'input.json'),int(sys.argv[2]))
    except (ValueError,KeyError,TypeError,OSError) as error:
        save(d/'result.json',dict(status='failed',exit_code=1,error=str(error)))
