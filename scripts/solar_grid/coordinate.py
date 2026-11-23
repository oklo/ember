"""Two local lanes for the true-solar sequence. No provider calls."""
import fcntl
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import time
import sys

scripts = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(scripts))
root = Path(os.environ.get('EMBER_GRID_ROOT', scripts.parent / 'out' / 'solar-grid-z01648-nov23')).resolve()
out=root.parent
reference=out/'solar-boundary-m0070-z01648-nov23'
solver=out/'boundary-corrected-item1-oct13/layer_refinement_oct14/ember-layer-refined'
campaign=json.loads((root/'campaign.json').read_text())
solar=campaign['initial_Z']
masses=[6,8,5,9,4,10,*[i for i in range(1,31) if i not in (4,5,6,7,8,9,10)]]
def read(p):return json.loads(p.read_text())
def write(p,obj):
    tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(obj,indent=2)+'\n');tmp.replace(p)
def last(p):
    with p.open('rb') as f:
        f.seek(0,2);f.seek(max(f.tell()-32768,0))
        for line in reversed(f.read().splitlines()):
            try:return json.loads(line)
            except ValueError:continue
    return {}
def dirs(m):
    return sorted((root/f'm{m:02}').glob('s[0-9][0-9][0-9][0-9]'))
def locked(path):
    with path.open('a') as f:
        try:fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return True
    return False
def alternate_birth(m):
    # Reuse only numerical starting guesses from an accepted same-mass birth.
    cases=out/'solar-1k-grid-nov17'/f'm{10*m:04}'/'cases'
    for f in sorted(cases.glob('*/evolution/history.jsonl')):
        with f.open() as stream:
            first=stream.readline();second=stream.readline()
        if not second:continue
        if json.loads(first)['years']==0 and json.loads(second)['years']>0:
            config=f.parents[1]/'config.txt'
            if config.is_file():return config
    return None

def next_step(m):
    cases=dirs(m)
    if not cases:return 'birth',None
    case=cases[-1]
    if locked(case/'worker.lock'):return 'running',case
    evo=case/'evolution';p=evo/'report.json'
    if len(cases)==1 and not (evo/'latest.checkpoint').is_file():
        if alternate_birth(m):return 'alternate birth guess',case
    if not p.is_file():
        ck=evo/'latest.checkpoint'
        if ck.is_file() and Path(str(ck)+'.lithium').is_file() and (evo/'history.jsonl').is_file():
            lines=ck.read_text().splitlines();idx=9+int(lines[8])
            state=lines[idx].split()
            row=last(evo/'history.jsonl')
            if (lines[-1]=='END' and row.get('accepted',0)>0 and int(state[4])==row['accepted']
                    and abs(float(state[2])/31557600-row['years'])<.1
                    and row['years']>read(case/'launch.json')['start_age_years']):
                return 'resume interrupted checkpoint',case
        if not (evo/'history.jsonl').exists() and (case/'run.log').is_file():
            tail=(case/'run.log').read_text()[-2000:]
            if 'pressure-matched atmosphere' in tail and not any(read(c/'launch.json').get('action')=='atmosphere birth control' for c in cases):
                return 'atmosphere birth control',case
            if 'atmosphere' in tail.lower():return 'review birth atmosphere',case
        return 'inspect interrupted run',case
    r=read(p);reason=r['stop_reason'].lower()
    if not (evo/'history.jsonl').exists():return 'inspect initialization',case
    v=last(evo/'history.jsonl')
    if v.get('Teff_K',math.inf)<=campaign['endpoint_Teff_K']:return 'review 1 K endpoint',case
    if not all((evo/f).is_file() for f in ['final.checkpoint','final.checkpoint.lithium']):return 'inspect checkpoint',case
    if reason in ('cpu budget','accepted-step budget','requested age') and r['accepted_this_invocation']>0:
        j=read(case/'launch.json')
        if r['age_years']-j['start_age_years']>=max(1e-4,1e-8*r['age_years']):return 'continue',case
    diagnosis=case/'gravity-review.json'
    if any(x in reason for x in ('atmosphere','pressure family','envelope surface','thermal reservoir')):
        if not diagnosis.exists():return 'assess atmosphere',case
        review=read(diagnosis)
        if review.get('eligible'):return 'gravity control',case
        if review.get('reason')=='one independent atmosphere column needed':
            column=case/'native-atmosphere'
            if not column.exists():return 'compute atmosphere',case
            if locked(column/'worker.lock'):return 'atmosphere running',case
            return 'assess atmosphere',case
    if any(x in reason for x in ('outside','source support','coverage','domain','table')):return 'check table coverage',case
    if any(x in reason for x in ('convergence','line search','repeated step rejection','thermal reservoir')):
        if not read(case/'launch.json')['short_step']:return 'short-step control',case
    return 'review physics or solver',case
def start(m,action,case):
    dest=root/f'm{m:02}'/f's{len(dirs(m)):04}';dest.mkdir(parents=True)
    cfg=(reference/'solar-birth.cfg').read_text() if case is None else (case/'config.txt').read_text()
    guess=None
    if action=='alternate birth guess':
        guess=alternate_birth(m)
        if guess is None:raise ValueError('birth guess unavailable')
        prior=dict(shlex.split(l) for l in guess.read_text().splitlines() if l.strip())
        for key in ('initial_radius_Rsun','initial_Teff_K','initial_entropy_loss','initial_step_years'):
            if key in prior:cfg=re.sub(r'(?m)^'+key+r' ".*"$',f'{key} "{prior[key]}"',cfg)
        case=None
    cfg=re.sub(r'(?m)^mass_Msun ".*"$',f'mass_Msun "{m/100:.3f}"',cfg)
    cfg=re.sub(r'(?m)^zone_threads ".*"$','zone_threads "1"',cfg)
    if action=='gravity control':cfg=(case/'gravity-control/config.txt').read_text()
    cfgpath=dest/'config.txt';cfgpath.write_text(cfg)
    age=0 if case is None else last(case/'evolution/history.jsonl')['years']
    target=max(1e11,2*age)
    cmd=[str(solver),'--lifetime',str(cfgpath),str(target),str(dest/'evolution'),'--cpu-seconds',str(campaign['segment_cpu_seconds'])]
    short=action=='short-step control'
    if case:
        checkpoint=case/('gravity-control/start.checkpoint' if action=='gravity control' else 'evolution/final.checkpoint')
        if action=='resume interrupted checkpoint':checkpoint=case/'evolution/latest.checkpoint'
        cmd+=['--restart',str(checkpoint)]
        if action=='gravity control':cmd+=['--restart-step-years','100000']
        if short:
            prev=last(case/'evolution/history.jsonl')['step_years']
            cmd+=['--restart-step-years',str(max(1,min(prev/4,1e7)))]
    with (dest/'worker.lock').open('a') as lock,(dest/'run.log').open('w') as log:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        identity=dict(mass_Msun=m/100,initial_Z=solar,start_age_years=age,action=action,short_step=short,
                      command=cmd,source=str(case) if case else None,initial_guess_source=str(guess) if guess else None)
        write(dest/'launch.json',identity)
        child=subprocess.Popen(cmd,stdout=log,stderr=log,stdin=subprocess.DEVNULL,pass_fds=(lock.fileno(),),
                               start_new_session=True,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
        write(dest/'launch.json',dict(identity,pid=child.pid))
def watch():
    mix=dict(shlex.split(l) for l in (reference/'solar-birth.cfg').read_text().splitlines() if l.strip())
    values=list(map(float,Path(mix['composition']).read_text().split()))
    assert sorted(masses+[7])==[round(t['mass_Msun']*100) for t in campaign['tracks']]
    assert math.isclose(sum(values[3:8]),solar,abs_tol=1e-15)
    assert read(reference/'reduced-cooling/sequence-verification.json')['initial_Z']==solar
    handle=(root/'coordinator.lock').open('a')
    try:fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:return
    while True:
        states={}
        for m in masses:
            try:states[m]=next_step(m)
            except (OSError,ValueError,KeyError,IndexError) as error:states[m]=(f'inspect case records: {error}',None)
        active=sum(action in ('running','atmosphere running') for action,_ in states.values())
        column_busy=any(action=='atmosphere running' for action,_ in states.values())
        hold=next((s for s in ('STOP','PAUSE_NATIVE') if (root/s).exists()),None)
        if shutil.disk_usage(root).free<campaign['minimum_free_GiB']*1024**3:hold=hold or 'disk reserve'
        if (root/'.supervision/operator.json').exists() and read(root/'.supervision/operator.json').get('action')=='stop':hold='operator stop'
        queue=[m for m,(s,_) in states.items() if s in ('birth','continue','short-step control','assess atmosphere','gravity control','resume interrupted checkpoint','compute atmosphere','alternate birth guess','atmosphere birth control') and (s!='compute atmosphere' or not column_busy)]
        from ember_campaign_recovery import tick
        # Always count adopted/background jobs before allocating the two lanes.
        # Pending evolution takes priority over launching a new research task.
        available=campaign['lanes']-active
        try:
            repairs,recovery=tick(root,states,available,not hold and not queue)
        except (OSError,ValueError,KeyError,RuntimeError,subprocess.SubprocessError) as error:
            repairs=campaign['lanes']
            recovery={'observation_error':str(error)}
        active+=repairs
        queue.sort(key=lambda m:(len(dirs(m)),masses.index(m)))
        for m in queue[:max(0,campaign['lanes']-active)] if not hold else []:
            action,case=states[m]
            if action=='compute atmosphere':
                if column_busy:continue
                from atmosphere_repair import launch_column
                launch_column(case,m,read(case/'gravity-review.json'))
                column_busy=True;states[m]=('atmosphere running',case)
                continue
            if action=='assess atmosphere':
                from atmosphere_repair import prepare,score_column
                try:
                    col=case/'native-atmosphere'
                    if col.exists() and not (col/'assessment.json').exists():
                        write(col/'assessment.json',score_column(col))
                    review=prepare(case,m)
                except (OSError,ValueError,KeyError,AssertionError,subprocess.SubprocessError) as error:
                    review=dict(eligible=False,reason=f'could not validate local overlap: {error}')
                write(case/'gravity-review.json',review)
                continue
            try:
                if action=='atmosphere birth control':
                    from birth_control import run
                    run(m)
                else:start(m,action,case)
                states[m]=('running',dirs(m)[-1])
            except (OSError,ValueError,KeyError) as error:
                states[m]=(f'inspect launch: {error}',case)
        rows={}
        for m,(s,c) in states.items():
            row=dict(status=s)
            if c:
                row['case']=str(c.relative_to(root))
                if (c/'evolution/history.jsonl').exists():
                    v=last(c/'evolution/history.jsonl')
                    row.update({k:v[k] for k in ('years','Teff_K','accepted') if k in v})
                if (c/'evolution/report.json').exists():row['reason']=read(c/'evolution/report.json')['stop_reason']
            rows[f'{m/100:.2f}']=row
        rows['0.07']=dict(status='existing conditional sequence',record=str(reference/'reduced-cooling/sequence-verification.json'))
        science=any(t['status'] in ('running','atmosphere running') for t in rows.values()) or 'running' in recovery.values()
        blocked=not hold and not science and not queue and 'ready' not in recovery.values()
        phase=('paused' if hold else 'observation error' if 'observation_error' in recovery else
               'needs intervention' if blocked else 'working')
        if phase in ('needs intervention','observation error'):
            groups={}
            for mass,r in rows.items():
                groups.setdefault(r['status'],[]).append(mass)
            notice=dict(phase=phase,groups=groups,recovery=recovery,
                        next_action='Inspect failed recovery receipts, add a distinct bounded investigation to recovery-plan.json, or supply a checked restart.')
            path=root/'needs-attention.json'
            if not path.exists() or read(path)!=notice:
                write(path,notice)
                print('Grid needs attention. See '+str(path),flush=True)
                subprocess.Popen(['/usr/bin/osascript','-e','display notification "No further automatic science work is available. See needs-attention.json." with title "Ember solar grid"'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        else:(root/'needs-attention.json').unlink(missing_ok=True)
        write(root/'status.json',dict(pid=os.getpid(),updated_unix=time.time(),initial_Z=solar,hold=hold,phase=phase,recovery=recovery,tracks=rows))
        if hold=='STOP' and not active:break
        time.sleep(15)
if __name__=='__main__':watch()
