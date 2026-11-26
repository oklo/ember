"""Refresh the local release board; optionally request one bounded Luna scheduling decision."""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

from release import board, read, save


BRIEF_VERSION = 4
MODEL = 'gpt-6-luna'


def eligible(work):
    return [key for key, row in work.items() if row['disposition'] == 'open' and row['ready']]


def supervise(root, config, snapshot):
    out = root / 'shepherd'
    work = snapshot['work_orders']
    choices = eligible(work)
    if not current(root, snapshot):
        return 'held, stale or complete'
    if not choices:
        choices = [k for k,o in work.items() if o['disposition'] == 'open']
    if not choices:
        return 'no open decision'
    # Distinguish a failed job from unstarted descendants. Send the frontier
    # and its evidence locations instead of hundreds of repeated graph edges.
    summary = {}
    for key, row in work.items():
        blocked = row['failed_prerequisites']
        summary[key] = dict(question=row['question'], ready=row['ready'], ready_columns=row['ready_columns'],
                            running=row['running'], criterion=row['criterion'], next_action=row['next_action'],
                            target_states=row['target_states'],
                            target_dependencies={n:row['unfinished_dependencies'].get(n,[]) for n in row['targets']},
                            failed_jobs={n:row['owners'][n] for n in blocked[:6]}, total_failed_jobs=len(blocked))
    actions = {key+':'+job:(key,job) for key in choices for job in work[key]['ready']}
    if not actions:
        actions = {key+':review':(key,None) for key in choices}
    payload = dict(brief_version=BRIEF_VERSION, lanes=snapshot['lanes'], occupied_native=snapshot['occupied_native'],
                   occupied_recovery=snapshot['occupied_recovery'], occupied_columns=snapshot['occupied_columns'],
                   actions=list(actions), work=summary, complete=snapshot['complete'])
    signature = hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()
    path = out / 'luna-calls.json'
    calls = read(path) if path.exists() else []
    now = time.time()
    recent = [c for c in calls if c['at'] > now - 86400]
    if len(recent) >= 2 or (recent and recent[-1].get('brief_version') == BRIEF_VERSION and now - recent[-1]['at'] < 3600):
        return 'within call budget cooldown'
    if calls and calls[-1].get('completed') and calls[-1]['signature'] == signature:
        return 'same dispatch evidence'
    call = dict(at=now, model=MODEL, signature=signature, brief_version=BRIEF_VERSION)
    calls.append(call)
    save(path, calls)  # reserve the spend even if killed or timed out
    schema = dict(type='object', properties=dict(action=dict(type='string',enum=list(actions)),
        reason=dict(type='string',maxLength=600), flaw=dict(type='string',maxLength=600)),
        required=['action','reason','flaw'], additionalProperties=False)
    prompt = ('Select one exact listed action for the NEXT free lane. A :review action means staff adjudication. A failed job is terminal for that attempt; dependency failed means a descendant could not run. Prefer closing the '
              '0.06 review or a shared table edge, while advancing independent work when a column '
              'is already running. Your choice affects only the next free slot. The dispatcher prefers a ready non-column when a column is running. Each target has its own prerequisites. Passed jobs are omitted from unfinished_dependencies. Targets listed together are jointly reviewed. When no order has ready jobs, nominate one blocked decision for investigation. Some targets may have failed dependencies. Identify a concrete '
              'process flaw. You may select a priority only. Campaign data follows.\n'+json.dumps(payload))
    binary = config.get('binary') or shutil.which('codex')
    try:
        with tempfile.TemporaryDirectory() as tmp:
            cwd = Path(tmp)
            schema_path = cwd/'response-schema.json'
            answer = cwd/'answer.json'
            save(schema_path,schema)
            cmd = [binary,'exec','--model',MODEL,'--sandbox','read-only','--ignore-user-config',
                   '--ephemeral','--skip-git-repo-check','--disable','plugins','--disable','multi_agent',
                   '-c','model_reasoning_effort="low"','--output-schema',str(schema_path),
                   '--output-last-message',str(answer),'--json','-']
            proc = subprocess.run(cmd, input='Use only the supplied queue data. No tool calls. '+prompt,
                                  capture_output=True,text=True,cwd=cwd,timeout=120,check=True)
            vote = read(answer)
            if vote['action'] not in actions or not all(isinstance(vote.get(k),str) for k in schema['required']):
                raise ValueError('invalid agent decision')
            events = [json.loads(line) for line in proc.stdout.splitlines() if line.startswith('{')]
            usage = next((e.get('usage') for e in reversed(events) if e.get('usage')),None)
        order, job = actions[vote['action']]
        call.update(completed=True, usage=usage, order=order, job=job, **vote)
        if current(root, snapshot):
            save(out/'recommendation.json', dict(call, release_sha256=snapshot['release_sha256'], expires=now+3600))
    except (ValueError, TypeError, OSError, KeyError, subprocess.SubprocessError) as exc:
        call['error'] = type(exc).__name__ + ': ' + str(getattr(exc,'stderr','') or str(exc))[:400]
    save(path, calls)
    return call


def current(root, snapshot):
    status = read(root/'status.json')
    return (not status['hold'] and status['phase'] in ('working','needs intervention') and snapshot['complete'] < 30
            and 0 <= time.time()-status['updated_unix'] < 120
            and not ((root/'.supervision/operator.json').exists() and read(root/'.supervision/operator.json').get('action')=='stop')
            and not any((root/name).exists() for name in ('STOP','PAUSE_NATIVE')))


def run(root):
    out = root/'shepherd'
    out.mkdir(exist_ok=True)
    with (out/'lock').open('a') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return
        status = read(root/'status.json')
        if time.time()-status['updated_unix'] >= 120:
            save(out/'health.json',dict(phase='stale scheduler heartbeat', checked=time.time(),
                                        heartbeat=status['updated_unix']))
            return
        try:
            snapshot = board(root, status)
        except (OSError, ValueError, KeyError, IndexError, TypeError) as exc:
            save(out/'health.json',dict(phase='board unavailable', checked=time.time(), error=str(exc)))
            raise
        save(out/'health.json',dict(phase='observed',checked=time.time(),scheduler_pid=status['pid']))
        save(root/'release-status.json',snapshot)
        config = root/'shepherd-config.json'
        if config.exists() and read(config).get('enabled') is True:
            print(supervise(root,read(config),snapshot),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign',type=Path)
    run(parser.parse_args().campaign.resolve())
