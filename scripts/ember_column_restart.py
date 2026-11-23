"""Decide whether an interrupted native column has a new, usable numerical seed."""
import hashlib
import json
import math
import os
from pathlib import Path
import shutil

def read(path):
    return json.loads(path.read_text())

def seed(directory):
    spec=read(directory/'spec.json')
    suffix=spec['anchor']+spec['tag']+'.json'
    if (directory/('hhe_arm_'+suffix)).exists():return None
    path=directory/('latest_guess_'+suffix)
    start=directory/('start_'+suffix)
    if not path.is_file() or not start.is_file():return None
    raw=path.read_bytes()
    g=json.loads(raw);initial=read(start)
    if g.get('accepted') is not False:return None
    n=spec['levels'];p=g['pressure_bar'];t=g['temperature_K']
    if len(p)!=n or len(t)!=n or len(initial['temperature_K'])!=n:return None
    if any(not math.isfinite(x) or x<=0 for x in p+t):return None
    if any(b<=a for a,b in zip(p,p[1:])) or p!=initial['pressure_bar']:return None
    if t==initial['temperature_K'] or g.get('outer_calls',0)<1:return None
    top=initial['rcb_guess']
    if not isinstance(top,int) or not 3<=top<n-1:return None
    return dict(path=str(path),sha256=hashlib.sha256(raw).hexdigest(),
                profile_sha256=hashlib.sha256(json.dumps([p,t],separators=(',',':')).encode()).hexdigest(),top=top)

def retry_job(root,job,parent,checkpoint,number):
    dst=root/'recovery'/job['name']/f'column-{number:02}'
    cfg=job['resume_column']
    spec=read(parent/'spec.json')
    receipt=dst/('hhe_arm_'+spec['anchor']+spec['tag']+'.json')
    return dict(name=job['name']+f'-resume-{number:02}',column_directory=str(dst),
                command=[cfg['python'],str(dst/'solve.py'),spec['anchor'],'spec.json',str(cfg['cpu_seconds'])],
                cwd=str(dst),wall_seconds=job.get('wall_seconds',cfg['cpu_seconds']*2),
                outputs=[str(receipt)],seed=checkpoint,parent=str(parent))

def prepare(job):
    dest=Path(job['column_directory']);parent=Path(job['parent'])
    dest.mkdir(parents=True,exist_ok=True);(dest/'grid_hhe').mkdir(exist_ok=True)
    raw=Path(job['seed']['path']).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=job['seed']['sha256']:raise ValueError('checkpoint changed before dispatch')
    frozen=dest/'resume-guess.json';frozen.write_bytes(raw)
    spec=read(parent/'spec.json')
    spec.update(resume_guess=str(frozen),resume_region_top=job['seed']['top'])
    (dest/'spec.json').write_text(json.dumps(spec,indent=2))
    shutil.copyfile(parent/'solve.py',dest/'solve.py')
    if (parent/'prediction.json').is_file():shutil.copyfile(parent/'prediction.json',dest/'prediction.json')
