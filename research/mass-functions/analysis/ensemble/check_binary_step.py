"""Compare two unsoftened force-split steps from one 100-Myr binary snapshot."""
from pathlib import Path
import json,os,subprocess,shutil,time,resource
import numpy as np
import analysis_functions as a
W=Path(__file__).resolve().parent;src=W/'runs/hyades_r2.5_assort_s41'
base=W/'step_validation';base.mkdir(exist_ok=True)
par={x.split()[1]:x.split()[2] for x in (src/'input.par').read_text().splitlines()}
meta=json.loads((src/'initial.json').read_text());rows=[]
for name,step in [('control',float(par['s'])),('double',2*float(par['s']))]:
    d=base/name;d.mkdir(exist_ok=True);snap=src/'data.2'
    with snap.open() as stream:assert float(stream.readline().split()[2])==100
    existing=[f for f in d.glob('probe.*') if f.name[6:].isdigit()]
    if existing:
        f=max(existing,key=lambda f:int(f.suffix[1:]));h,p=a.read_snapshot(f)
        if h.time==105 and len(p.mass)==h.n:
            row,_,_,_=a.structure(h,p,meta)
            record=dict(step=step,exit_code=0,wall_seconds=None,cpu_seconds=None,**row)
            (d/'receipt.json').write_text(json.dumps(record,indent=2,default=lambda x:x.item()))
            rows.append(record);continue
    for f in list(src.glob('input.par*'))+list(src.glob('data.2*')):shutil.copy2(f,d/f.name)
    scale=(step/float(par['s']))**(2/3)
    cmd=[str(W/'install/bin/petar'),'-p','input.par','-a','0','-f','probe','-r',str(float(par['r'])*scale),'--r-search-min',str(float(par['r-search-min'])*scale),'--r-bin','0','-s',str(step),'-t','105','-o','5','-w','1','data.2']
    (d/'command.json').write_text(json.dumps({'command':cmd,'threads':2},indent=2))
    env=os.environ.copy();env.update(OMP_NUM_THREADS='2',OMP_STACKSIZE='128M')
    before=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.monotonic()
    with (d/'run.log').open('w') as f:
        p=subprocess.Popen(cmd,cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT);(d/'pid').write_text(str(p.pid));code=p.wait()
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    record=dict(step=step,exit_code=code,wall_seconds=time.monotonic()-start,cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime)
    files=[f for f in d.glob('probe.*') if f.name[6:].isdigit()]
    f=max(files,key=lambda f:int(f.suffix[1:]));h,p=a.read_snapshot(f)
    row,_,_,_=a.structure(h,p,meta);record.update(row)
    (d/'receipt.json').write_text(json.dumps(record,indent=2,default=lambda x:x.item()));rows.append(record)
    print(name,record['wall_seconds'],row['mass_within_jacobi_Msun'],row['rhalf_within_jacobi_pc'],flush=True)
delta={k:rows[1][k]/rows[0][k]-1 for k in ['mass_within_jacobi_Msun','rhalf_within_jacobi_pc','r10_ratio','jacobi_ratio']}
(base/'comparison.json').write_text(json.dumps({'runs':rows,'fractional_differences':delta},indent=2,default=lambda x:x.item())+'\n')
print(delta,flush=True)
