"""PeTar's documented timestep/changeover-radius tuning, private restart tests."""
from pathlib import Path
import subprocess,shutil,os,json,time,resource
W=Path(__file__).resolve().parent;src=W/'runs/hyades_r2_s41';base=W/'step_tuning_core';base.mkdir(exist_ok=True)
rows=[]
for k,step in enumerate([.015625,.03125]):
    d=base/f'probe{k}';d.mkdir(exist_ok=True)
    for f in list(src.glob('input.par*'))+list(src.glob('data.4*')):shutil.copy2(f,d/f.name)
    cmd=[str(W/'install/bin/petar'),'-p','input.par','-a','0','-f','probe','-r',str(.03699029302058526*(step/.00390625)**(2/3)),'--r-search-min',str(.05*(step/.00390625)**(2/3)),'--r-bin','0','-s',str(step),'-t',str(125),'-o','5','-w','1','data.4']
    env=os.environ.copy();env.update(OMP_NUM_THREADS='2',OMP_STACKSIZE='128M')
    t=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
    with (d/'run.log').open('w') as f:
        try:res=subprocess.run(cmd,cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=240);code=res.returncode
        except subprocess.TimeoutExpired:code=-999
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    rows.append(dict(step=step,exit_code=code,wall=time.monotonic()-t,cpu=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime));print(rows[-1],flush=True)
(base/'receipts.json').write_text(json.dumps(rows,indent=2))
