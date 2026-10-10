"""Bounded independent field runs, at most four two-thread jobs together."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess,time,json,hashlib,datetime,sys
W=Path(__file__).resolve().parent
BIN=W.parents[1]/'out/luhman16-independent-oct8-v1/build/apps/ember-evolve'
def run(case):
    d=W/'tracks'/case;cfg=d/'configuration.txt'
    if (d/'launch.json').exists():raise RuntimeError(f'Refusing to overwrite {case}')
    scope=json.loads((d/'scope.json').read_text());binary=Path(scope.get('executable',BIN))
    cmd=[str(binary),'--lifetime',str(cfg),str(scope['target_years']),str(d/'evolution'),
         '--max-steps',str(scope['maximum_steps']),'--cpu-seconds',str(scope['maximum_cpu_seconds'])]
    rec=dict(command=cmd,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
         binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
         configuration_sha256=hashlib.sha256(cfg.read_bytes()).hexdigest())
    start=time.monotonic()
    with (d/'stdout.log').open('w') as out,(d/'stderr.log').open('w') as err:
        p=subprocess.Popen(cmd,stdout=out,stderr=err)
        rec['pid']=p.pid;(d/'launch.json').write_text(json.dumps(rec,indent=2)+'\n')
        code=p.wait()
    rec.update(returncode=code,wall_seconds=time.monotonic()-start,
               ended_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    report=d/'evolution/report.json'
    if report.exists():rec['report']=json.loads(report.read_text())
    (d/'completion.json').write_text(json.dumps(rec,indent=2)+'\n')
    print(case,code,round(rec['wall_seconds'],1),flush=True)
if __name__=='__main__':
    cases=sys.argv[1:] or json.loads((W/'cases.json').read_text())
    with ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(run,cases))
