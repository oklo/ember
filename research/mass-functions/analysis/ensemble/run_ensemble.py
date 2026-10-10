"""Independent birth realizations; no parameters chosen from brown-dwarf counts.

Reuse the validated private PeTar driver and binary model unchanged. Six jobs,
two OpenMP threads each. Each output retains all surviving escapers.
"""
from pathlib import Path
import argparse, concurrent.futures, hashlib, json, os, subprocess, sys, time

W=Path(__file__).resolve().parent
V1=W.parent/'cluster-dissolution-oct9-v1'
sys.path.insert(0,str(V1))
import run_cluster as driver

def one(case):
    cmd=[sys.executable,str(__file__),'--one',case['name']]
    with (W/(case['name']+'.log')).open('w') as log:
        p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT)
        (W/(case['name']+'.pid')).write_text(str(p.pid))
        return {'name':case['name'],'exit_code':p.wait()}

def cases():
    return [dict(name=f'{cl}_r{rh:g}_assort_s{seed}',cluster=cl,rh=rh,seed=seed,
                 n=3000,end=0,threads=2,output=25 if cl=='pleiades' else 50,
                 binary_fraction=.25,binary_mode='assortative',step=0,
                 theta=.3,stellar_evolution=1,group=512)
            for cl,rh in [('hyades',2.5),('pleiades',3.)] for seed in [41,42,43]]

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--one');args=p.parse_args()
    if args.one:
        driver.W=W
        driver.run(argparse.Namespace(**next(c for c in cases() if c['name']==args.one)))
    else:
        assert not (W/'ensemble_plan.json').exists(), 'Do not relaunch this ensemble.'
        (W/'install').symlink_to(V1/'install',target_is_directory=True)
        plan={'cases':cases(),'started_unix':time.time(),
              'driver_sha256':hashlib.sha256((V1/'run_cluster.py').read_bytes()).hexdigest(),
              'selection':'Sizes chosen to seek present cluster mass/size, not fit BD ratios; seed spread is essential.',
              'physics':'Exactly the v1 similar-mass binary population, no initial gas or mass segregation.',
              'maximum_simultaneous_threads':12}
        (W/'ensemble_plan.json').write_text(json.dumps(plan,indent=2)+'\n')
        start=time.monotonic()
        with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
            result=list(pool.map(one,cases()))
        (W/'ensemble_receipt.json').write_text(json.dumps({'results':result,'wall_seconds':time.monotonic()-start},indent=2)+'\n')
        print(result,flush=True)
