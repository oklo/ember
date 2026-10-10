"""Explicit lower-field normalization, not an IMF value fitted to clusters.

R=.55 is the lower edge of the independent field sensitivity interval. The
existing 4-pc Pleiades configuration was closest in absolute reference-star
counts. Three seeds quantify draw sensitivity; Hyades remains at R=.63.
"""
from pathlib import Path
import argparse,concurrent.futures,json,subprocess,sys,time
import numpy as np
W=Path(__file__).resolve().parent;V1=W.parent/'cluster-dissolution-oct9-v1';sys.path.insert(0,str(V1))
import run_cluster as driver

def low_field_masses(n,rng):
    bdamp=.55*driver.integral(.1,.2,1.3)/driver.integral(.03,.08,.4)
    seg=[(.01,.08,.4,bdamp),(.08,.5,1.3,1.),(.5,50.,2.3,.5)]
    prob=np.array([c*driver.integral(a,b,s) for a,b,s,c in seg]);prob/=prob.sum()
    k=rng.choice(3,n,p=prob);m=np.empty(n)
    for i,(a,b,s,c) in enumerate(seg):
        z=rng.uniform(size=(k==i).sum());m[k==i]=(a**(1-s)+z*(b**(1-s)-a**(1-s)))**(1/(1-s))
    return m

def one(seed):
    cmd=[sys.executable,str(__file__),'--seed',str(seed)]
    with (W/f'field_low_s{seed}.log').open('w') as f:
        p=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT)
        (W/f'field_low_s{seed}.pid').write_text(str(p.pid));return dict(seed=seed,exit_code=p.wait())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--seed',type=int);args=p.parse_args()
    if args.seed:
        driver.W=W;driver.masses=low_field_masses
        driver.run(argparse.Namespace(name=f'pleiades_r4_fieldlow_s{args.seed}',cluster='pleiades',rh=4.,seed=args.seed,n=3000,end=0,threads=2,output=25,binary_fraction=.25,binary_mode='assortative',step=0,theta=.3,stellar_evolution=1,group=512,expected_birth_ratio=.55))
    else:
        assert not (W/'field_sensitivity_plan.json').exists()
        (W/'field_sensitivity_plan.json').write_text(json.dumps(dict(expected_birth_ratio=.55,seeds=[41,42,43],rh_pc=4,started_unix=time.time(),reason=__doc__),indent=2))
        start=time.monotonic()
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:r=list(pool.map(one,[41,42,43]))
        (W/'field_sensitivity_receipt.json').write_text(json.dumps(dict(runs=r,wall_seconds=time.monotonic()-start),indent=2));print(r,flush=True)
