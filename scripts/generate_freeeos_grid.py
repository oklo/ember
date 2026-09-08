#!/usr/bin/env python3
"""Generate a raw fixed-X H/He EOS1 grid using an external FreeEOS 3.0 probe.

Build freeeos_probe.f90 as documented in docs/FREEEOS.md. The probe is an
offline source-data generator, not a runtime dependency of ember.
"""
import argparse
import concurrent.futures
import gzip
import json
import math
from pathlib import Path
import subprocess


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('probe',type=Path)
    ap.add_argument('output',type=Path)
    ap.add_argument('--step',type=float,default=.0125)
    ap.add_argument('--hydrogen',type=float,default=.7)
    a=ap.parse_args()
    if not math.isfinite(a.step) or a.step<.005 or a.step>.1:
        raise ValueError('step must be .005 .. .1 dex')
    ts=[3.2+i*a.step for i in range(round(3.8/a.step)+1)]
    qs=[-5+i*a.step for i in range(round(7.4/a.step)+1)]
    if not math.isfinite(a.hydrogen) or not 0 <= a.hydrogen <= .98:
        raise ValueError('hydrogen must be 0 .. .98')
    # Preserve the original .7 reference's exact floating-point inputs.
    helium = .3 if a.hydrogen == .7 else 1-a.hydrogen
    eps=[a.hydrogen/1.00782503,helium/4.00260325]+[0.]*18
    header=' '.join(map(str,eps))+'\n3 1 -2\n'
    def run(part):
        text=header+''.join(f'{math.log(10)*(q+1.5*(t-6)):.17g} {math.log(10)*t:.17g}\n'
                            for t in part for q in qs)
        p=subprocess.run([str(a.probe.resolve())],input=text,text=True,capture_output=True,timeout=600)
        rows=[[float(x) for x in line.split()] for line in p.stdout.splitlines()]
        if p.returncode or len(rows)!=len(part)*len(qs):
            raise RuntimeError((p.returncode,len(rows),p.stderr[-1000:]))
        if any(len(r)!=22 or r[0]!=0 or not all(math.isfinite(x) for x in r) for r in rows):
            raise ValueError('failed or invalid direct FreeEOS source state: '+repr([(part[k//len(qs)],qs[k%len(qs)],r[:6]) for k,r in enumerate(rows) if len(r)!=22 or r[0]!=0 or not all(math.isfinite(x) for x in r)]))
        return rows
    # A fresh process per isotherm avoids carrying the source library's
    # cached ionization state across the large dense-to-dilute row reset.
    # Within a row queries increase smoothly in density.
    parts=[[t] for t in ts]
    rows=[]
    with concurrent.futures.ThreadPoolExecutor(4) as pool:
        for part,result in zip(parts,pool.map(run,parts)):
            rows.extend(result)
            print(f'completed logT={part[-1]:.4f}: {len(rows)} nodes',flush=True)
    obj=dict(version='FreeEOS 3.0.0',options=[3,1,-2],eps=eps,logT=ts,logQ=qs,
             source_archive_sha256='4ab1c15a51385a3eab3b08c6f3f240739c0105d92ec828d635ac95720edefb09',
             execution='Fresh probe process per isotherm; increasing density; no fallback',
             columns=['info','iterations','rho','T','P','E','S','chiRho','chiT','E_r','E_t','S_r','S_t',
                      'cp','grad_ad','delta','mu_inv','free_e','H2','Hplus','Heplus','Heplusplus'],data=rows)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('wb') as f:
        with gzip.GzipFile(fileobj=f,filename='',mode='wb',mtime=0) as g:
            g.write((json.dumps(obj,separators=(',',':'),allow_nan=False)+'\n').encode())


if __name__=='__main__':
    main()
