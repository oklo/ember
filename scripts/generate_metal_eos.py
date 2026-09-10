#!/usr/bin/env python3
"""Generate restartable metal-bearing baryonic FreeEOS material source planes."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess
from metal_eos_composition import mixture


def sha(data):return hashlib.sha256(data).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('probe',type=Path);p.add_argument('work',type=Path)
    p.add_argument('--hydrogen',type=float,nargs='+',default=[.3,.4,.5,.6,.7,.75])
    p.add_argument('--helium3',type=float,nargs='+',default=[0,.12])
    p.add_argument('--step',type=float,default=.0125)
    p.add_argument('--jobs',type=int,default=4)
    a=p.parse_args()
    if not 1<=a.jobs<=8 or not .005<=a.step<=.05:raise ValueError('invalid jobs or material step')
    if any(sorted(set(axis))!=axis for axis in [a.hydrogen,a.helium3]):raise ValueError('axes must increase')
    ts=[3.5+i*a.step for i in range(round(3.6/a.step)+1)]
    qs=[-1.5+i*a.step for i in range(round(4/a.step)+1)]
    source_sha=sha(a.probe.read_bytes());a.work.mkdir(parents=True,exist_ok=True)
    spec={'hydrogen':a.hydrogen,'helium3':a.helium3,'logT':ts,'logQ':qs,'probe_sha256':source_sha,
          'source_archive_sha256':'4ab1c15a51385a3eab3b08c6f3f240739c0105d92ec828d635ac95720edefb09'}
    manifest=a.work/'specification.json'
    if manifest.exists() and json.loads(manifest.read_text())!=spec:raise ValueError('changed source/settings require a new work directory')
    manifest.write_text(json.dumps(spec,indent=2)+'\n')
    planes=[]
    for i,(x,y) in enumerate(itertools.product(a.hydrogen,a.helium3)):
        d=a.work/f'plane-{i:03d}';d.mkdir(exist_ok=True)
        m=mixture(x,y);(d/'mixture.json').write_text(json.dumps(m,indent=2)+'\n');planes.append((d,m))

    def run(job):
        d,m,it,t=job;path=d/f'temperature-{it:03d}.json.gz'
        scale=m['source_mass_scale']
        request=' '.join(map(str,m['eps']))+'\n3 1 -2\n'+''.join(
            f'{math.log(scale)+math.log(10)*(q+1.5*(t-6)):.17g} {math.log(10)*t:.17g}\n' for q in qs)
        fingerprint=sha((source_sha+request).encode())
        if path.exists():
            saved=json.loads(gzip.decompress(path.read_bytes()))
            if saved['input_sha256']!=fingerprint:raise ValueError('cached source input mismatch')
            return
        result=subprocess.run([str(a.probe.resolve())],input=request,text=True,capture_output=True,timeout=600)
        rows=[list(map(float,line.split())) for line in result.stdout.splitlines()]
        if result.returncode or len(rows)!=len(qs) or any(len(r)!=22 or not all(math.isfinite(v) for v in r) for r in rows):
            (d/f'temperature-{it:03d}.failure.log').write_text(result.stdout+'\n'+result.stderr)
            raise ValueError(f'{d.name} logT={t}: failed source state')
        # Nonconverged states remain explicitly flagged in the raw source.
        # The potential importer excludes their entire derivative stencil;
        # no invented thermodynamic values can enter a supported cell.
        failed=sum(r[0]!=0 for r in rows)
        # Material and radiation energies per atomic source gram both acquire
        # the same scale; P and dimensionless responses do not change.
        for r in rows:
            r[2]/=scale
            for k in [5,6,9,10,11,12,13,16,17]:r[k]*=scale
        record={'input_sha256':fingerprint,'input':request,'stderr':result.stderr,'data':rows}
        path.write_bytes(gzip.compress((json.dumps(record,separators=(',',':'),allow_nan=False)+'\n').encode(),mtime=0))
        print(f"XH={m['hydrogen']:g} X3={m['helium3']:g} logT={t:g} complete; {failed} source failures retained for masking",flush=True)

    jobs=[(d,m,it,t) for it,t in enumerate(ts) for d,m in planes]
    with ThreadPoolExecutor(a.jobs) as pool:list(pool.map(run,jobs))
    for d,m in planes:
        data=[]
        for it in range(len(ts)):data.extend(json.loads(gzip.decompress((d/f'temperature-{it:03d}.json.gz').read_bytes()))['data'])
        raw={**m,'version':'FreeEOS 3.0.0','options':[3,1,-2],'logT':ts,'logQ':qs,
             'source_archive_sha256':spec['source_archive_sha256'],'probe_sha256':source_sha,'data':data}
        (d/'source.json.gz').write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'),allow_nan=False)+'\n').encode(),mtime=0))
        print('assembled',d,flush=True)


if __name__=='__main__':main()
