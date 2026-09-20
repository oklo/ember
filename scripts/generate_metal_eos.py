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
    p.add_argument('--electron-integrals',choices=['fitted','numerical'],default='fitted',
                   help='numerical uses FreeEOS option 223: integrated electrons, with radiation added by Ember')
    p.add_argument('--grid-from',type=Path,
                   help='reuse exact temperature/density coordinates from a source specification or manifest')
    p.add_argument('--precision-fallback',type=Path,
                   help='pinned tighter-quadrature probe, used only if the nominal source process fails')
    a=p.parse_args()
    if not 1<=a.jobs<=8 or not .005<=a.step<=.05:raise ValueError('invalid jobs or material step')
    if any(sorted(set(axis))!=axis for axis in [a.hydrogen,a.helium3]):raise ValueError('axes must increase')
    ts=[3.5+i*a.step for i in range(round(3.6/a.step)+1)]
    qs=[-1.5+i*a.step for i in range(round(4/a.step)+1)]
    source_sha=sha(a.probe.read_bytes());a.work.mkdir(parents=True,exist_ok=True)
    grid_reference={}
    if a.grid_from:
        grid=json.loads(a.grid_from.read_text())
        if (grid['probe_sha256']!=source_sha
                or grid['source_archive_sha256']!='4ab1c15a51385a3eab3b08c6f3f240739c0105d92ec828d635ac95720edefb09'):
            raise ValueError('reference grid uses different source physics')
        ts,qs=grid['logT'],grid['logQ']
        for axis in [ts,qs]:
            if (len(axis)<5 or not all(math.isfinite(v) for v in axis)
                    or any(v>=w for v,w in zip(axis,axis[1:]))
                    or any(abs((w-v)/a.step-1)>1e-10 for v,w in zip(axis,axis[1:]))):
                raise ValueError('reference grid does not match the requested uniform material step')
        grid_reference={'grid_reference':str(a.grid_from.resolve()),
                        'grid_reference_sha256':sha(a.grid_from.read_bytes())}
    spec={'hydrogen':a.hydrogen,'helium3':a.helium3,'logT':ts,'logQ':qs,'probe_sha256':source_sha,
          'source_archive_sha256':'4ab1c15a51385a3eab3b08c6f3f240739c0105d92ec828d635ac95720edefb09',**grid_reference}
    options=[3,223,-2] if a.electron_integrals=='numerical' else [3,1,-2]
    if options!=[3,1,-2]:
        spec['source_options']=options
        spec['source_radiation_included']=False
    fallback=None
    if a.precision_fallback:
        fallback=json.loads(a.precision_fallback.read_text())
        if (options!=[3,223,-2] or fallback['nominal_probe_sha256']!=source_sha
                or fallback['source_archive_sha256']!=spec['source_archive_sha256']
                or fallback['nominal_relative_integral_error']!=1e-9
                or fallback['fallback_relative_integral_error']!=1e-11
                or fallback['changed_Fortran_files']!=['src/fermi_dirac_direct.f90']
                or sha(Path(fallback['probe']).read_bytes())!=fallback['probe_sha256']
                or sha(Path(fallback['library']).read_bytes())!=fallback['library_sha256']
                or sha(Path(fallback['source_control_report']).read_bytes())!=fallback['source_control_sha256']):
            raise ValueError('unverified source precision fallback')
        control=json.loads(Path(fallback['source_control_report']).read_text())
        if control['queries']<100 or control['maximum_scaled_physical_difference']>=1e-8:
            raise ValueError('precision fallback lacks matching source controls')
        spec['precision_fallback']={**fallback,'receipt_sha256':sha(a.precision_fallback.read_bytes())}
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
        request=' '.join(map(str,m['eps']))+'\n'+' '.join(map(str,options))+'\n'+''.join(
            f'{math.log(scale)+math.log(10)*(q+1.5*(t-6)):.17g} {math.log(10)*t:.17g}\n' for q in qs)
        fingerprint=sha((source_sha+request).encode())
        if path.exists():
            saved=json.loads(gzip.decompress(path.read_bytes()))
            if saved['input_sha256']!=fingerprint:raise ValueError('cached source input mismatch')
            if saved.get('precision_fallback') and saved['precision_fallback']!=spec.get('precision_fallback'):
                raise ValueError('cached source precision differs')
            if len(saved['data'])!=len(qs) or any(len(r)!=22 or not all(math.isfinite(v) for v in r) for r in saved['data']):
                raise ValueError('invalid cached source responses')
            return
        result=subprocess.run([str(a.probe.resolve())],input=request,text=True,capture_output=True,timeout=600)
        def responses(result):
            try:rows=[list(map(float,line.split())) for line in result.stdout.splitlines()]
            except ValueError:return None
            if result.returncode or len(rows)!=len(qs) or any(len(r)!=22 or not all(math.isfinite(v) for v in r) for r in rows):return None
            return rows
        rows=responses(result);override={}
        if rows is None and fallback:
            nominal=d/f'temperature-{it:03d}.nominal-failure.log'
            nominal.write_text(result.stdout+'\n'+result.stderr)
            result=subprocess.run([fallback['probe']],input=request,text=True,capture_output=True,timeout=600)
            rows=responses(result)
            override={'precision_fallback':spec['precision_fallback'],
                      'actual_probe_input_sha256':sha((fallback['probe_sha256']+request).encode()),
                      'nominal_failure_sha256':sha(nominal.read_bytes())}
        if rows is None:
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
        record={'input_sha256':fingerprint,'input':request,'stderr':result.stderr,'data':rows,**override}
        path.write_bytes(gzip.compress((json.dumps(record,separators=(',',':'),allow_nan=False)+'\n').encode(),mtime=0))
        print(f"XH={m['hydrogen']:g} X3={m['helium3']:g} logT={t:g} complete; {failed} source failures retained for masking",flush=True)

    jobs=[(d,m,it,t) for it,t in enumerate(ts) for d,m in planes]
    with ThreadPoolExecutor(a.jobs) as pool:list(pool.map(run,jobs))
    for d,m in planes:
        data=[];overrides=[]
        for it in range(len(ts)):
            cache=json.loads(gzip.decompress((d/f'temperature-{it:03d}.json.gz').read_bytes()))
            data.extend(cache['data'])
            if cache.get('precision_fallback'):
                overrides.append({'temperature_index':it,'logT':ts[it],
                                  **{key:cache[key] for key in ['actual_probe_input_sha256','nominal_failure_sha256']}})
        raw={**m,'version':'FreeEOS 3.0.0','options':options,'logT':ts,'logQ':qs,
             'source_archive_sha256':spec['source_archive_sha256'],'probe_sha256':source_sha,'data':data}
        if fallback:
            raw['precision_fallback']=spec['precision_fallback']
            raw['precision_fallback_isotherms']=overrides
        (d/'source.json.gz').write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'),allow_nan=False)+'\n').encode(),mtime=0))
        print('assembled',d,flush=True)
    if fallback and (sha(Path(fallback['probe']).read_bytes())!=fallback['probe_sha256']
                     or sha(Path(fallback['library']).read_bytes())!=fallback['library_sha256']):
        raise ValueError('precision fallback source changed during calculation')


if __name__=='__main__':main()
