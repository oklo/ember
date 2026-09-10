#!/usr/bin/env python3
"""Check independent TOPS compositions with Ember's actual table interpolator.

Temperature and density are sampled from a frozen stellar profile and treated
as source coordinates. There is no isotope remapping, structural re-solve or
estimate of stellar luminosity/lifetime error. AESOPUS is held fixed in its
blend with TOPS. The comparison isolates TOPS composition interpolation.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

from audit_tops_heldout import KEV_TO_K, active_top_weight
from import_tops_composition import read


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def header(path):
    with path.open() as f:
        nx,nt,nr=map(int,f.readline().split()[:3])
        rr=list(map(float,f.readline().split()));tt=list(map(float,f.readline().split()))
        xs=[]
        for _ in range(nx):
            x,z=map(float,f.readline().split());xs.append(x)
            for _ in range(nt):f.readline()
    if len(rr)!=nr or len(tt)!=nt:raise ValueError('invalid runtime axes')
    return xs,tt,rr,z


def probe(executable,low,high,queries):
    p=subprocess.run([str(executable.resolve()),str(low),str(high)],input=queries,
                     text=True,capture_output=True,check=True)
    rows=[list(map(float,line.split())) for line in p.stdout.splitlines()]
    if len(rows)!=len(queries.splitlines()) or any(len(r)!=3 or r[0]<=0 or not all(math.isfinite(v) for v in r) for r in rows):
        raise ValueError('invalid runtime probe output')
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['probe','family_directory','heldout_manifest','profile','work','output']:
        p.add_argument(name,type=Path)
    a=p.parse_args()
    if a.work.exists():raise FileExistsError('use a new diagnostic work directory')
    a.work.mkdir(parents=True)
    track=json.loads(a.profile.read_text())
    profile=[dict(zip(track['profile_columns'],row,strict=True)) for row in track['profile']]
    records=json.loads(a.heldout_manifest.read_text())['planes'];comparisons=[];table_hashes={}
    for record in records:
        x,z=record['X'],record['Z'];prefix=f'tops_gs98_mixture_z{round(z*1000):03d}'
        family=[a.family_directory/(prefix+f'_{label}.dat') for label in ['low','high']]
        axes=[header(f) for f in family]
        if any(x in h[0] or not h[0][0]<x<h[0][-1] or z!=h[3] for h in axes):
            raise ValueError('heldout is not independently bracketed')
        source=a.heldout_manifest.parent/record['file'];request=a.heldout_manifest.parent/record['request']
        if sha(request)!=record['request_sha256']:raise ValueError('heldout request changed')
        ts,rs,cells,excluded=read(source,record)
        reference=[]
        for label,path,(_,temps,densities,_) in zip(['low','high'],family,axes,strict=True):
            tt=[min(ts,key=lambda t:abs(math.log10(t*KEV_TO_K)-q)) for q in temps]
            rr=[min(rs,key=lambda r:abs(math.log10(r)-q)) for q in densities]
            if any(abs(math.log10(t*KEV_TO_K)-q)>1e-12 for t,q in zip(tt,temps)) or any(abs(math.log10(r)-q)>1e-12 for r,q in zip(rr,densities)):
                raise ValueError('runtime axes disagree with original source')
            if any((t,r) in excluded for t in tt for r in rr):raise ValueError('heldout has substituted cells in family rectangle')
            target=a.work/f'x{x:g}_z{z:g}_{label}.dat'
            lines=[f'1 {len(tt)} {len(rr)} independent original source plane',
                   ' '.join(format(q,'.17g') for q in densities),' '.join(format(q,'.17g') for q in temps),f'{x:.17g} {z:.17g}']
            lines+=[' '.join(format(math.log10(cells[t,r]),'.17g') for r in rr) for t in tt]
            target.write_text('\n'.join(lines)+'\n');reference.append(target);table_hashes[str(path)]=sha(path)
        rectangles={label:(10**h[2][0],10**h[2][-1]) for label,h in zip(['low','high'],axes)}
        selected=[]
        for row in profile:
            T,rho=row['temperature_K'],row['density_g_cm3']
            weight=active_top_weight(T/KEV_TO_K,rho,rectangles,[-8,6])
            if weight is not None:selected.append((T,rho,weight))
        if not selected:raise ValueError('no supported profile states')
        query=''.join(f'{x:.17g} {T:.17g} {rho:.17g}\n' for T,rho,w in selected)
        predicted=probe(a.probe,*family,query);truth=probe(a.probe,*reference,query);rows=[]
        for (T,rho,w),b,c in zip(selected,predicted,truth,strict=True):
            d=math.log(b[0]/c[0]);u=(math.log10(T)-4.4)/.1
            dw=6*u*(1-u)/(.1*math.log(10)) if 0<u<1 else 0.
            rows.append({'T_K':T,'rho_g_cm3':rho,'blended_relative_difference':math.expm1(w*d),
                         'delta_dlnk_dlnT':w*(b[1]-c[1])+dw*d,'delta_dlnk_dlnrho':w*(b[2]-c[2])})
        comparisons.append({'X_atomic':x,'Z_atomic':z,'profile_points':len(rows),'heldout':record,
                            'worst':max(rows,key=lambda r:abs(r['blended_relative_difference'])),
                            'max_abs_difference_dlnk_dlnT':max(abs(r['delta_dlnk_dlnT']) for r in rows),
                            'max_abs_difference_dlnk_dlnrho':max(abs(r['delta_dlnk_dlnrho']) for r in rows)})
    report={'scope':__doc__,'probe_sha256':sha(a.probe),'profile_sha256':sha(a.profile),
            'heldout_manifest_sha256':sha(a.heldout_manifest),'runtime_table_sha256':table_hashes,
            'gas_logR_support':[-8,6],'comparisons':comparisons}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps([{'X':r['X_atomic'],'Z':r['Z_atomic'],'worst_percent':100*r['worst']['blended_relative_difference']} for r in comparisons]))


if __name__=='__main__':main()
