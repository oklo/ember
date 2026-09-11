#!/usr/bin/env python3
"""Compare each saved stellar zone with direct FreeEOS at its actual mixture.

These checks concern the selected profile and interpolation of the adopted EOS.
They do not bound the physical uncertainty of FreeEOS or unvisited mixtures.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
import math
from pathlib import Path
import shlex
import subprocess

from metal_eos_composition import mixture


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['track','family','ember_probe','source_probe','output']:
        p.add_argument(name,type=Path)
    p.add_argument('--extra-queries',type=Path)
    p.add_argument('--jobs',type=int,default=4)
    a=p.parse_args()
    if not 1<=a.jobs<=8:raise ValueError('invalid job count')
    track=json.loads(a.track.read_text());columns=track['profile_columns']
    profile=[dict(zip(columns,row,strict=True)) for row in track['profile']]
    queries=[[r['X'],r['Y3'],r['temperature_K'],r['density_g_cm3']] for r in profile]
    extra=json.loads(a.extra_queries.read_text()) if a.extra_queries else []
    queries+=extra
    inputs={str(f.resolve()):digest(f) for f in [a.track,a.family,a.ember_probe,a.source_probe,Path(__file__)]}
    if a.extra_queries:inputs[str(a.extra_queries.resolve())]=digest(a.extra_queries)
    for line in a.family.read_text().splitlines()[3:]:
        f=a.family.parent/shlex.split(line)[0];inputs[str(f.resolve())]=digest(f)
    groups={}
    for i,(x,y,T,rho) in enumerate(queries):groups.setdefault((x,y),[]).append((i,T,rho))
    def source(group):
        (x,y),points=group;m=mixture(x,y);scale=m['source_mass_scale']
        request=' '.join(map(str,m['eps']))+'\n3 1 -2\n'+''.join(
            f'{math.log(scale*rho):.17g} {math.log(T):.17g}\n' for _,T,rho in points)
        r=subprocess.run([str(a.source_probe.resolve())],input=request,text=True,capture_output=True,check=True)
        rows=[list(map(float,line.split())) for line in r.stdout.splitlines()]
        if len(rows)!=len(points):raise ValueError('incomplete direct source result')
        result=[]
        for (i,T,rho),row in zip(points,rows,strict=True):
            if (len(row)!=22 or not all(math.isfinite(v) for v in row) or row[0]!=0
                    or abs(row[2]/(scale*rho)-1)>1e-9 or abs(row[3]/T-1)>1e-10):
                raise ValueError('direct source did not converge at requested profile state')
            result.append((i,{'P':row[4],'E':row[5]*scale,'cv':row[12]*scale,'cp':row[13]*scale,'grad_ad':row[14]}))
        return result,{'mixture':[x,y],'input':request,'stdout':r.stdout,'stderr':r.stderr}
    with ThreadPoolExecutor(a.jobs) as pool:calculated=list(pool.map(source,groups.items()))
    direct={i:r for rows,_ in calculated for i,r in rows}
    raw=a.output.with_suffix('.source.json.gz')
    raw.write_bytes(gzip.compress(json.dumps([r for _,r in calculated],allow_nan=False).encode(),mtime=0))
    request=''.join(' '.join(format(v,'.17g') for v in q)+'\n' for q in queries)
    response=subprocess.run([str(a.ember_probe.resolve()),str(a.family.resolve())],
                            input=request,text=True,capture_output=True,check=True)
    records=[]
    for i,(q,line) in enumerate(zip(queries,response.stdout.splitlines(),strict=True)):
        row=list(map(float,line.split()))
        if len(row)!=18 or row[:4]!=q or not all(math.isfinite(v) for v in row):
            raise ValueError('unexpected Ember EOS response')
        differences={k:row[j]/direct[i][k]-1 for k,j in [('P',4),('E',5),('cv',10),('cp',9),('grad_ad',11)]}
        records.append({'query':q,'profile_zone':i if i<len(profile) else None,
                        'direct':direct[i],'ember':row,'relative_difference':differences})
    limits={'P':.001,'E':.002,'cv':.003,'cp':.003,'grad_ad':.003}
    maximum={k:max(abs(r['relative_difference'][k]) for r in records[:len(profile)]) for k in limits}
    passed=all(maximum[k]<limits[k] for k in limits)
    if any(digest(f)!=expected for f,expected in inputs.items()):raise ValueError('input changed during audit')
    report={'scope':__doc__,'profile_passed':passed,'profile_zones':len(profile),'diagnostic_queries':len(extra),
            'maximum_profile_relative_difference':maximum,'source_comparison_limits':limits,
            'raw_source_sha256':digest(raw),'input_sha256':inputs,'records':records}
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['records','input_sha256']}))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
