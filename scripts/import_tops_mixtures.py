#!/usr/bin/env python3
"""Import verified TOPS (X,Z) planes, retaining only common original cells."""
import argparse
import hashlib
import json
import math
from pathlib import Path
from import_tops_composition import read


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('manifest',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();spec=json.loads(a.manifest.read_text());a.output.mkdir(parents=True,exist_ok=True)
    planes=spec['planes'];zs=sorted({v['Z'] for v in planes})
    for v in planes:
        request=a.manifest.parent/v['request']
        if hashlib.sha256(request.read_bytes()).hexdigest()!=v['request_sha256']:
            raise ValueError('source request checksum mismatch')
    if len(zs)<2:raise ValueError('at least two metallicities required')
    manifests={label:[f'EMBER_OPACITY_MIXTURE 1 {len(zs)} logRho TOPS_ATOMIC_GS98_{label}'] for label in ['low','high']}
    for z in zs:
        subset=sorted([v for v in planes if v['Z']==z],key=lambda v:v['X']);xs=[v['X'] for v in subset]
        if len(xs)<2 or len(set(xs))!=len(xs):raise ValueError('invalid composition grid')
        sources=[read(a.manifest.parent/v['file'],v) for v in subset]
        tt,rr=sources[0][:2]
        if any(t!=tt or r!=rr for t,r,_,_ in sources):raise ValueError('unaligned source grids')
        for label,temps in [('low',tt),('high',[t for t in tt if t>=.025])]:
            densities=[]
            for r in rr:
                if any((t,r) in excluded for _,_,_,excluded in sources for t in temps):break
                densities.append(r)
            if len(densities)<4:raise ValueError('insufficient common density support')
            lines=[f'{len(xs)} {len(temps)} {len(densities)} TOPS ATOMIC GS98 X/Z family; native density; all substituted cells excluded',
                ' '.join(f'{math.log10(r):.17g}' for r in densities),
                ' '.join(f'{math.log10(t*1e3*1.602176634e-12/1.380649e-16):.17g}' for t in temps)]
            for x,(_,_,cells,_) in zip(xs,sources,strict=True):
                lines.append(f'{x:.17g} {z:.17g}')
                lines.extend(' '.join(f'{math.log10(cells[t,r]):.17g}' for r in densities) for t in temps)
            name=f'tops_gs98_mixture_z{round(z*1000):03d}_{label}.dat'
            (a.output/name).write_text('\n'.join(lines)+'\n');manifests[label].append(f'{z:.17g} "{name}"')
            print(f'Z={z:g} {label}: {len(xs)} X, {len(temps)} T, {len(densities)} rho')
    for label,lines in manifests.items():(a.output/f'tops_gs98_mixture_{label}.dat').write_text('\n'.join(lines)+'\n')
if __name__=='__main__':main()
