#!/usr/bin/env python3
"""Import an (X,Z) family of untouched AESOPUS 2.1 GS98 source cells."""
import argparse
import hashlib
import math
from pathlib import Path
import re
from zipfile import ZipFile
from import_aesopus import SOURCE_SHA256
SUBSET_SHA256="0bfad00740be06cfa13b44d4375b9a893876d9196ed32e48fe85c928cd3e6255"


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('archive',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--metals',type=float,nargs='+',default=[.01,.02,.03],
                   help='ordered source Z coordinates; every requested source cell must exist')
    p.add_argument('--hydrogen-endpoints',action='store_true',
                   help='also import each Z plane\'s actual hydrogen-rich endpoint from the full archive')
    a=p.parse_args()
    digest=hashlib.sha256(a.archive.read_bytes()).hexdigest()
    if digest not in (SOURCE_SHA256,SUBSET_SHA256):raise ValueError('unrecognized source archive')
    if a.hydrogen_endpoints and digest!=SOURCE_SHA256:raise ValueError('hydrogen endpoints require the full pinned source archive')
    a.output.mkdir(parents=True,exist_ok=True)
    xs=[0,.1,.2,.35,.5,.7,.8,.9,.95];zs=a.metals
    if (sorted(set(zs))!=zs or not all(math.isfinite(z) and 0<=z<1-max(xs) for z in zs)):
        raise ValueError('invalid source metal coordinates')
    tt=list(range(2000,3701,10))+list(range(3720,4501,20));rr=list(range(-8000,6001,200))
    manifest=[f'EMBER_OPACITY_MIXTURE 1 {len(zs)} logR AESOPUS_2.1_GS98_gas'];planes=0
    with ZipFile(a.archive) as source:
        for z in zs:
            # The zero-Z member is named X0.9999999, while its six-decimal
            # header rounds X to 1.000000. Retain the member's coordinate;
            # do not manufacture a pure-H cell or change its opacity values.
            endpoint=.9999999 if z==0 else 1-z
            plane_xs=xs+[endpoint] if a.hydrogen_endpoints else xs
            planes+=len(plane_xs)
            lines=[f'{len(plane_xs)} {len(tt)} {len(rr)} AESOPUS 2.1 GS98 original cells, Z={z}',
                ' '.join(f'{r/1000:.3f}' for r in rr),' '.join(f'{t/1000:.3f}' for t in tt)]
            for x in plane_xs:
                cells={}
                for region in ['lowT','highT','highR']:
                    x_name='0.9999999' if x==.9999999 and z==0 else f'{x:g}'
                    member=f'03-GS98/GS98_a0.0_OPALZ_{region}/aesopus2.0_gasbroad_GS98_Z{z:.6f}_X{x_name}.tab'
                    text=source.read(member).decode('ascii')
                    if 'Calculations performed with AESOPUS 2.1' not in text:raise ValueError('unexpected source version')
                    header=next(line for line in text.splitlines() if line.startswith('# TABLE'))
                    def field(name):return float(re.search(r'\b'+name+r'=\s*([\d.E+-]+)',header)[1])
                    header_x=1. if z==0 and x==.9999999 else x
                    if field('X')!=header_x or field('Zref')!=z or abs(field('Z')-z)>1.1e-8:raise ValueError('wrong source composition')
                    density=[r for r in rr if (r>1000)==(region=='highR')]
                    axis=re.search(r'log10\(R\) range: nre=\s*(\d+)\s+values from\s+([\d.+-]+)\s+to\s+([\d.+-]+)\s*\n#\s+in steps of\s+([\d.]+)',text)
                    if not axis or (int(axis[1]),round(float(axis[2])*1000),round(float(axis[3])*1000),round(float(axis[4])*1000))!=(len(density),density[0],density[-1],200):raise ValueError('wrong density axis')
                    for line in text.splitlines():
                        if not line.strip() or line.startswith('#'):continue
                        row=line.split();temp=round(float(row[0])*1000)
                        if len(row)!=len(density)+1:raise ValueError('truncated row')
                        for r,value in zip(density,row[1:],strict=True):
                            if (temp,r) in cells or not math.isfinite(float(value)) or abs(float(value))>50:raise ValueError('invalid cell')
                            cells[temp,r]=value
                if set(cells)!={(t,r) for t in tt for r in rr}:raise ValueError('missing source cells')
                lines.append(f'{x:.17g} {z:.17g}')
                lines.extend(' '.join(cells[t,r] for r in rr) for t in tt)
            name=f'aesopus21_gs98_mixture_z{round(z*1000):03d}.dat'
            (a.output/name).write_text('\n'.join(lines)+'\n');manifest.append(f'{z:.17g} "{name}"')
    (a.output/'aesopus21_gs98_mixture.dat').write_text('\n'.join(manifest)+'\n')
    print(f'Imported {planes} (X,Z) planes x 211 T x 71 density original AESOPUS cells')

if __name__=='__main__':main()
