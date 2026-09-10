#!/usr/bin/env python3
"""Retain original Ioffe 2021 conductivity cells for light ions (Z=1..12)."""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('sources',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--metal-coverage',action='store_true',help='retain original source ions through Zn (Z=30), bracketing every GS98 element')
    a=p.parse_args();manifest=json.loads((a.sources/'manifest.json').read_text());a.output.mkdir(parents=True,exist_ok=True)
    charges=[1,2,3,4,6,8,12];masses=[1,4,7,9,12,16,24]
    if a.metal_coverage:charges += [16,20,26,30];masses += [32,40,56,64]
    for file,meta in manifest.items():
        packed=(a.sources/file).read_bytes();raw=gzip.decompress(packed)
        if hashlib.sha256(packed).hexdigest()!=meta['gzip_sha256'] or hashlib.sha256(raw).hexdigest()!=meta['raw_sha256']:
            raise ValueError('unrecognized source checksum')
        rows=raw.decode('ascii').splitlines()[1:]
        if len(rows)!=15*65:raise ValueError('unexpected table dimensions')
        temps=list(map(float,rows[0].split()[1:]));densities=[float(rows[i].split()[0]) for i in range(1,65)]
        if len(temps)!=19 or densities[0]!=-6 or densities[-1]!=9.75:raise ValueError('wrong axes')
        out=[f'EMBER_CONDUCTIVITY 1 {len(charges)} 19 64 {file.removesuffix(".dat.gz")}',
            ' '.join(map(str,temps)),' '.join(map(str,densities))]
        for i,(z,mass) in enumerate(zip(charges,masses,strict=True)):
            header=list(map(float,rows[i*65].split()))
            if header!=[z]+temps:raise ValueError('source charge/temperature mismatch')
            out.append(f'{z} {mass}')
            for j,rho in enumerate(densities):
                row=rows[i*65+j+1].split()
                if len(row)!=20 or float(row[0])!=rho or any(not math.isfinite(float(v)) for v in row):raise ValueError('invalid source row')
                out.append(' '.join(row[1:]))
        target=file.removesuffix('.dat.gz')+('_metals' if a.metal_coverage else '')+'.dat'
        (a.output/target).write_text('\n'.join(out)+'\n')
        print(target, f'{len(charges)} x 19 x 64 untouched source conductivity values')
if __name__=='__main__':main()
