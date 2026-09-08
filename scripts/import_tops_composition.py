#!/usr/bin/env python3
"""Import composition planes from checksum-pinned, individually verified TOPS data.

Only original, un-substituted cells in the common density support are used.
The output uses the existing fixed-Z, multiple-X, native-density format.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re


def read(path, expected):
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected['sha256']:
        raise ValueError(f'unexpected checksum: {path}')
    text=raw.decode()
    if 'Number of T =  50  Number of rho =  71  Number of materials =  21' not in text:
        raise ValueError('unexpected dimensions')
    warning=text.split('Temp        Den Req     Den Used\n')[1].split('Normalized composition')[0]
    excluded={tuple(map(float,row.split()[:2])) for row in warning.splitlines() if row.strip()}
    rows=text.split('No. Fraction Mass Fraction  At. No.  Chem. Sym.  Mat ID.\n')[1].split('Temperature grid')[0]
    elements={r.split()[3]:float(r.split()[1]) for r in rows.splitlines() if r.strip()}
    x=expected['X']
    if len(elements)!=21 or abs(elements['H']-x)>1e-6 or abs(elements['He']-(.98-x))>1e-6:
        raise ValueError('wrong source mixture')
    # Verify every metal against the original GS98 request, not just sum Z.
    for element,value in expected['metals'].items():
        if abs(elements[element]-value)>max(1e-9,5e-5*value):
            raise ValueError(f'wrong source metal abundance: {element}')
    cells={}
    for part in text.split('Density     Ross opa    Planck opa  No. Free    Av Sq Free  T=  ')[1:]:
        lines=part.splitlines();t=float(lines[0])
        for line in lines[1:]:
            values=line.split()
            if len(values)!=5 or not all(re.fullmatch(r'[0-9.E+-]+',v) for v in values):break
            rho,kappa=map(float,values[:2])
            if (t,rho) in cells or not rho>0 or not math.isfinite(kappa) or kappa<=0:
                raise ValueError('bad source cell')
            cells[t,rho]=kappa
    tt=sorted({t for t,r in cells});rr=sorted({r for t,r in cells})
    if len(tt)!=50 or len(rr)!=71 or set(cells)!={(t,r) for t in tt for r in rr}:
        raise ValueError('missing source cells')
    return tt,rr,cells,excluded


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('manifest',type=Path);ap.add_argument('output',type=Path)
    a=ap.parse_args();spec=json.loads(a.manifest.read_text())
    planes=spec['planes'];xs=[p['X'] for p in planes]
    if len(xs)<2 or any(b<=v for v,b in zip(xs,xs[1:])):raise ValueError('invalid X grid')
    sources=[read(a.manifest.parent/p['file'],dict(p,metals=spec['metals'])) for p in planes]
    tt,rr=sources[0][:2]
    if any(t!=tt or r!=rr for t,r,_,_ in sources):raise ValueError('unaligned source grids')
    a.output.mkdir(parents=True,exist_ok=True)
    for label,temps in [('low',tt),('high',[t for t in tt if t>=.025])]:
        densities=[]
        for r in rr:
            if any((t,r) in excluded for _,_,_,excluded in sources for t in temps):break
            densities.append(r)
        if len(densities)<4:raise ValueError('insufficient common density support')
        lines=[f'{len(xs)} {len(temps)} {len(densities)} TOPS ATOMIC GS98 variable X Z=.02 {label}; native log rho; all substituted cells excluded',
               ' '.join(f'{math.log10(r):.17g}' for r in densities),
               ' '.join(f'{math.log10(t*1e3*1.602176634e-12/1.380649e-16):.17g}' for t in temps)]
        for x,(_,_,cells,excluded) in zip(xs,sources):
            lines.append(f'{x:.17g} 0.02')
            lines.extend(' '.join(f'{math.log10(cells[t,r]):.17g}' for r in densities) for t in temps)
        (a.output/f'tops_gs98_composition_z020_{label}.dat').write_text('\n'.join(lines)+'\n')
        print(f'{label}: {len(xs)} X planes, {len(temps)} temperatures, {len(densities)} original densities')


if __name__=='__main__':main()
