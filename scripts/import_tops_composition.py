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


def validate_mixture(text, expected, *, dimensions=(50,71)):
    """Check every returned element, allowing omitted H only for exactly zero H.

    TOPS omits zero-abundance hydrogen from its normalized composition and
    reports 20 materials for our helium/GS98-metal endpoint. Positive-hydrogen
    mixtures still require all 21 elements; no absent metal is inferred.
    """
    if (len(dimensions)!=2 or any(type(n) is not int or n<1 for n in dimensions)
            or dimensions[1]>100):
        raise ValueError('invalid requested dimensions')
    header=re.search(r'Number of T\s*=\s*(\d+)\s+Number of rho\s*=\s*(\d+)\s+Number of materials\s*=\s*(\d+)',text)
    if not header or tuple(map(int,header.groups()[:2]))!=tuple(dimensions):
        raise ValueError('unexpected dimensions')
    rows=text.split('No. Fraction Mass Fraction  At. No.  Chem. Sym.  Mat ID.\n')[1].split('Temperature grid')[0]
    elements={}
    for row in rows.splitlines():
        if not row.strip():continue
        values=row.split();element=values[3];fraction=float(values[1])
        if element in elements or not math.isfinite(fraction) or fraction<0:
            raise ValueError('invalid source element row')
        elements[element]=fraction
    x=expected['X']
    required=set(expected['metals'])|{'He'}
    if x!=0 or 'H' in elements:required.add('H')
    hydrogen=elements.get('H',0.)
    hydrogen_ok=hydrogen==0 if x==0 else abs(hydrogen-x)<=min(1e-6,5e-5*x)
    if (len(elements)!=int(header.group(3)) or set(elements)!=required
            or not hydrogen_ok or abs(elements['He']-(1-expected.get('Z',.02)-x))>1e-6):
        raise ValueError('wrong source mixture')
    # Verify every metal against the original GS98 request, not just sum Z.
    for element,value in expected['metals'].items():
        if abs(elements[element]-value)>max(1e-9,5e-5*value):
            raise ValueError(f'wrong source metal abundance: {element}')
    return elements


def read(path, expected, *, dimensions=(50,71)):
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected['sha256']:
        raise ValueError(f'unexpected checksum: {path}')
    text=raw.decode()
    validate_mixture(text,expected,dimensions=dimensions)
    prefix=text.split('Normalized composition')[0]
    marker='Temp        Den Req     Den Used\n'
    if marker in prefix:
        warning=prefix.split(marker)[1]
        excluded={tuple(map(float,row.split()[:2])) for row in warning.splitlines() if row.strip()}
    elif re.search(r'warning|den req|den used|substitut',prefix,re.IGNORECASE):
        raise ValueError('unrecognized source density warning')
    else:
        # A supported single-temperature request can have no warning block.
        excluded=set()
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
    if len(tt)!=dimensions[0] or len(rr)!=dimensions[1] or set(cells)!={(t,r) for t in tt for r in rr}:
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
