#!/usr/bin/env python3
"""Independent off-grid FreeEOS queries for the baryonic H1/He3 EOS audit.

EPS is moles per gram of each element, summed over its isotopes. FreeEOS
does not distinguish helium isotope entropy, so this audit covers P, E and
thermal responses, not the added ideal isotope mixing/translation entropy.
Build the external probe using the instructions in docs/FREEEOS.md.
"""
import argparse
import math
from pathlib import Path
import subprocess


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('probe',type=Path);ap.add_argument('output',type=Path)
    a=ap.parse_args()
    points=[(3.371,-6.773),(3.731,-4.881),(4.173,-5.19),(4.461,-3.17),
            (5.173,-1.177),(6.119,.17),(6.713,2.31)]
    rows=[]
    for x,y3 in [(.625,0),(.675,.001),(.7,.004),(.725,.002)]:
        # Metals remain the He4 proxy, so effective He4 is 1-X-Y3.
        eps=[x,y3/3+(1-x-y3)/4]+[0.]*18
        for lt,lr in points:
            t=10**lt;r=10**lr
            request=' '.join(map(str,eps))+'\n3 1 -2\n'+f'{math.log(r):.17g} {math.log(t):.17g}\n'
            p=subprocess.run([str(a.probe.resolve())],input=request,text=True,capture_output=True,check=True,timeout=60)
            values=list(map(float,p.stdout.split()))
            if len(values)!=22 or values[0]!=0 or not all(math.isfinite(v) for v in values):
                raise ValueError('invalid FreeEOS source response')
            if abs(values[2]/r-1)>1e-10 or abs(values[3]/t-1)>1e-10:
                raise ValueError('FreeEOS did not reach requested density/temperature')
            rows.append([x,y3,t,r,values[4],values[5],values[12],values[13],values[14]])
    a.output.write_text(''.join(' '.join(f'{v:.17g}' for v in row)+'\n' for row in rows))
    print(f'{len(rows)} independent FreeEOS elemental-number-density states')


if __name__=='__main__':main()
