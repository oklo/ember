#!/usr/bin/env python3
"""Independent conductivity samples from the public Ioffe conduct21.f code.

Build scripts/conduction_reference_probe.f90 together with that external code.
The source SHA is recorded with the versioned reference in data/conduction/README.md.
"""
import argparse
import math
from pathlib import Path
import subprocess

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('probe',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    points=[(1,1,8.1e5,.1),(1,1,4.6e6,330),(1,1,6.2e6,500),(2,4,8.1e5,.2),(2,4,4.6e6,660),(2,4,6.2e6,1000),
        (6,12,4.6e6,660),(7,14,4.6e6,660),(8,16,4.6e6,660),(10,20,4.6e6,660)]
    rows=[]
    for mode in [-1,0,1]:
        for z,mass,t,rho in points:
            p=subprocess.run([str(a.probe.resolve())],input=f'{t} {rho} {z} {mass} {mode}\n',text=True,capture_output=True,check=True)
            k=float(p.stdout.strip())
            if not math.isfinite(k) or k<=0:raise ValueError('invalid conductivity response')
            rows.append([mode,z,mass,t,rho,k])
    a.output.write_text(''.join(' '.join(f'{v:.17g}' for v in row)+'\n' for row in rows))
    print(len(rows),'independent conductivity queries')
if __name__=='__main__':main()
