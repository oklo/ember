#!/usr/bin/env python3
"""Generate independent scalar HRW plasma-fit values for the runtime port.

Uses the existing offline implementation of Haft, Raffelt & Weiss (1994),
equations23--27, not the C++ loss implementation. This tests a formula port;
it is not a fresh exact plasma calculation or a bound on omitted channels.
"""
import argparse
from pathlib import Path
from audit_extended_track import plasmon


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);a=p.parse_args()
    rows=['# HRW plasma fit; T[K] rho[g/cm3] Ye eps[erg/g/s]',
          '# Source: https://arxiv.org/abs/astro-ph/9309014; scripts/generate_neutrino_reference.py']
    for ye in [.5, .85, 1.]:
        for T,rho in [(1e6,1), (4.5e6,360), (1e7,1e3), (1e7,1e5), (1e8,1e6), (1e9,1e7), (100,1e5)]:
            rows.append(' '.join(format(v,'.17g') for v in [T,rho,ye,plasmon(T,rho,ye)]))
    a.output.write_text('\n'.join(rows)+'\n')


if __name__=='__main__':main()
