#!/usr/bin/env python3
"""Extract untouched AMES-COND tau=100 cells from checksum-pinned MESA data.

The version-2 output explicitly declares a solar-mixture proxy. Its
composition is the permitted interior query, not a reconstructed COND mixture.
"""
import argparse
import hashlib
import math
from pathlib import Path
from stellar_composition import interior_composition

HASHES = {
    'tau100_Pgas.data': '77606d68e2f0c02602624e04e31792c5150c7a15940b53cc8a364ba3682f5900',
    'tau100_T.data': 'e9551fac9e2d8d036ab94c55ed87cb66811d2a0c54e82d27c846c1aabdea4c3a',
}


def read(path):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != HASHES[path.name]:
        raise ValueError(f'{path}: unexpected source checksum')
    rows = [[float(x) for x in line.split()] for line in data.decode().splitlines()
            if line.strip() and not line.startswith('#')]
    if len(rows) != 110 or any(len(row) != 14 for row in rows):
        raise ValueError('unexpected source dimensions')
    return {row[0]: row[1:] for row in rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source_dir', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    pressure, temperature = [read(args.source_dir / name) for name in HASHES]
    # create_tau100.f90 reads COND at g=2.5..6, T=100..3600;
    # CK overwrite starts at 3500, and transition smoothing at 3400.
    # This smaller cool-dwarf rectangle excludes every altered region.
    teffs = list(range(1800, 3301, 100))
    gravities = [3.5, 4., 4.5, 5., 5.5, 6.]
    lines = ['EMBER_ATMOSPHERE 2',
             'source "AMES-COND-2000 / Allard et al. 2001; MESA fd396fd73d3f936da8063ffdf9d92361882eb557; unmodified cool-dwarf tau100 cells"',
             'tau 100',
             'composition_proxy "Native GN93 solar atmosphere used for the specified X=.7 Z=.02 interior; detailed mixtures and helium abundance are not matched"',
             'composition ' + ' '.join(format(x, '.17g') for x in interior_composition()),
             'log_teff 16 ' + ' '.join(format(math.log10(t), '.17g') for t in teffs),
             'log_g 6 ' + ' '.join(map(str, gravities)), 'data']
    for t in teffs:
        for g in gravities:
            i = round(g*2)
            T, P = temperature[t][i], pressure[t][i]
            if not (math.isfinite(T) and math.isfinite(P) and T > 0 and P > 0):
                raise ValueError('invalid selected cell')
            lines.append(f'{math.log10(T):.17g} {math.log10(P):.17g}')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text('\n'.join(lines)+'\n')


if __name__ == '__main__':
    main()
