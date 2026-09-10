#!/usr/bin/env python3
"""Extract one plotted trial state from a native checkpoint, never a restart."""
import argparse
import hashlib
import json
import math
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('checkpoint', type=Path); p.add_argument('output', type=Path)
    a = p.parse_args(); lines = a.checkpoint.read_text().splitlines()
    if lines[0] != 'EMBER_EVOLUTION_CHECKPOINT 1':
        raise ValueError('unrecognized checkpoint')
    index = 8 + int(lines[7]); header = lines[index].split(); points = int(header[0])
    rows = [list(map(float, line.split())) for line in lines[index+1:index+1+points]]
    if len(rows) != points or any(len(r) != 15 or not all(math.isfinite(v) for v in r) or abs(sum(r[7:])-1) > 2e-12 for r in rows):
        raise ValueError('invalid checkpoint profile')
    if rows[-1][0] != float(header[1]):
        raise ValueError('surface mass differs')
    radius = math.exp(rows[-1][1]); luminosity = rows[-1][4]
    result = {'scope': 'Last saved accepted trial state; requested endpoint not reached. Converted only for plotting, not restart.',
        'input': str(a.checkpoint), 'input_sha256': hashlib.sha256(a.checkpoint.read_bytes()).hexdigest(),
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'points': points, 'accepted_steps': int(header[4]), 'rejected_steps': int(header[5]),
        'age_yr': float(header[2])/31557600, 'R_Rsun': radius/6.957e10, 'L_Lsun': luminosity/3.828e33,
        'Teff_K': (luminosity/(4*math.pi*radius**2*5.670374419e-5))**.25,
        'central_X': rows[0][7], 'central_Y3': rows[0][8],
        'surface_X': rows[-1][7], 'surface_Y3': rows[-1][8]}
    a.output.write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
