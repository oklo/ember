#!/usr/bin/env python3
"""Import two un-clamped rectangles from the archived TOPS GS98 response.

Density warnings are data: every warned (T,rho) cell must be excluded.
The native T[keV], rho[g/cm3] and kappa[cm2/g] values are converted to logs,
without resampling, extrapolation, or including substituted densities.
"""
import hashlib
import math
from pathlib import Path
import re
import sys

SHA256 = "6e99ef76f5409926168a7defe7205c73d12e42166d8a388daffe3a4d44542344"


def main():
    source, target = map(Path, sys.argv[1:])
    if hashlib.sha256(source.read_bytes()).hexdigest() != SHA256:
        raise ValueError("TOPS response differs from documented source")
    target.mkdir(parents=True, exist_ok=True)
    text = source.read_text()
    if "Number of T =  50  Number of rho =  71  Number of materials =  21" not in text:
        raise ValueError("Unexpected TOPS grid or mixture")
    warning = text.split("Temp        Den Req     Den Used\n")[1].split("Normalized composition")[0]
    excluded = {tuple(map(float, r.split()[:2])) for r in warning.splitlines() if r.strip()}
    # Validate the printed mass fractions, allowing only their display precision.
    rows = text.split("No. Fraction Mass Fraction  At. No.  Chem. Sym.  Mat ID.\n")[1].split("Temperature grid")[0].splitlines()
    elements = {r.split()[3]: float(r.split()[1]) for r in rows if r.strip()}
    if len(elements) != 21 or elements['H'] != .7 or elements['He'] != .28 or abs(
            sum(v for k, v in elements.items() if k not in ('H', 'He')) - .02) > 1e-6:
        raise ValueError("Returned mixture is not the requested X=.7, Z=.02")
    cells = {}
    for section in text.split("Density     Ross opa    Planck opa  No. Free    Av Sq Free  T=  ")[1:]:
        lines = section.splitlines(); T = float(lines[0])
        for line in lines[1:]:
            fields = line.split()
            if len(fields) != 5 or not all(re.fullmatch(r"[0-9.E+-]+", v) for v in fields): break
            rho, kappa = map(float, fields[:2])
            if (T, rho) in cells or not (rho > 0 and kappa > 0 and math.isfinite(kappa)):
                raise ValueError("Invalid or duplicate source cell")
            cells[T, rho] = kappa
    tt = sorted({t for t, r in cells}); rr = sorted({r for t, r in cells})
    if len(cells) != 50 * 71 or set(cells) != {(t, r) for t in tt for r in rr}:
        raise ValueError("Missing or unexpected TOPS cells")
    for label in ('low', 'high'):
        temps = tt if label == 'low' else [t for t in tt if t >= .025]
        densities = [r for r in rr if r <= 251.19] if label == 'low' else rr
        if excluded & {(t, r) for t in temps for r in densities}:
            raise ValueError("Requested rectangle would include server-clamped values")
        out = [f"1 {len(temps)} {len(densities)} TOPS ATOMIC GS98 X=.7 Z=.02 {label} rectangle; log rho axis",
               ' '.join(f'{math.log10(r):.17g}' for r in densities),
               ' '.join(f'{math.log10(t*1e3*1.602176634e-12/1.380649e-16):.17g}' for t in temps),
               '0.7 0.02']
        out += [' '.join(f'{math.log10(cells[t,r]):.17g}' for r in densities) for t in temps]
        (target / f'tops_gs98_x070_z020_{label}.dat').write_text('\n'.join(out)+'\n')
        print(f'{label}: {len(temps)} x {len(densities)} original cells; all clamped points excluded')


if __name__ == '__main__':
    main()
