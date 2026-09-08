#!/usr/bin/env python3
"""Extract original CMS19 TP density/entropy cells; do not repair energy.

The static EOS uses pressure and entropy. The source energy has join defects
and is deliberately not enabled for time-dependent energy calculations.
"""
import hashlib
import math
from pathlib import Path
import sys
import tarfile

SHA256 = "736de2a0b26c02b897fbd906504bb6aa144f08a6d0f3cff36da62970054d149b"


def main():
    source, target = map(Path, sys.argv[1:])
    if hashlib.sha256(source.read_bytes()).hexdigest() != SHA256:
        raise ValueError("CMS19 archive differs from documented source")
    target.mkdir(parents=True, exist_ok=True)
    with tarfile.open(source) as archive:
        for species in ("H", "HE"):
            text = archive.extractfile(f"DirTABLES-EOS2019/TABLE_{species}_TP_v1").read().decode("ascii")
            rows = [r.split() for r in text.splitlines() if r.strip() and not r.startswith("#")]
            if len(rows) != 121 * 441:
                raise ValueError("Unexpected TP grid size")
            output = [f"121 441 CMS19 {species}: original log10 rho and log10 S in published units",
                      " ".join(f"{2 + .05 * i:.2f}" for i in range(121)),
                      " ".join(f"{-9 + .05 * j:.2f}" for j in range(441))]
            for n, row in enumerate(rows):
                if len(row) != 10 or not all(math.isfinite(float(v)) for v in row):
                    raise ValueError("Malformed source cell")
                if abs(float(row[0]) - (2 + .05 * (n // 441))) > 1e-10 or abs(
                        float(row[1]) - (-9 + .05 * (n % 441))) > 1e-10:
                    raise ValueError("Source axes disagree with expected grid")
                # Preserve even unphysical source corners; runtime validity
                # masks exclude them. No interpolation/filling during import.
                output.append(f"{row[2]} {row[4]}")
            (target / f"cms19_{species.lower()}_tp.dat").write_text("\n".join(output) + "\n")


if __name__ == "__main__":
    main()
