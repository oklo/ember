#!/usr/bin/env python3
"""Join original AESOPUS 2.1 GS98 Z=.02 cells into ember's rectangular format.

Input: https://stev.oapd.inaf.it/aesopus_2.1/tables/SOLAR/03-GS98.zip
The archive's filenames retain 'aesopus2.0'; the headers identify version 2.1.
These are gas tables, without grain opacity. No missing cells are filled.
"""
import hashlib
import math
from pathlib import Path
import re
import sys
import zipfile

SOURCE_SHA256 = "9379a46e89b0e2c134d018986e8499e27cf6b45038de493f57b2f766ffaaa830"


def main():
    source, target = map(Path, sys.argv[1:])
    if hashlib.sha256(source.read_bytes()).hexdigest() != SOURCE_SHA256:
        raise ValueError("Source differs from the documented AESOPUS 2.1 GS98 archive")
    xs = [0, 0.1, 0.2, 0.35, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98]
    # Integer coordinates avoid mistaking float rounding for missing nodes.
    temperatures = list(range(2000, 3701, 10)) + list(range(3720, 4501, 20))
    densities = list(range(-8000, 6001, 200))
    output = [f"{len(xs)} {len(temperatures)} {len(densities)} AESOPUS 2.1 gas GS98 Z=0.020; original cells",
              " ".join(f"{r / 1000:.3f}" for r in densities),
              " ".join(f"{t / 1000:.3f}" for t in temperatures)]
    with zipfile.ZipFile(source) as archive:
        for x in xs:
            cells = {}
            for region in ("lowT", "highT", "highR"):
                member = (f"03-GS98/GS98_a0.0_OPALZ_{region}/"
                          f"aesopus2.0_gasbroad_GS98_Z0.020000_X{x:g}.tab")
                text = archive.read(member).decode("ascii")
                if "Calculations performed with AESOPUS 2.1" not in text:
                    raise ValueError("Unexpected source version")
                header = next(line for line in text.splitlines() if line.startswith("# TABLE"))
                # The source prints Zref=.02000000, but some actual-Z fields
                # are .02000001 (a 0.5 ppm normalization difference). Preserve
                # every opacity value and identify the grid by nominal Zref.
                if float(re.search(r"\bX=\s*([\d.E+-]+)", header)[1]) != x or float(
                        re.search(r"\bZref=\s*([\d.E+-]+)", header)[1]) != .02 or float(
                        re.search(r"\bZ=\s*([\d.E+-]+)", header)[1]) not in (.02, .02000001):
                    raise ValueError("Composition disagrees with archive member name")
                rr = [r for r in densities if (r > 1000) == (region == "highR")]
                axis = re.search(r"log10\(R\) range: nre=\s*(\d+)\s+values from\s+([\d.+-]+)\s+to\s+([\d.+-]+)\s*\n#\s+in steps of\s+([\d.]+)", text)
                if not axis or (int(axis[1]), round(float(axis[2]) * 1000),
                                round(float(axis[3]) * 1000), round(float(axis[4]) * 1000)) != (
                                    len(rr), rr[0], rr[-1], 200):
                    raise ValueError("Unexpected density axis")
                for line in text.splitlines():
                    if not line.strip() or line.startswith("#"):
                        continue
                    row = line.split()
                    if len(row) != len(rr) + 1:
                        raise ValueError("Truncated source row")
                    t = round(float(row[0]) * 1000)
                    for r, value in zip(rr, row[1:], strict=True):
                        if (t, r) in cells or not math.isfinite(float(value)) or abs(float(value)) > 50:
                            raise ValueError("Duplicate or invalid source cell")
                        cells[t, r] = value
            if set(cells) != {(t, r) for t in temperatures for r in densities}:
                raise ValueError("Source grid has missing or unexpected cells")
            output.append(f"{x:.4f} 0.0200")
            output.extend(" ".join(cells[t, r] for r in densities) for t in temperatures)
    target.write_text("\n".join(output) + "\n")


if __name__ == "__main__":
    main()
