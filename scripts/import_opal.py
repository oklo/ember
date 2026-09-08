#!/usr/bin/env python3
"""Extract unmodified, complete Z=0.020 OPAL cells from Boothroyd's GS98hz.gz.

curl -fLO https://www.cita.utoronto.ca/~boothroy/data/kappa/GS98hz.gz
python3 scripts/import_opal.py GS98hz.gz data/opacity/opal_gs98_z020.dat
No interpolation, extrapolation, or replacement of missing cells.
"""
import gzip
import hashlib
import pathlib
import re
import sys

SOURCE_SHA256 = "dbb0ebb9a231e30201f7347f2033c72783ed8acc1f6abcaa57bb89f2df0b2ca3"


def main():
    source, target = map(pathlib.Path, sys.argv[1:])
    raw = gzip.decompress(source.read_bytes())
    if hashlib.sha256(raw).hexdigest() != SOURCE_SHA256:
        raise ValueError("Source differs from the documented OPAL GS98 distribution")
    tables = []
    for block in re.split(r"(?m)^TABLE #", raw.decode("ascii"))[1:]:
        header, *lines = block.splitlines()
        x, z = (float(re.search(rf"{key}=([0-9.]+)", header)[1]) for key in ("X", "Z"))
        if z != 0.020:
            continue
        axis = next(line for line in lines if line.startswith("logT")).split()[1:]
        rows = [line.split() for line in lines if re.match(r"^[0-9]\.[0-9]", line)]
        rows = [row for row in rows if 4.0 <= float(row[0]) <= 7.1]
        if len(rows) != 52 or any(len(row) != 20 or "9.999" in row for row in rows):
            raise ValueError("The selected subset contains missing cells")
        tables.append((x, axis, rows))
    tables.sort()
    if [x for x, _, _ in tables] != [0.0, 0.1, 0.2, 0.35, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98]:
        raise ValueError("Unexpected hydrogen axis")
    temps = [row[0] for row in tables[0][2]]
    axis = tables[0][1]
    output = [f"{len(tables)} {len(temps)} {len(axis)} OPAL GS98 Z=0.020; original cells, logT=4.0..7.1",
              " ".join(axis), " ".join(temps)]
    for x, density_axis, rows in tables:
        if density_axis != axis or [row[0] for row in rows] != temps:
            raise ValueError("Inconsistent axes")
        output.append(f"{x:.4f} 0.0200")
        output.extend(" ".join(row[1:]) for row in rows)
    target.write_text("\n".join(output) + "\n")


if __name__ == "__main__":
    main()
