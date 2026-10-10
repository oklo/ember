"""Normalize downloaded catalogues without assigning stellar masses.

Run with the project's science Python and this study's deps on PYTHONPATH.
Gaia source IDs stay strings. All published membership columns are retained;
no probabilities from different selection methods are combined here.
"""
from pathlib import Path
import csv
from collections import Counter
import hashlib
import json

import numpy as np
import pandas as pd
from astropy.table import Table

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
DEST = ROOT / "catalogues"
DEST.mkdir(exist_ok=True)


def read_asu(path):
    lines = [line for line in path.read_text().splitlines()
             if line.strip() and not line.startswith("#")]
    # Only single-table ASU exports are accepted here. Unit and ruler rows
    # are not data. A malformed/multitable file must fail rather than slip in.
    names = lines[0].split("\t")
    assert len(names) == len(set(names)), path
    assert all(not cell.strip("- ") for cell in lines[2].split("\t")), path
    rows = []
    for line in lines[3:]:
        cells = [cell.strip() for cell in line.split("\t")]
        assert len(cells) == len(names), (path, line[:80])
        rows.append(dict(zip(names, cells)))
    return names, rows


def write_csv(path, names, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)


receipt = {"stage": "catalogue extraction; no mass inference or completeness correction",
           "inputs": {}, "tables": {}, "dance_scaling_validation": {}}
exports = {
    "hyades_liu2025_members": 976,
    "pleiades_liu2025_members": 1763,
    "hyades_kos2024": 4902,
    "pleiades_kos2024": 3392,
    "pleiades_olivares2018": 5444,
    "pleiades_liu2025_binaries": 6041,
}
for stem, expected in exports.items():
    source = RAW / (stem + ".tsv")
    names, rows = read_asu(source)
    assert len(rows) == expected, (stem, len(rows), expected)
    if "GaiaDR3" in names:
        assert all(row["GaiaDR3"].isdigit() for row in rows)
        assert len({row["GaiaDR3"] for row in rows}) == len(rows), stem
    write_csv(DEST / (stem + ".csv"), names, rows)
    receipt["inputs"][source.name] = hashlib.sha256(source.read_bytes()).hexdigest()
    receipt["tables"][stem] = {"rows": len(rows), "columns": names}
    if "Memb" in names:
        receipt["tables"][stem]["tidal_radius_flags"] = {
            flag: sum(row["Memb"] == flag for row in rows)
            for flag in sorted({row["Memb"] for row in rows})}

# The VizieR TAP join returns several Bouy columns in storage units despite
# labeling them with physical units. The Olivares probability columns are
# already decoded by the VOTable reader. Validate against an independently
# retrieved ASU table, whose columns are in physical units.
source = RAW / "pleiades_dance_join.vot"
dance = Table.read(source, format="votable")
join_counts = Counter(map(str, dance["BouyID"]))
# Rounded position-based identifiers are not unique in these catalogues.
# Preserve each joined row and flag the ambiguity; never count joined rows
# directly as stars. Later crossmatching needs the original coordinates.
dance["join_rows_for_id"] = [join_counts[str(x)] for x in dance["BouyID"]]
scales = {}
for name in ("RAJ2000", "DEJ2000", "e_RAJ2000", "e_DEJ2000"):
    scales[name] = 3.6e6
for name in ("pmRA", "pmDE", "e_pmRA", "e_pmDE", "chi2r"):
    scales[name] = 100
for name in ("Epoch", "b_Epoch", "B_Epoch", "Pmb"):
    scales[name] = 1e4
for band in ("u", "g", "r", "i", "z", "Y", "J", "H", "K"):
    for prefix in ("", "e_"):
        scales[prefix + band + "mag"] = 1000
assert float(dance["RAJ2000"][0]) > 360, "TAP format changed; reassess scaling"
for name, divisor in scales.items():
    dance[name] = dance[name].astype(float) / divisor
assert np.all((dance["RAJ2000"] >= 0) & (dance["RAJ2000"] < 360))
assert np.all((dance["DEJ2000"] >= -90) & (dance["DEJ2000"] <= 90))
assert np.all((dance["pc"] >= 0) & (dance["pc"] <= 1))
_, reference_rows = read_asu(RAW / "pleiades_bouy2015.tsv")
reference = {row["DANCe"]: row for row in reference_rows}
matches = [(row, reference[str(row["DANCe"])]) for row in dance
           if str(row["DANCe"]) in reference]
assert len(matches) >= 20
tolerances = {"RAJ2000": 6e-7, "DEJ2000": 6e-7, "pmRA": 1e-9,
              "pmDE": 1e-9, "Jmag": 1e-9, "Kmag": 1e-9,
              "Pmb": 1e-9, "Lbol": 1e-6}
for name, tolerance in tolerances.items():
    differences = [abs(float(row[name]) - float(ref[name]))
                   for row, ref in matches if ref[name] and not np.ma.is_masked(row[name])]
    assert differences and max(differences) < tolerance, (name, differences)
    receipt["dance_scaling_validation"][name] = {
        "comparisons": len(differences), "max_abs_difference": max(differences)}
dance.write(DEST / "pleiades_dance_physical.csv", format="ascii.csv", overwrite=True)
# Replace the previously created but incorrectly scaled convenience export.
dance.write(RAW / "pleiades_dance_join_physical.csv", format="ascii.csv", overwrite=True)
receipt["inputs"][source.name] = hashlib.sha256(source.read_bytes()).hexdigest()
receipt["tables"]["pleiades_dance_physical"] = {
    "rows": len(dance), "storage_divisors": scales,
    "unique_position_ids": len(join_counts),
    "ambiguous_position_ids": {k: v for k, v in join_counts.items() if v > 1},
    "pc_gt_0.84": int(np.sum(dance["pc"] > .84)),
    "finite_lbol": int(np.sum(~np.ma.getmaskarray(dance["Lbol"]))) }

# The original numeric Sarro ID resolves the ambiguous rounded names exactly.
prob = pd.read_csv(DEST / "pleiades_olivares2018.csv", dtype={"SarroID": str})
phot = pd.read_csv(DEST / "pleiades_dance_physical.csv", dtype={"__SBB2014_": str})
phot = phot.drop_duplicates("recno")
phot = phot.drop(columns=["BouyID", "pc", "pemb", "s_pc", "s_pemb", "MissingCI"])
resolved = prob.merge(phot, left_on=["BouyID", "SarroID"],
                      right_on=["DANCe", "__SBB2014_"], validate="one_to_one")
assert len(resolved) == len(prob) == resolved.SarroID.nunique() == 5444
resolved.to_csv(DEST / "pleiades_dance_resolved.csv", index=False)
receipt["tables"]["pleiades_dance_resolved"] = {
    "rows": len(resolved), "identity": "original Sarro numeric ID plus Bouy name",
    "unmatched": 0}

# Preserve all columns, including original literature masses, but do not use
# those masses for Ember inference. Standardize documented NaNs and the
# additional large negative UKIDSS missing-value sentinel.
source = RAW / "pleiades_hobart2026.txt"
lines = [line for line in source.read_text().splitlines()
         if line.strip() and not line.startswith("#")]
names = lines[0].split()
rows = []
missing = 0
for line in lines[1:]:
    cells = line.split()
    assert len(cells) == len(names)
    clean = []
    for value in cells:
        try:
            number = float(value)
            if not np.isfinite(number) or number < -1e8:
                value = ""
                missing += 1
        except ValueError:
            pass
        clean.append(value)
    rows.append(dict(zip(names, clean)))
write_csv(DEST / "pleiades_hobart2026.csv", names, rows)
receipt["inputs"][source.name] = hashlib.sha256(source.read_bytes()).hexdigest()
receipt["tables"]["pleiades_hobart2026"] = {"rows": len(rows), "missing_cells": missing}

source = RAW / "pleiades_boyle2025_table2_mrt.txt"
boyle = Table.read(source, format="ascii.cds")
assert len(boyle) == 3091
boyle.write(DEST / "pleiades_boyle2025.csv", format="ascii.csv", overwrite=True)
receipt["inputs"][source.name] = hashlib.sha256(source.read_bytes()).hexdigest()
receipt["tables"]["pleiades_boyle2025"] = {"rows": len(boyle)}

(DEST / "extraction_checks.json").write_text(json.dumps(receipt, indent=2) + "\n")
print(json.dumps({k: v["rows"] for k, v in receipt["tables"].items()}, indent=2))
