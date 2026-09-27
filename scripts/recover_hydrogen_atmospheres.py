#!/usr/bin/env python3
"""Restore measured hydrogen-rich atmospheres from retained source receipts.

No new source column, extrapolation, or missing-cell fill is introduced. The
output is separate from every running calculation's immutable input package.
"""
import hashlib
import itertools
import json
import math
from pathlib import Path
import shlex

from audit_atmosphere_retention import read


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/atmosphere/lifetime_hydrogen"
PROVENANCE = {}
COUNTS = {}


def source(path):
    path = ROOT / path
    PROVENANCE[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path


def report(name):
    return json.loads(source("docs/results/" + name + ".json").read_text())


def put(cells, key, value):
    key = tuple(round(float(x), 11) for x in key)
    if key in cells:
        if max(abs(a - b) for a, b in zip(cells[key], value)) > 5e-12:
            raise ValueError(f"conflicting measured state at {key}: {cells[key]} / {value}")
    else:
        cells[key] = tuple(value)


def grid(name, cells, labels, header, logarithmic_temperature_axis):
    axes = [sorted({key[i] for key in cells}) for i in range(len(labels))]
    lines = list(header)
    for i, (label, axis) in enumerate(zip(labels, axes)):
        values = [math.log10(x) for x in axis] if i == logarithmic_temperature_axis else axis
        lines.append(label + " " + str(len(axis)) + " " + " ".join(format(x, ".17g") for x in values))
    lines.append("data")
    for key in itertools.product(*axes):
        lines.append("1 " + " ".join(format(v, ".17g") for v in cells[key]) if key in cells else "0")
    (OUT / name).write_text("\n".join(lines) + "\n")
    COUNTS[name] = dict(measured=len(cells), missing=math.prod(map(len, axes))-len(cells), axes=axes)


def composition(name, z, cells, metals):
    grid(name, cells, ("hydrogen", "helium3", "log_teff", "log_g"), [
        "EMBER_COMPOSITION_ATMOSPHERE 2",
        'source "Retained TLUSTY/SYNSPEC matching states; provenance in recovery.json"',
        'approximation "LTE gas, GS98 metal ratios, effective isotope masses; no grains or irradiation"',
        "basis baryon_mass", "tau 100",
        "metals " + " ".join(format(x*z/.02, ".17g") for x in metals)], 2)


def response(name, za, zb, cells):
    grid(name, cells, ("hydrogen", "log_teff", "log_g"), [
        "EMBER_METAL_ATMOSPHERE_RESPONSE 1",
        'source "Retained paired gas-atmosphere matching states; provenance in recovery.json"',
        'approximation "Linear-in-Z logarithmic T/Pgas response, separable He3 <= .003, GS98 ratios"',
        "basis baryon_mass", f"reference_Z {za:.17g}", f"source_Z {zb:.17g}",
        "maximum_helium3 .003", "tau 100"], 1)


def chain(name, rows):
    (OUT / name).write_text(f"EMBER_METAL_ATMOSPHERE_CHAIN 1 {len(rows)}\n" +
                           "".join(f"{z:.17g} {json.dumps(p)}\n" for z, p in rows))


def measured(rows, z=None):
    for r in rows:
        c = r["coordinates"]
        if len(c) == 5:
            zz, h, y, t, g = c
        else:
            h, y, t, g = c
            zz = r.get("Z", z)
        if zz is None:
            raise ValueError("unidentified source metallicity")
        s = r["state"]
        if "tau_bracket" in s and not s["tau_bracket"][0] <= 100 <= s["tau_bracket"][1]:
            raise ValueError("source does not bracket the matching depth")
        yield zz, (h, y, t, g), (math.log10(s["T"]), math.log10(s["Pgas"]))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    base = source("data/atmosphere/nongrey_gs98_z020_x090_late_v1.dat")
    header, _, old = read(base)
    metals = list(map(float, header["metals"]))
    high = {}
    # The retained warm family begins at 4000 K; earlier temperatures remain
    # in the broad main-sequence reference used below the hydrogen join.
    for (h, y, lt, g), s in old.items():
        if s[0] and abs(10**lt - 4000) < 1e-6:
            put(high, (h, y, 4000, g), s[1:])
    _, _, warm = read(source("data/atmosphere/lifetime/main_sequence_reference.dat"))
    for (h, y, lt, g), s in warm.items():
        t = round(10**lt)
        if s[0] and h in (.7, .85) and y == 0 and t in (4200, 4400) and g in (4.9, 5.15):
            put(high, (h, y, t, g), s[1:])
    for name, key in (("x090_warm_atmosphere_reference_v1", "new_source_states"),
                      ("x090_warm_atmosphere_reference_v2", "new_source_states"),
                      ("h094_warm_atmosphere_reference_v1", "new_source_states"),
                      ("highg_material_assembly_v1", "source_states")):
        for z, c, s in measured(report(name)[key], .02):
            if z == .02:
                put(high, c, s)
    if len(high) != 72:
        raise ValueError(f"expected 72 retained high-Z reference states, got {len(high)}")
    composition("reference_z020.dat", .02, high, metals)

    high_chain = []
    for za, zb, name in ((.02, .01, "metal_z020_z010.dat"), (.01, .005, "metal_z010_z005.dat")):
        (OUT / name).write_bytes(source("data/atmosphere/lifetime/" + name).read_bytes())
        high_chain.append((za, name))
    for za, zb, name, suffix in ((.005, .0025, "metal_z005_z0025.dat", "z0025_response.dat"),
                                (.0025, .001, "metal_z0025_z001.dat", "z001_response.dat"),
                                (.001, .0005, "metal_z001_z0005.dat", "z0005_response.dat")):
        p = source("data/atmosphere/lifetime/" + name)
        lines = p.read_text().splitlines()
        axes = [list(map(float, shlex.split(next(x for x in lines if x.startswith(label+' ')))[2:]))
                for label in ("hydrogen", "log_teff", "log_g")]
        cells = {}
        for c, line in zip(itertools.product(*axes), lines[lines.index("data")+1:], strict=True):
            row = list(map(float, line.split()))
            if row[0]:
                put(cells, (c[0], 10**c[1], c[2]), row[1:])
        for receipt in ("h094_metal_atmosphere_assembly_v1", "highg_material_assembly_v1"):
            table = next(x for x in report(receipt)["tables"] if x["table"].endswith('/'+suffix))
            for row in table["added"]:
                if row["reference_Z"] != za or row["source_Z"] != zb:
                    raise ValueError("response metallicity mismatch")
                put(cells, row["coordinates"], row["delta"])
        response(name, za, zb, cells)
        high_chain.append((za, name))
    chain("chain_z020.dat", high_chain)

    # The low-Z reference is independent: keeping Z=.02 here would require
    # a negative helium abundance at the largest hydrogen fractions.
    reports = [report(n) for n in ("low_z_material_assembly_v1", "z0001_material_assembly_v1",
                "low_z_warmer_material_assembly_v1", "h0995_material_assembly_v1")]
    states = {}
    for r in reports:
        for z, c, s in measured(r["source_states"]):
            put(states, (z, *c), s)
    low = {c[1:]: s for c, s in states.items() if c[0] == .001}
    if len(low) != 82:
        raise ValueError(f"expected 82 retained low-Z states, got {len(low)}")
    composition("reference_z001.dat", .001, low, metals)
    low_chain = []
    for za, zb, name, count in ((.001, .0005, "response_z001_z0005.dat", 32),
            (.0005, .0001, "response_z0005_z0001.dat", 26),
            (.0001, .00001, "response_z0001_z00001.dat", 8)):
        cells = {}
        # Original declared pairs choose their own temperature and gravity
        # axes. Do not insert extra source knots into an incomplete rectangle.
        if za == .001:
            pairs = reports[0]["actual_response_pairs"]
        elif za == .0005:
            pairs = reports[1]["actual_response_pairs"]
        else:
            pairs = reports[3]["lowest_metal_response_pairs"]
        for row in pairs:
            put(cells, row["coordinates"], row["delta"])
        if za != .0001:
            for r in reports[2:]:
                for z, c, _ in measured(r["source_states"]):
                    h, y, t, g = c
                    if z != za or y != 0:
                        continue
                    akey = tuple(round(x,11) for x in (za,h,y,t,g))
                    bkey = tuple(round(x,11) for x in (zb,h,y,t,g))
                    if bkey in states:
                        a, b = states[akey], states[bkey]
                        put(cells, (h,t,g), tuple((v-u)*math.log(10) for u,v in zip(a,b)))
        if len(cells) != count:
            raise ValueError(f"unexpected measured response count {name}: {len(cells)} != {count}")
        response(name, za, zb, cells)
        low_chain.append((za, name))
    chain("chain_z001.dat", low_chain)
    chain("lower_z0005.dat", low_chain[1:2])
    chain("lower_z0001.dat", low_chain[2:])

    def retained(checksum, name):
        path = source(f"out/material-input-recovery-v1/closed/{checksum[:2]}/{checksum}.dat")
        if hashlib.sha256(path.read_bytes()).hexdigest() != checksum:
            raise ValueError("retained table digest mismatch")
        (OUT / name).write_bytes(path.read_bytes())
    retained("2164f572768684c8197f5d624efd2bb5fa07ae7aa2d0af74a6637809579e287c", "reference_z0001.dat")
    retained("98fd271a90878910ad2f98c5bc42d40401668523618197b8a7841205c31b34f0", "near_response_z0001_z00001.dat")
    retained("781f56be194ab1f89b43b5ec89d2c2b6b341336822d60f93f695c60369c69057", "near_response_z00001_z000001.dat")
    chain("chain_z0001.dat", [(.0001,"near_response_z0001_z00001.dat"), (.00001,"near_response_z00001_z000001.dat")])
    (OUT / "interval_z020.dat").write_text('EMBER_HYDROGEN_ATMOSPHERE_INTERVAL 1\nreference "reference_z020.dat"\nchain "chain_z020.dat"\nreference_Z .02\nhydrogen .85 .9\n')
    (OUT / "interval_z001.dat").write_text('EMBER_HYDROGEN_ATMOSPHERE_INTERVAL 2\nlower_interval "interval_z020.dat"\nlower_metal_chain "lower_z0005.dat"\nlower_reference_Z .0005\nreference "reference_z001.dat"\nchain "chain_z001.dat"\nreference_Z .001\nhydrogen .9 .94\n')
    (OUT / "interval_z0001.dat").write_text('EMBER_HYDROGEN_ATMOSPHERE_INTERVAL 2\nlower_interval "interval_z001.dat"\nlower_metal_chain "lower_z0001.dat"\nlower_reference_Z .0001\nreference "reference_z0001.dat"\nchain "chain_z0001.dat"\nreference_Z .0001\nhydrogen .99 .995\n')
    result = dict(input_sha256=PROVENANCE, tables=COUNTS,
                  limitations=["Recovery of existing accepted sources, not new physical validation.",
                    "No selected stellar trajectory yet; native corridor and derivative checks required.",
                    "Gas only at tau 100, fixed GS98 ratios, bounded separable He3 response.",
                    "Highest hydrogen fraction .9995; later trace-He and pure-H families are separate."],
                  output_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.glob('*.dat')})
    (OUT / "recovery.json").write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(COUNTS))


if __name__ == "__main__":
    main()
