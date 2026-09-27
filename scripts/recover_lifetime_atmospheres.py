#!/usr/bin/env python3
"""Recover retained atmosphere values from their numerical source reports.

This performs no atmosphere calculation or extrapolation. Missing cells remain
masked. Each recovered number is traced to a saved, previously accepted source.
"""
import hashlib
import itertools
import json
import math
from pathlib import Path

from audit_atmosphere_retention import read


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    root = Path(__file__).resolve().parents[1]
    output = root / "data/atmosphere/lifetime"
    output.mkdir(parents=True, exist_ok=True)
    provenance = {}

    def report(name):
        path = root / "docs/results" / (name + ".json")
        provenance[str(path.relative_to(root))] = digest(path)
        return json.loads(path.read_text())

    source = root / "data/atmosphere/nongrey_gs98_z020_x085_t4000_v1.dat"
    provenance[str(source.relative_to(root))] = digest(source)
    header, axes, states = read(source)
    preserved = dict(states)
    additions = []
    for name in ("warm_metal_atmosphere_reference_v2", "warm_metal_atmosphere_reference_v3"):
        for row in report(name)["added_actual_states"]:
            h, y, t, g = row["coordinates"]
            raw = (h, y, math.log10(t), g)
            # Match only roundoff in coordinates of the existing source grid.
            key = tuple(min(axis, key=lambda x: abs(x - q)) for axis, q in zip(axes, raw))
            if any(abs(a - b) > 1e-12 for a, b in zip(key, raw)):
                raise ValueError("saved source is not on an existing grid coordinate")
            state = row["state"]
            value = (1, math.log10(state["T"]), math.log10(state["Pgas"]))
            if states[key][0] and states[key] != value:
                raise ValueError("refusing to replace a retained atmosphere value")
            states[key] = value
            additions.append(dict(coordinates=key, report=name, work=row["work"]))
    lines = source.read_text().splitlines()
    rows = lines[:lines.index("data") + 1]
    rows += [" ".join(format(v, ".17g") for v in states[key]) for key in itertools.product(*axes)]
    (output / "main_sequence_reference.dat").write_text("\n".join(rows) + "\n")
    assert all(states[key] == value for key, value in preserved.items() if value[0])

    recovered = []
    for ref_z, source_z, name, filename in (
        (.005, .0025, "depleted_metal_atmosphere_response_v2", "metal_z005_z0025.dat"),
        (.0025, .001, "z001_atmosphere_response_v1", "metal_z0025_z001.dat"),
        (.001, .0005, "z0005_atmosphere_response_v1", "metal_z001_z0005.dat"),
    ):
        pairs = report(name)["comparisons"]
        cells = {}
        for row in pairs:
            h, t, g = row["coordinates"]
            key = (h, math.log10(t), g)
            value = (row["delta_logT"], row["delta_logPgas"])
            # The stored responses must agree with the independently recorded
            # matching states, including the sign of the lower-Z response.
            expected = (math.log(row["source_T"] / row["reference_T"]),
                        math.log(row["source_Pgas"] / row["reference_Pgas"]))
            if any(abs(a - b) > 2e-14 for a, b in zip(value, expected)):
                raise ValueError("saved response disagrees with its source states")
            if key in cells:
                raise ValueError("duplicate recovered response coordinate")
            cells[key] = value
        response_axes = [sorted({key[k] for key in cells}) for k in range(3)]
        lines = ["EMBER_METAL_ATMOSPHERE_RESPONSE 1",
                 'source "Measured LTE gas-atmosphere pairs retained in ' + name + '.json"',
                 'approximation "Linear-in-Z logarithmic T/Pgas response; separable He3 <= .003; fixed GS98 ratios"',
                 "basis baryon_mass", f"reference_Z {ref_z:.17g}", f"source_Z {source_z:.17g}",
                 "maximum_helium3 .003", "tau 100"]
        for label, axis in zip(("hydrogen", "log_teff", "log_g"), response_axes):
            lines.append(label + " " + str(len(axis)) + " " + " ".join(f"{v:.17g}" for v in axis))
        lines.append("data")
        for key in itertools.product(*response_axes):
            lines.append("1 " + " ".join(f"{v:.17g}" for v in cells[key]) if key in cells else "0")
        (output / filename).write_text("\n".join(lines) + "\n")
        recovered.append(dict(file=filename, reference_Z=ref_z, source_Z=source_z,
                              paired_sources=len(cells), original_report=name))
    chain = [(".02", "metal_z020_z010.dat"), (".01", "metal_z010_z005.dat")]
    chain += [(format(row["reference_Z"], ".17g"), row["file"]) for row in recovered]
    (output / "metal_chain_extended.dat").write_text("EMBER_METAL_ATMOSPHERE_CHAIN 1 5\n" +
        "".join(z + " " + json.dumps(name) + "\n" for z, name in chain))
    receipt = dict(input_sha256=provenance, unchanged_present_states=sum(v[0] for v in preserved.values()),
                   recovered_reference_states=additions, recovered_responses=recovered,
                   limitations=["No new atmosphere solution or independent physical validation.",
                                "Gas only, tau=100, fixed GS98 ratios, separable small-He3 response.",
                                "The reference H axis ends at .85; later H-rich families remain separate.",
                                "Missing source values remain masked; derivative changes at newly supported edges need checking."],
                   output_sha256={p.name: digest(p) for p in output.glob("*.dat")})
    (output / "recovery.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"unchanged_present_states": receipt["unchanged_present_states"],
                      "recovered_reference_states": len(additions), "responses": recovered}))


if __name__ == "__main__":
    main()
