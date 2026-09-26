#!/usr/bin/env python3
"""Compare retained metal-depletion spectra at identical temperature and density.

This performs no new source calculations and changes no installed material
table. The means include source absorption only, excluding scattering/grains.
"""
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path

import numpy as np

from nongrey_opacity import read_table, validate_table
from prepare_nongrey_sources import digest
from write_scientific_result import write_result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("preserve previous result")
    work = Path("/tmp/ember-metal-depletion-atmosphere-v1")
    reference = Path("/tmp/ember-primary-nongrey-x085-opacity-v1/plane-000/opacity/fort.63")
    inputs = {str(Path(__file__).resolve()): digest(__file__)}
    records = {}
    base = None
    for label, path in [("reference", reference),
                        ("cn_only", work / "cn_only/opacity/fort.63"),
                        ("all_metals", work / "all_metals/opacity/fort.63")]:
        table = read_table(path)
        inputs[str(path.resolve())] = digest(path)
        temperatures = np.exp(table["log_temperature"])
        densities = np.exp(table["log_density"])
        if base is None:
            base = table
            refcheck = Path("/tmp/ember-primary-nongrey-x085-warm-v1/x850-y000-t3800-g490/provenance.json")
            assert digest(path) == json.loads(refcheck.read_text())["opacity_sha256"]
            inputs[str(refcheck)] = digest(refcheck)
        else:
            for key in ["shape", "frequency", "log_temperature", "log_density", "flags"]:
                assert table[key] == base[key], key
            control_path = work / label / "composition.json"
            receipt = work / label / "opacity_complete.json"
            control = json.loads(control_path.read_text())
            assert digest(path) == json.loads(receipt.read_text())["table_sha256"]
            validate_table(path, control["number_abundances_relative_to_H"], temperatures, densities)
            for p in [control_path, receipt]:
                inputs[str(p)] = digest(p)
        frequencies = np.asarray(table["frequency"])[::-1]
        log_opacity = np.asarray(table["log_opacity"], dtype=float).reshape(table["shape"])[::-1]
        means, fractions = [], []
        for index, temperature in enumerate(temperatures):
            # Same physical constants as the retained atmosphere source.
            x = 6.6256e-27 * frequencies / (1.38054e-16 * temperature)
            weight = x**4 * np.exp(-x) / (-np.expm1(-x))**2
            norm = np.trapezoid(weight, x)
            mean = norm / np.trapezoid(weight[:, None] * np.exp(-log_opacity[:, :, index]), x, axis=0)
            assert np.isfinite(mean).all() and np.min(mean) > 0
            means.append(mean.tolist())
            fractions.append(float(norm / (4 * math.pi**4 / 15)))
        # Source electron-density storage order: temperature, density.
        records[label] = dict(absorption_mean_cm2_g=means,
                              log_electron_density=np.asarray(table["log_electron_density"]).reshape(len(temperatures), len(densities)).tolist(),
                              rosseland_weight_fraction=fractions)
    reference_mean = np.array(records["reference"]["absorption_mean_cm2_g"])
    reference_electrons = np.array(records["reference"]["log_electron_density"])
    for record in records.values():
        record["opacity_ratio_to_reference"] = (np.array(record["absorption_mean_cm2_g"]) / reference_mean).tolist()
        record["electron_density_ratio_to_reference"] = np.exp(np.array(record["log_electron_density"]) - reference_electrons).tolist()
    samples = []
    for target_t in [3500, 5000, 6500, 10000]:
        ti = int(np.argmin(abs(np.log(temperatures / target_t))))
        for target_rho in [1e-7, 1e-5]:
            ri = int(np.argmin(abs(np.log(densities / target_rho))))
            samples.append(dict(temperature_K=float(temperatures[ti]), density_g_cm3=float(densities[ri]),
                                values={label: {key: record[key][ti][ri] for key in
                                        ["absorption_mean_cm2_g", "opacity_ratio_to_reference", "electron_density_ratio_to_reference"]}
                                        for label, record in records.items()}))
    assert all(digest(path) == value for path, value in inputs.items())
    result = dict(utc=datetime.now(timezone.utc).isoformat(), outcome="comparison_completed",
                  selected=False, scope=__doc__, new_source_jobs=0,
                  temperature_K=temperatures.tolist(), density_g_cm3=densities.tolist(),
                  records=records, representative_source_nodes=samples, input_sha256=inputs,
                  limitations=["Source absorption mean; scattering and grain absorption are excluded.",
                               "At fixed hydrogen X=.85, helium3=0. Removed metal mass is replaced by helium4.",
                               "No stellar effective temperature or settling rate follows from these material-state ratios alone.",
                               "A ratio change does not identify the responsible molecule without individual absorption contributions."])
    write_result(args.output, result)
    print(json.dumps(dict(outcome=result["outcome"], representative_source_nodes=samples), indent=2))


if __name__ == "__main__":
    main()
