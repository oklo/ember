#!/usr/bin/env python3
"""Estimate trace-metal removal scales at the current convective-envelope base.

Uses the hydrogen-background force in Heinonen et al. (2020), equations 22
and 17, and the retained screened-ion collision integrals. This is an
order-of-magnitude diagnostic, not independent validation of those integrals
or a multicomponent stellar transport solution.
"""
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re

import numpy as np
from scipy.interpolate import CubicSpline

from ion_collision_integrals import KB, MU, E2, screening_length
from prepare_nongrey_sources import digest
from write_scientific_result import write_result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("preserve previous result")
    root = Path(__file__).resolve().parents[1]
    checkpoint = Path("/tmp/ember-cn-main-sequence-v22/checkpoint.json")
    source = json.loads(checkpoint.read_text())
    model = np.array(source["model_record"]["model"])
    index = source["numerical_history_row"]["convective_regions"][-1][0]
    face = (model[index - 1] + model[index]) / 2
    mass, radius, density, temperature = face[:4]
    envelope_mass = model[-1, 0] - mass
    gravity = 6.6743e-8 * mass / radius**2
    header = root / "include/ember/gs98_mixture.hpp"
    metals = [tuple(map(float, row.split(","))) for row in
              re.findall(r"^\s*\{([^{}]+)\}, //", header.read_text(), re.M)]
    assert len(metals) == 19 and abs(sum(row[3] for row in metals) - 1) < 1e-12
    masses = np.array([1., 3., 4.] + [row[1] for row in metals])
    charges = np.array([1., 2., 2.] + [row[0] for row in metals])
    fractions = np.r_[face[5:7], 1 - face[5] - face[6] - .02,
                      [.02 * row[3] for row in metals]]
    assert min(fractions) >= 0 and abs(sum(fractions) - 1) < 1e-12
    number = density / MU * fractions / masses
    screen = screening_length(number, charges, temperature)
    report_path = root / "docs/results/yukawa_collision_table_v1.json"
    receipt = json.loads(report_path.read_text())
    table_path = Path(receipt["table"])
    assert receipt["status"] == "pass" and digest(table_path) == receipt["table_sha256"]
    table = json.loads(table_path.read_text())
    moments = CubicSpline(table["log10_strength"], np.log(table["dimensionless_integrals"]),
                          axis=0, extrapolate=False)
    inputs = {str(p.resolve()): digest(p) for p in
              [Path(__file__), checkpoint, header, report_path, table_path,
               root / "scripts/ion_collision_integrals.py"]}
    # Ideal ionic pressure on each side of the base. No electron pressure is
    # added to this term in the published trace-force expression.
    def ionic_pressure(row):
        ion_number_per_mass = (row[5] + row[6] / 3 + (1-row[5]-row[6]-.02) / 4
                               + .02 * sum(v[3] / v[1] for v in metals)) / MU
        return row[2] * KB * row[3] * ion_number_per_mass
    gradients = []
    for width in [1, 2]:
        left, right = model[index-width], model[index+width-1]
        gradient = math.log(ionic_pressure(right)/ionic_pressure(left)) / (right[1]-left[1])
        assert gradient < 0
        gradients.append(dict(stencil_width=width, dln_ionic_pressure_dr_cm_inverse=gradient))
    records = []
    for element, a, z in [("C",12,6), ("N",14,7), ("O",16,8), ("Na",23,11),
                          ("Mg",24,12), ("Si",28,14), ("Ca",40,20), ("Fe",56,26)]:
        for charge_fraction in [1., .5]:
            charge = z * charge_fraction
            # Test-particle friction against H, He3 and He4. Keeping the
            # abundant helium drag improves on a strictly pure-H resistance.
            reduced_mass = MU * a * masses[:3] / (a+masses[:3])
            interaction = charge * charges[:3] * E2
            strength = interaction / (KB * temperature * screen["length_cm"])
            k11 = np.exp(moments(np.log10(strength)))[:, 0]
            assert np.isfinite(k11).all() and min(k11) > 0
            omega = np.sqrt(2*math.pi/reduced_mass)*interaction**2 / (KB*temperature)**1.5 * k11
            diffusivity = 3 * KB * temperature / (16 * np.sum(number[:3] * reduced_mass * omega))
            for gradient in gradients:
                force = ((charge-a)*MU*gravity/(KB*temperature)
                         + (charge-1)*gradient["dln_ionic_pressure_dr_cm_inverse"])
                velocity = diffusivity * force
                assert velocity < 0
                timescale = envelope_mass / (-4*math.pi*radius**2*density*velocity)
                records.append(dict(element=element, mass_number=a, charge=charge,
                                    charge_fraction=charge_fraction, **gradient,
                                    trace_diffusivity_cm2_s=float(diffusivity),
                                    force_over_kT_cm_inverse=float(force),
                                    velocity_cm_s=float(velocity),
                                    fixed_reservoir_efolding_Gyr=float(timescale/(31557600*1e9)),
                                    collision_strength=strength.tolist()))
    assert all(digest(p) == value for p, value in inputs.items())
    result = dict(utc=datetime.now(timezone.utc).isoformat(), outcome="diagnostic_completed",
                  selected=False, scope=__doc__,
                  force_reference="https://arxiv.org/pdf/2005.05891, equations 17 and 22",
                  base=dict(face_index=index, temperature_K=temperature, density_g_cm3=density,
                            radius_cm=radius, gravity_cm_s2=gravity, envelope_mass_g=envelope_mass,
                            hydrogen_mass_fraction=face[5], hydrogen_ion_number_fraction=float(number[0]/sum(number)),
                            screening=screen), records=records, input_sha256=inputs,
                  limitations=["Hydrogen-background force with measured ionic-pressure gradient; H/He multicomponent electric field is approximated.",
                               "Trace mobility uses retained screened-ion integrals, not independent transport coefficients.",
                               "Half and full ionic charge are illustrative sensitivities, not ionization predictions or uncertainty bounds.",
                               "Thermal diffusion, concentration gradients, electron drag and Coulomb chemical forces are omitted.",
                               "Fixed reservoir estimate; changing convective mass and replenishment can alter the integrated depletion.",
                               "No abundance, opacity, atmosphere or stellar trajectory is changed."])
    write_result(args.output, result)
    print(json.dumps(dict(base=result["base"], samples=[r for r in records if r["stencil_width"]==1]), indent=2))


if __name__ == "__main__":
    main()
