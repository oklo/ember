#!/usr/bin/env python3
"""Check Thoul reference units and estimate envelope settling with its coefficients.

The independent dimensional Burgers solver uses published collision equations.
Actual-gravity controls retain the measured ideal-ion pressure gradients.
Classical electron heat flow is not a physical model for a degenerate core;
omitting heat flow is a separate approximation, not an uncertainty bound.
No composition or selected stellar physics is changed.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re

import numpy as np

from diffusion_burgers import KB, MU, ME, solve
from ion_collision_integrals import E2
from prepare_nongrey_sources import digest
from write_scientific_result import write_result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("preserve previous results")
    root = Path(__file__).resolve().parents[1]
    source = root / "docs/results/thoul_diffusion_comparison_v1.json"
    data = json.loads(source.read_text())
    assert data["outcome"] == "comparison_completed" and not data["failures"]
    header = root / "include/ember/gs98_mixture.hpp"
    metals = [tuple(map(float, row.split(","))) for row in
              re.findall(r"^\s*\{([^{}]+)\}, //", header.read_text(), re.M)]
    masses = np.array([1., 4., 3.] + [v[1] for v in metals] + [ME/MU])
    charges = np.array([1., 2., 2.] + [v[0] for v in metals] + [-1.])
    paths = [Path(__file__), source, header, root/"scripts/diffusion_burgers.py",
             root/"scripts/ion_collision_integrals.py"]
    inputs = {str(p.resolve()): digest(p) for p in paths}
    records = []
    for row in data["records"]:
        checkpoint = Path(f"/tmp/ember-cn-main-sequence-v{row['version']}/checkpoint.json")
        assert digest(checkpoint) == data["input_sha256"][str(checkpoint)]
        inputs[str(checkpoint)] = digest(checkpoint)
        cp = json.loads(checkpoint.read_text())
        c = np.array(row["concentrations"])
        rho, temperature = row["density_g_cm3"], row["temperature_K"]
        n = c*rho/(MU*(masses[:-1]@c[:-1]))
        pressure = sum(n)*KB*temperature
        dc = np.array(row["concentration_gradients"])
        # He4 is the dependent concentration in the exported routine. Use
        # charge neutrality for its differential instead of a finite-log
        # approximation when checking that routine's algebra and units.
        dc[1] = 0.
        dc[1] = -(c*charges)@dc/(c[1]*charges[1])
        gp, gt = row["ideal_pressure_gradient"], row["temperature_gradient"]
        partial_gradients = n[:-1]*KB*temperature*(gp+dc[:-1]-(c@dc)/sum(c))
        reduced_mass = MU*np.outer(masses,masses)/(masses[:,None]+masses[None,:])
        # Thoul et al. equations 11 and 13, with the corrected source log.
        resistance = (4/3*np.sqrt(2*np.pi*reduced_mass)*E2**2*np.outer(n,n)
                      *np.outer(charges**2,charges**2)/(KB*temperature)**1.5
                      *np.array(row["coulomb_logarithms"]))
        resistance = (resistance+resistance.T)/2
        kwargs = dict(number_density=n, mass_u=masses, charge=charges,
                      ion_pressure_gradient=partial_gradients,
                      gravity=-pressure*gp/rho, temperature=temperature,
                      temperature_gradient=temperature*gt, resistance=resistance,
                      z=np.full(resistance.shape,.6), zprime=np.full(resistance.shape,1.3),
                      zdoubleprime=np.full(resistance.shape,2.), heat_flow="classical")
        dimensional = solve(**kwargs)
        exported = np.array(row["reference"]["ideal_pressure_gradient"]["velocity_cm_s"])
        error = float(max(abs(dimensional.velocity_cm_s-exported))/max(abs(exported)))
        # The export's rounded solar normalization and its electron-inclusive
        # mass convention differ slightly from the dimensional cgs equations.
        assert error < .01
        lo, hi = np.array(row["face_input"]).reshape(2,6)
        radius = row["native"]["radius"]
        gravity = 6.6743e-8*.5*(lo[0]+hi[0])/radius**2
        def fractions(state):
            return np.r_[state[4], .98-state[4]-state[5], state[5], [.02*v[3] for v in metals]]
        # These are the actual ideal-ion pressure gradients. Do not replace
        # total degenerate pressure by an ideal gas when obtaining gravity.
        dln_pi = (np.log(hi[2]/lo[2])+np.log(hi[3]/lo[3])
                  +np.log(fractions(hi)/fractions(lo)))/row["dr_mass_cm"]
        actual_dp = n[:-1]*KB*temperature*dln_pi
        controls = {}
        base_face = cp["numerical_history_row"]["convective_regions"][-1][0]-1
        envelope_mass = cp["model_record"]["model"][-1][0]-.5*(lo[0]+hi[0])
        for heat in ["classical", "suppressed"]:
            answer = solve(**{**kwargs, "gravity":gravity,
                              "ion_pressure_gradient":actual_dp, "heat_flow":heat})
            v = answer.velocity_cm_s
            speed = max(abs(v))
            baryon_error = float(abs((n[:-1]*masses[:-1])@v[:-1])/(rho/MU*speed))
            current_error = float(abs((n*charges)@v)/(sum(n*abs(charges))*speed))
            assert max(baryon_error,current_error,answer.scaled_backward_error) < 1e-10
            settling = {}
            if row["face"] == base_face:
                for i in range(3,len(v)-1):
                    settling[data["species"][i]] = dict(velocity_cm_s=float(v[i]),
                        fixed_reservoir_efolding_Gyr=(float(envelope_mass/(-4*np.pi*radius**2*rho*v[i])
                                                          /(31557600*1e9)) if v[i]<0 else None))
            controls[heat] = dict(velocity_cm_s=v.tolist(),
                H_velocity_over_native=float(v[0]/row["native"]["native_H_velocity"]),
                baryon_residual=baryon_error, current_residual=current_error,
                scaled_backward_error=answer.scaled_backward_error,
                base_metal_settling=settling)
        records.append(dict(version=row["version"], age_years=row["age_years"], face=row["face"],
            q=row["q"], eta=row["native"]["eta"], convective_face=row["convective_face"],
            envelope_base=row["face"]==base_face, temperature_K=temperature, density_g_cm3=rho,
            gravity_cm_s2=gravity, envelope_mass_g=envelope_mass if row["face"]==base_face else None,
            dimensional_vs_exported_velocity_error=error,
            formal_ideal_over_actual_gravity=kwargs["gravity"]/gravity,
            native_H_velocity_cm_s=row["native"]["native_H_velocity"],
            native_over_formal_reference_H=row["native_over_classical_reference_H"],
            actual_gravity_controls=controls))
    assert all(digest(p)==value for p,value in inputs.items())
    result = dict(utc=datetime.now(timezone.utc).isoformat(),outcome="assessment_completed",
        selected=False, scope=__doc__, records=records, species=data["species"],
        maximum_dimensional_vs_exported_velocity_error=max(r["dimensional_vs_exported_velocity_error"] for r in records),
        input_sha256=inputs,
        references=["https://arxiv.org/pdf/astro-ph/9304005",
                    "https://www.sns.ias.edu/~jnb/SNdata/Export/Diffusion/README",
                    "https://arxiv.org/abs/1710.08424"],
        limitations=["All ions are fully stripped; heavy-element ionization needs a separate assessment.",
            "Formal export comparison tests units and classical equations; it is not valid degenerate-core gravity.",
            "Actual-gravity controls use measured ideal-ion gradients and omit Coulomb chemical forces.",
            "Classical residual electron heat flow is invalid in a strongly degenerate core; suppressed heat is a separate approximation, not a bound.",
            "Fixed-metal initial profiles are used; actual composition gradients and the moving convective boundary can change subsequent depletion.",
            "No physical prescription has been selected or stellar checkpoint changed."])
    write_result(args.output,result)
    print(json.dumps(dict(outcome=result["outcome"], rows=len(records),
        maximum_dimensional_error=result["maximum_dimensional_vs_exported_velocity_error"],
        base_results=[r for r in records if r["envelope_base"]]),indent=2))


if __name__ == "__main__":
    main()
