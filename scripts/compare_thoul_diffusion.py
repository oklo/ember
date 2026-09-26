#!/usr/bin/env python3
"""Compare retained stellar faces with the exported Thoul diffusion routine.

The classical reference includes freely moving metals. Ember uses degenerate
electron transport, nonideal material forces, and stationary metal collision
carriers. Ideal-force controls help distinguish thermodynamic differences.
This calculation does not select a new transport model or evolve a star.
"""
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
import subprocess

import numpy as np

from ion_collision_integrals import KB, MU, ME
from prepare_nongrey_sources import digest
from write_scientific_result import write_result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("preserve completed results")
    root = Path(__file__).resolve().parents[1]
    stage = Path("/tmp/ember-thoul-reference-v1")
    eos = Path("/tmp/ember-hydrogen-x090-eos-family-v1/freeeos300_gs98_z020.dat")
    collision = Path("/tmp/ember-collision-eta32-lowb-192-v1.dat")
    header = root / "include/ember/gs98_mixture.hpp"
    metals = [tuple(map(float, row.split(","))) for row in
              re.findall(r"^\s*\{([^{}]+)\}, //", header.read_text(), re.M)]
    names = ["H1", "He4", "He3"] + re.findall(r"}, // (\w+)", header.read_text()) + ["electron"]
    masses = np.array([1., 4., 3.] + [v[1] for v in metals] + [ME/MU])
    charges = np.array([1., 2., 2.] + [v[0] for v in metals] + [-1.])
    count = len(masses)
    assert count == len(names) == 23
    paths = [Path(__file__), header, eos, collision, stage/"source.json", stage/"build.json",
             stage/"native_build.json", stage/"native_face_probe", stage/"probe_real4", stage/"probe_real8"]
    inputs = {str(p.resolve()): digest(p) for p in paths}
    for line in eos.read_text().splitlines()[3:]:
        p = (eos.parent/json.loads(line)).resolve()
        inputs[str(p)] = digest(p)
    samples = []
    def abundance(row):
        return np.r_[row[5], .98-row[5]-row[6], row[6], [.02*v[3] for v in metals], 0.]
    def concentration(fractions):
        zxa = np.sum(charges[:-1]*fractions[:-1]/masses[:-1])
        return np.r_[fractions[:-1]/masses[:-1]/zxa, 1.]
    for version in [6, 10, 19, 22]:
        cp = Path(f"/tmp/ember-cn-main-sequence-v{version}/checkpoint.json")
        d = json.loads(cp.read_text()); inputs[str(cp)] = digest(cp)
        model = np.array(d["model_record"]["model"])
        regions = d["numerical_history_row"]["convective_regions"]
        q = .5*(model[:-1,0]+model[1:,0])/model[-1,0]
        indices = sorted({int(np.argmin(abs(q-target))) for target in [.0001,.01,.1,.4,.7,.9]} |
                         {regions[-1][0]-1})
        for index in indices:
            lo,hi = model[index:index+2]
            if min(lo[3],hi[3]) < 2e6:
                continue
            fractions = (abundance(lo)+abundance(hi))/2
            c = concentration(fractions)
            temperature = math.sqrt(lo[3]*hi[3]); density = (lo[2]+hi[2])/2
            # Source README constants and its corrected 1.2 exponent are
            # retained exactly. No floor is imposed on the Coulomb logarithm.
            ne = density/(1.6726e-24*float(masses@c))
            ni = ne*sum(c[:-1]); a0 = (0.23873/ni)**(1/3)
            debye = 6.9010*math.sqrt(temperature/(ne*float(c@(charges**2))))
            screening = max(debye,a0)
            xij = 2.3939e3*temperature*screening/abs(np.outer(charges,charges))
            cl = .81245*np.log1p(.18769*xij**1.2)
            radius = .5*(lo[1]+hi[1]); dr = (hi[0]-lo[0])/(4*np.pi*radius**2*density)
            def ideal_pressure(row):
                x=abundance(row)
                return row[2]/MU*KB*row[3]*np.sum(x[:-1]*(1+charges[:-1])/masses[:-1])
            pressure = .5*(ideal_pressure(lo)+ideal_pressure(hi))
            gravity = 6.6743e-8*.5*(lo[0]+hi[0])/radius**2
            samples.append(dict(version=version, age_years=d["age_seconds"]/31557600,
                                face=index, q=float(q[index]), temperature_K=temperature,
                                density_g_cm3=density, mass_fractions=fractions.tolist(),
                                concentrations=c.tolist(), coulomb_logarithms=cl.tolist(),
                                source_screening_length_cm=screening,
                                dr_mass_cm=dr, ideal_pressure_gradient=math.log(ideal_pressure(hi)/ideal_pressure(lo))/dr,
                                hydrostatic_ideal_pressure_gradient=-density*gravity/pressure,
                                temperature_gradient=math.log(hi[3]/lo[3])/dr,
                                concentration_gradients=(np.log(concentration(abundance(hi))/concentration(abundance(lo)))/dr).tolist(),
                                convective_face=any(a<=index and index+1<b for a,b in regions),
                                face_input=[float(row[k]) for row in [lo,hi] for k in [0,1,2,3,5,6]]))
    native_input="".join(" ".join(format(v,".17g") for v in r["face_input"])+"\n" for r in samples)
    p=subprocess.run([str(stage/"native_face_probe"),str(eos),str(collision)],input=native_input,
                     text=True,capture_output=True,timeout=120,check=True)
    native=[json.loads(line) for line in p.stdout.splitlines()]
    (stage/"face_inputs.txt").write_text(native_input)
    (stage/"native_responses.jsonl").write_text(p.stdout)
    assert len(native)==len(samples)
    query=[]
    for row in samples:
        query.append(str(count))
        for v in [masses,charges,row["mass_fractions"],*row["coulomb_logarithms"]]:
            query.append(" ".join(format(x,".17g") for x in v))
    text="\n".join(query)+"\n"; (stage/"reference_inputs.txt").write_text(text)
    answers={}
    for precision in ["real4","real8"]:
        p=subprocess.run([str(stage/("probe_"+precision))],input=text,text=True,capture_output=True,
                         timeout=30,check=True)
        (stage/(precision+"_responses.txt")).write_text(p.stdout)
        values=np.array([list(map(float,line.split())) for line in p.stdout.splitlines()])
        assert values.shape==(len(samples)*count,count+2) and np.isfinite(values).all()
        answers[precision]=values.reshape(len(samples),count,count+2)
    precision_error=[]; records=[]; failures=[]
    for index,(sample,response) in enumerate(zip(samples,native,strict=True)):
        if "error" in response:
            failures.append(dict(version=sample["version"],face=sample["face"],error=response["error"]))
            continue
        a4,a8=answers["real4"][index],answers["real8"][index]
        precision_error.append(float(np.max(abs(a4[:3]-a8[:3]))/np.max(abs(a8[:3]))))
        reference={}
        # The source defines gradients per solar radius and speeds in
        # R_sun/(6e13 yr). Keep its published normalization explicit.
        scale=(sample["temperature_K"]/1e7)**2.5/(sample["density_g_cm3"]/100)*6.957e10**2/(6e13*31557600)
        for gradient in ["ideal_pressure_gradient","hydrostatic_ideal_pressure_gradient"]:
            pressure_term=a8[:,0]*sample[gradient]
            thermal_term=a8[:,1]*sample["temperature_gradient"]
            concentration_term=a8[:,2:]@np.array(sample["concentration_gradients"])
            velocity=(pressure_term+thermal_term+concentration_term)*scale
            c=np.array(sample["concentrations"]); speed=np.max(abs(velocity))
            mass_residual=abs((c*masses)@velocity)/(sum(c*masses)*speed)
            current_residual=abs((c*charges)@velocity)/(sum(c*abs(charges))*speed)
            assert max(mass_residual,current_residual)<1e-8
            reference[gradient]=dict(velocity_cm_s=velocity.tolist(),hydrogen_pressure_cm_s=float(pressure_term[0]*scale),
                                     hydrogen_thermal_cm_s=float(thermal_term[0]*scale),
                                     hydrogen_concentration_cm_s=float(concentration_term[0]*scale),
                                     mass_residual=float(mass_residual),current_residual=float(current_residual))
        recompute=np.array(response["recomputed_rate"]);actual=np.array(response["native_rate"])
        recompute_error=float(np.max(abs(recompute-actual))/max(np.max(abs(actual)),1e-100))
        assert recompute_error<1e-5
        v=reference["ideal_pressure_gradient"]["velocity_cm_s"][0]
        records.append(dict(**sample,native=response,reference=reference,
                            native_recomputed_rate_error=recompute_error,
                            native_over_classical_reference_H=response["native_H_velocity"]/v,
                            ideal_force_over_classical_reference_H=response["ideal_force_H_velocity"]/v,
                            reference_coefficients=a8.tolist()))
    assert all(digest(p)==value for p,value in inputs.items())
    result=dict(utc=datetime.now(timezone.utc).isoformat(),outcome="comparison_completed" if not failures else "incomplete_comparison",
                selected=False,scope=__doc__,species=names,records=records,failures=failures,
                maximum_real4_real8_HHe_coefficient_difference=max(precision_error),input_sha256=inputs,
                reference_source="https://www.sns.ias.edu/~jnb/SNdata/Export/Diffusion/routine.f",
                reference_instructions="https://www.sns.ias.edu/~jnb/SNdata/Export/Diffusion/README",
                limitations=["Thoul reference uses classical ideal electrons and ions, its published collision logarithms, and all freely moving metals.",
                             "Native ideal-force control retains Ember's electron degeneracy and collision operator. It isolates a thermodynamic force change, not a complete classical limit.",
                             "The hydrostatic reference substitutes an ideal pressure gradient consistent with actual gravity; it is a separate forcing diagnostic.",
                             "Solar normalization R_sun=6.957e10 cm, tau0=6e13 yr, T0=1e7 K and rho0=100 g/cm3 are explicit.",
                             "Rates on faces within a mixed convective region are diagnostics, not rates actually used to separate that region.",
                             "A classical comparison cannot independently validate strongly degenerate core transport."])
    write_result(args.output,result)
    print(json.dumps(dict(outcome=result["outcome"],maximum_precision_difference=max(precision_error),failures=failures,
                          summary=[dict(version=r["version"],q=r["q"],eta=r["native"]["eta"],
                                        native=r["native"]["native_H_velocity"],reference=r["reference"]["ideal_pressure_gradient"]["velocity_cm_s"][0],
                                        ratio=r["native_over_classical_reference_H"],ideal_force_ratio=r["ideal_force_over_classical_reference_H"])
                                   for r in records]),indent=2))


if __name__=="__main__":
    main()
