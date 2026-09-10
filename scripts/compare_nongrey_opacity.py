#!/usr/bin/env python3
"""Compare source absorption tables at identical material/frequency nodes.

Requires NumPy only for this scientific audit, not generation/import/build.
Reports monochromatic differences and Planck/Rosseland absorption means.
These sampled means are diagnostics, not atmosphere flux or accuracy tests.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from nongrey_opacity import read_table
from prepare_nongrey_sources import digest


def compare(candidate, reference, interpolate_reference=False):
    a, b = read_table(candidate), read_table(reference)
    keys = ["abundance","ifmol","tmolim","flags","frequency"]
    if not interpolate_reference:
        keys += ["shape","log_temperature","log_density"]
    for key in keys:
        if a[key] != b[key]:
            raise ValueError(f"source tables differ in {key}")
    frequency = np.array(a["frequency"])
    weights = np.abs(np.gradient(frequency)); weights[[0,-1]] *= .5
    ka = np.exp(np.array(a["log_opacity"],dtype=float).reshape(a["shape"]))
    logb = np.array(b["log_opacity"],dtype=float).reshape(b["shape"])
    if interpolate_reference:
        # Reproduce OPCTAB's bilinear interpolation of ln(kappa) in
        # ln(T), ln(rho), with strict support instead of trial clipping.
        sampled = np.empty(a["shape"])
        def bracket(axis, value):
            if len(axis) < 2 or not axis[0] <= value <= axis[-1]:
                raise ValueError("independent state outside interpolation support")
            j = min(max(int(np.searchsorted(axis,value))-1,0),len(axis)-2)
            return j,(value-axis[j])/(axis[j+1]-axis[j])
        for it, t in enumerate(a["log_temperature"]):
            jt,wt = bracket(b["log_temperature"],t)
            for ir, r in enumerate(a["log_density"]):
                jr,wr = bracket(b["log_density"],r)
                sampled[:,ir,it] = ((1-wt)*((1-wr)*logb[:,jr,jt]+wr*logb[:,jr+1,jt])
                    + wt*((1-wr)*logb[:,jr,jt+1]+wr*logb[:,jr+1,jt+1]))
        logb = sampled
    kb = np.exp(logb)
    rows = []
    for it,t in enumerate(np.exp(a["log_temperature"])):
        # Retain the source constants for a source-to-source comparison.
        z = 6.6256e-27*frequency/(1.38054e-16*t)
        planck = frequency**3/np.expm1(z)
        rosseland = planck*z/(-np.expm1(-z))
        for ir,rho in enumerate(np.exp(a["log_density"])):
            ac,bc = ka[:,ir,it],kb[:,ir,it]
            def means(k):
                return [float(np.sum(weights*planck*k)/np.sum(weights*planck)),
                        float(np.sum(weights*rosseland)/np.sum(weights*rosseland/k))]
            ma,mb = means(ac),means(bc)
            rows.append({"T_K":float(t),"rho_baryon_g_cm3":float(rho),
                "candidate_means_cm2_g":ma,"reference_means_cm2_g":mb,
                "means_relative_change":[x/y-1 for x,y in zip(ma,mb)],
                "monochromatic_relative_change_quantiles":np.quantile(ac/bc-1,[0,.01,.5,.99,1]).tolist()})
    return {"candidate_sha256":digest(candidate),"reference_sha256":digest(reference),
            "reference_material_interpolation":interpolate_reference,
            "means_order":["Planck absorption","Rosseland absorption"],
            "quantiles":[0,.01,.5,.99,1],"rows":rows}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("candidate",type=Path);p.add_argument("reference",type=Path);p.add_argument("output",type=Path)
    p.add_argument("--interpolate-reference",action="store_true",
                   help="compare independent material states with the reference table interpolant")
    a=p.parse_args();a.output.write_text(json.dumps(compare(a.candidate,a.reference,a.interpolate_reference),indent=2)+"\n")


if __name__ == "__main__":main()
