#!/usr/bin/env python3
"""Compare the complete 0.1 Msun atmosphere family with direct models and EOS."""
import argparse
from bisect import bisect_right
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess
import tempfile

from import_nongrey_grid import import_grid, read_text, source_state


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def interpolate(spec, states, point):
    choices=[]
    for name,value in zip(["hydrogen","helium3","teff_K","log_g"],point):
        axis=spec[name]
        if not axis[0] <= value <= axis[-1]:
            raise ValueError("held-out state outside the atmosphere grid")
        j=min(bisect_right(axis,value)-1,len(axis)-2)
        transform=math.log if name=="teff_K" else lambda x:x
        fraction=(transform(value)-transform(axis[j]))/(transform(axis[j+1])-transform(axis[j]))
        choices.append([(axis[j],1-fraction),(axis[j+1],fraction)])
    result={k:0. for k in ["T","Pgas","source_density"]}
    corners=[]
    for cell in itertools.product(*choices):
        key=tuple(v[0] for v in cell)
        weight=math.prod(v[1] for v in cell)
        if not weight:
            continue
        corners.append({"coordinates":key,"weight":weight})
        for name in result:
            result[name]+=weight*math.log(states[key][name])
    return {k:math.exp(v) for k,v in result.items()},corners


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest",type=Path)
    parser.add_argument("eos_probe",type=Path)
    parser.add_argument("--gs98",action="store_true",help="compare with the matched metal-bearing EOS")
    parser.add_argument('--eos-only',action='store_true',help='omit heldout gas comparisons when auditing different atmosphere physics')
    parser.add_argument("--output",type=Path,default=Path("docs/results/nongrey_family_audit.json"))
    args=parser.parse_args()
    spec=json.loads(args.manifest.read_text())
    from import_nongrey_grid import CONDENSATE_CALCULATION
    if spec['calculation']==CONDENSATE_CALCULATION and not args.eos_only:
        raise ValueError('gas heldouts do not match condensate physics; use --eos-only and separate condensate heldouts')
    metals=[.0040277310873548045,4.525540547589668e-5,.001113471054023328,
            .00921354245314597,.0056]
    if len(spec["metals"])!=5 or any(abs(a-b)>1e-14 for a,b in zip(spec["metals"],metals)):
        raise ValueError("EOS probe requires the fixed 0.1 Msun Z=.02 metal inventory")
    with tempfile.TemporaryDirectory() as temporary:
        states=import_grid(args.manifest,Path(temporary)/"grid.dat")
    inputs="".join(" ".join(str(v) for v in [*key[:2],s["T"],s["Pgas"],s["source_density"]])+"\n"
                   for key,s in states.items())
    output=subprocess.check_output([str(args.eos_probe.resolve())]+(['--gs98'] if args.gs98 else []),input=inputs,text=True)
    eos=[]
    for (key,state),line in zip(states.items(),output.splitlines(),strict=True):
        values=list(map(float,line.split()))
        if len(values)!=7 or not all(math.isfinite(v) for v in values) or values[5]<=0:
            raise ValueError("invalid EOS probe response")
        expected=[*key[:2],state["T"],state["Pgas"],state["source_density"]]
        if values[:5]!=expected:
            raise ValueError("EOS probe returned a different state")
        eos.append({"coordinates":key,"T":state["T"],"Pgas":state["Pgas"],
                    "source_density":state["source_density"],"ember_density":values[5],
                    "relative_density_difference":values[6]})
    checks={} if args.eos_only else json.loads(Path("docs/results/nongrey_atmosphere_validation.json").read_text())
    archive=Path("data/atmosphere/sources/nongrey_validation")
    comparisons={}
    cases=[] if args.eos_only else [("composition","heldout-composition-final"),("four_axis","heldout-four-axis-final")]
    for label,name in cases:
        record=next(r for r in checks["records"] if r["name"]==name)
        for filename,digest in record["files"].items():
            if sha(archive/filename)!=digest:
                raise ValueError("held-out source checksum mismatch")
        point=tuple(record[k] for k in ["XH","X3","Teff_K","log_g"])
        direct=source_state(read_text(archive/f"{name}-run.log.gz"),
                            read_text(archive/f"{name}-fort.9.gz"),*point[2:],record["opacity"],spec["tau"])
        estimated,corners=interpolate(spec,states,point)
        comparisons[label]={"coordinates":point,"direct_record":name,
            "direct":{k:direct[k] for k in estimated},"interpolated":estimated,
            "relative_difference":{k:estimated[k]/direct[k]-1 for k in estimated},"corners":corners}
    report={"description":"Complete-family matching-state audit; local comparisons, not a global physical error bound",
            "metal_inventory":"gs98" if args.gs98 else "carried_isotopes; metal-as-He EOS proxy",
            "manifest_sha256":sha(args.manifest),"models":len(states),
            "eos_probe_sha256":sha(args.eos_probe),"eos_probe_source_sha256":sha(Path(__file__).with_name("nongrey_eos_probe.cpp")),
            "interpolation":comparisons,"density_comparison":eos,
            "density_relative_difference_range":[min(r["relative_density_difference"] for r in eos),max(r["relative_density_difference"] for r in eos)],
            "note":"Source density interpolation is diagnostic only. Runtime density comes from Ember's EOS at the actual composition."}
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"models":len(states),"density_range":report["density_relative_difference_range"],
                      "interpolation":{k:v["relative_difference"] for k,v in comparisons.items()}},indent=2))


if __name__=="__main__":
    main()
