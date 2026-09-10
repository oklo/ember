#!/usr/bin/env python3
"""Check the native RYBENE convection matrix against finite differences.

Uses the actual CONVEC and RYBENE routines with an analytic ideal molecular
EOS to isolate matrix assembly. Pressure, opacity and the divergence density
prefactor stay fixed. Requires gfortran and a prepared TLUSTY source directory.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess


def probe(source, work):
    source=Path(source).resolve(); work=Path(work).resolve()
    work.mkdir(parents=True,exist_ok=True)
    text=(source/"tlusty208.f").read_text()
    parts=[]
    for name in ["CONVEC","RYBENE"]:
        match=re.search(r"^      SUBROUTINE "+name+r"\b.*?^      END\s*$",text,re.M|re.S)
        if match is None: raise ValueError("missing source convection routine")
        parts.append(match.group())
    kernel=work/"convection.f"; kernel.write_text("\n\n".join(parts)+"\n")
    driver=Path(__file__).with_name("nongrey_convection_probe.f").resolve()
    executable=work/"probe"
    subprocess.run(["gfortran","-O2","-fno-automatic","-std=legacy",
        "-fallow-argument-mismatch","-fcheck=bounds","-I"+str(source),
        str(driver),str(kernel),"-o",str(executable)],check=True)
    output=subprocess.check_output([str(executable)],text=True)
    (work/"output.txt").write_text(output)
    rows=[list(map(float,line.split())) for line in output.splitlines()]
    if len(rows)!=9 or any(len(row)!=5 for row in rows):
        raise ValueError("incomplete source derivative probe")
    return {"columns":["pressure_case","temperature_index","source_matrix",
                       "finite_difference","relative_error"],
            "rows":rows,"maximum_relative_error":max(abs(row[-1]) for row in rows),
            "source_sha256":hashlib.sha256((source/"tlusty208.f").read_bytes()).hexdigest(),
            "kernel_sha256":hashlib.sha256(kernel.read_bytes()).hexdigest(),
            "driver_sha256":hashlib.sha256(driver.read_bytes()).hexdigest()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",type=Path)
    parser.add_argument("work",type=Path)
    parser.add_argument("--report",type=Path)
    args=parser.parse_args(); result=probe(args.source,args.work)
    text=json.dumps(result,indent=2)+"\n"
    if args.report: args.report.write_text(text)
    print(text,end="")
    if result["maximum_relative_error"]>1e-4:
        raise SystemExit("Convection matrix fails the independent derivative check")


if __name__=="__main__":main()
