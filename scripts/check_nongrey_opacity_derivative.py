#!/usr/bin/env python3
"""Probe the native Rybicki opacity derivative and state restoration."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess


def probe(source, work):
    source = Path(source).resolve()
    work = Path(work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    original = source / "tlusty208.f"
    match = re.search(r"^      SUBROUTINE OPACTR\b.*?^      END\s*$",
                      original.read_text(), re.M | re.S)
    if match is None:
        raise ValueError("missing source OPACTR routine")
    kernel = work / "opacity.f"
    kernel.write_text(match.group()+"\n")
    parameter=re.search(r'PARAMETER\s*\(DELT=([0-9.DEd e+\-]+)\)',match.group())
    if parameter is None:raise ValueError('missing opacity derivative step')
    step=float(parameter[1].replace('D','E').replace('d','e'))
    if not 0<step<.1:raise ValueError('invalid opacity derivative step')
    template = Path(__file__).with_name("nongrey_opacity_derivative_probe.f").resolve()
    driver=work/'driver.f'
    text=template.read_text()
    anchor='            EXPECT=AK*(1.01D0**1.4D0-1.D0)/(.01D0*TBASE(ID))'
    if text.count(anchor)!=1:raise ValueError('unrecognized analytic opacity probe')
    number=format(step,'.12e').replace('e','D')
    text=text.replace(anchor,f'            STEP={number}\n'
        '            EXPECT=AK*((1.D0+STEP)**1.4D0-1.D0)/\n'
        '     *             (STEP*TBASE(ID))')
    driver.write_text(text)
    executable = work / "probe"
    subprocess.run(["gfortran", "-O2", "-fno-automatic", "-std=legacy",
                    "-fallow-argument-mismatch", "-fcheck=bounds", "-I"+str(source),
                    str(driver), str(kernel), "-o", str(executable)], check=True)
    output = subprocess.check_output([str(executable)], text=True)
    (work / "output.txt").write_text(output)
    rows = [list(map(float, line.split())) for line in output.splitlines()]
    if len(rows) != 6 or any(len(row) != 7 or not all(math.isfinite(v) for v in row) for row in rows):
        raise ValueError("incomplete opacity derivative probe")
    return {"columns": ["frequency_index", "depth_index", "derivative_relative_error",
                        "temperature_relative_change", "density_relative_change",
                        "molecular_weight_relative_change", "gas_pressure_relative_change"],
            "rows": rows, "maximum_absolute_error": max(abs(v) for row in rows for v in row[2:]),
            "opacity_derivative_relative_step":step,
            "analytic_secant_truncation_relative_error":math.expm1(1.4*math.log1p(step))/(1.4*step)-1,
            "source_sha256": hashlib.sha256(original.read_bytes()).hexdigest(),
            "kernel_sha256": hashlib.sha256(kernel.read_bytes()).hexdigest(),
            "driver_sha256": hashlib.sha256(driver.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("work", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    result = probe(args.source, args.work)
    text = json.dumps(result, indent=2)+"\n"
    if args.report:
        args.report.write_text(text)
    print(text, end="")
    if result["maximum_absolute_error"] > 1e-10:
        raise SystemExit("Opacity derivative/state restoration check failed")


if __name__ == "__main__":
    main()
