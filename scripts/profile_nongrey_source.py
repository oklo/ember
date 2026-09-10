#!/usr/bin/env python3
"""Time complete TLUSTY phases in an isolated, scientifically checked replay.

The source is copied and instrumented with CPU_TIME; the production executable
and grid cells are untouched. Start from a completed model's final structure.
"""
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess

from generate_nongrey_grid import atmosphere_inputs, composition, digest, execute
from import_nongrey_grid import source_inputs, source_state


def timed_call(text, statement, phase):
    if text.count(statement) != 1:
        raise ValueError(f"expected one {statement!r}")
    return text.replace(statement, "      CALL CPU_TIME(EMBER_T0)\n" + statement +
                        "\n      CALL CPU_TIME(EMBER_T1)\n" +
                        f"      WRITE(6,*) 'EMBER CPU {phase}: ',EMBER_T1-EMBER_T0")


def instrument(source):
    split = source.index("      BLOCK DATA")
    main, rest = source[:split], source[split:]
    main = main.replace("      OPEN(UNIT=91", "      DOUBLE PRECISION EMBER_T0,EMBER_T1,EMBER_CHEM\n" +
                        "      INTEGER EMBER_CALLS\n" +
                        "      COMMON/EMBERCPU/EMBER_CHEM,EMBER_CALLS\n" +
                        "      EMBER_CHEM=0.D0\n      EMBER_CALLS=0\n" +
                        "      OPEN(UNIT=91", 1)
    main = main.replace("   20 CONTINUE", "   20 CONTINUE\n" +
                        "      WRITE(6,*) 'EMBER KERNEL chemical_equilibrium: ',\n" +
                        "     *            EMBER_CHEM,EMBER_CALLS")
    for statement, phase in [("      CALL START", "initialization"),
                             ("      CALL RESOLV", "formal_solution")]:
        main = timed_call(main, statement, phase)
    start = rest.index("      SUBROUTINE RYBSOL\n")
    end = rest.index("      SUBROUTINE RYBMAT(IJ)", start)
    kernel = rest[start:end]
    kernel = kernel.replace("      dimension pop1", "      DOUBLE PRECISION EMBER_T0,EMBER_T1\n" +
                            "      dimension pop1", 1)
    marker = "      DO IJ=1,NFREQ"
    kernel = kernel.replace(marker, "      CALL CPU_TIME(EMBER_T0)\n" + marker, 1)
    marker = "      DO ID=1,ND\n         ABROSD(ID)="
    if kernel.count(marker) != 1:
        raise ValueError("radiative-loop endpoint changed")
    kernel = kernel.replace(marker, "      CALL CPU_TIME(EMBER_T1)\n" +
                            "      WRITE(6,*) 'EMBER CPU radiative_matrix: ',\n" +
                            "     *            EMBER_T1-EMBER_T0\n" + marker)
    for statement, phase in [("      CALL RYBENE", "energy_convection"),
                             ("      CALL LINEQS(WM,WR,CHANGT,ND,MDEPTH)", "dense_solve"),
                             ("      CALL RYBCHN(CHANGT)", "structure_update")]:
        kernel = timed_call(kernel, statement, phase)
    result = main + rest[:start] + kernel + rest[end:]
    declaration = "      SUBROUTINE RUSSEL(TEM,PG)"
    if result.count(declaration) != 1:
        raise ValueError("chemical kernel signature changed")
    result = result.replace(declaration, "      SUBROUTINE EMBER_RUSSEL(TEM,PG)")
    return result + """
      SUBROUTINE RUSSEL(TEM,PG)
      INCLUDE 'IMPLIC.FOR'
      DOUBLE PRECISION EMBER_T0,EMBER_T1,EMBER_CHEM
      INTEGER EMBER_CALLS
      COMMON/EMBERCPU/EMBER_CHEM,EMBER_CALLS
      CALL CPU_TIME(EMBER_T0)
      CALL EMBER_RUSSEL(TEM,PG)
      CALL CPU_TIME(EMBER_T1)
      EMBER_CHEM=EMBER_CHEM+EMBER_T1-EMBER_T0
      EMBER_CALLS=EMBER_CALLS+1
      RETURN
      END
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prepared", type=Path)
    parser.add_argument("specification", type=Path)
    parser.add_argument("model", type=Path)
    parser.add_argument("work", type=Path)
    args = parser.parse_args()
    prepared = json.loads(args.prepared.read_text())
    spec = json.loads(args.specification.read_text())
    record = json.loads((args.model / "validated.json").read_text())
    args.work.mkdir(parents=True, exist_ok=False)
    source_dir = Path(prepared["tlusty"]).parent
    build = args.work / "source"
    build.mkdir()
    for path in source_dir.glob("*.FOR"):
        shutil.copy2(path, build / path.name)
    original = source_dir / "tlusty208.f"
    source = build / "tlusty208.f"
    source.write_text(instrument(original.read_text()))
    executable = build / "tlusty.exe"
    flags = ["-O2", "-g", "-fno-automatic", "-std=legacy", "-fallow-argument-mismatch",
             "-fcheck=bounds", "-fbacktrace"]
    with (args.work / "build.log").open("w") as log:
        subprocess.run(["gfortran", *flags, "tlusty208.f", "-o", str(executable.resolve())],
                       cwd=build, stdout=log, stderr=subprocess.STDOUT, check=True)
    x, y, teff, logg = [record[k] for k in ["XH", "X3", "teff_K", "log_g"]]
    abundance, masses = composition(x, y, spec["metals"])
    model = args.work / "model"
    table = args.model.parent / "opacity" / "fort.63"
    atmosphere_inputs(model, prepared, spec, table, abundance, masses, teff, logg,
                      (args.model / "fort.7").read_text())
    execute(str(executable.resolve()), model, ["fort.7", "fort.9"])
    log = (model / "run.log").read_text()
    state = source_state(log, (model / "fort.9").read_text(), teff, logg,
                         {"temperature_K": [10**spec["log_temperature"][1], 10**spec["log_temperature"][2]],
                          "density_g_cm3": [10**spec["log_density"][1], 10**spec["log_density"][2]]},
                         spec["tau"])
    source_inputs({k: (model / f).read_text() for k, f in
                   [("atmosphere_input", "fort.5"), ("element_masses", "ember-masses.dat"),
                    ("parameters", "tas"), ("initial_structure", "fort.8")]},
                  spec, x, y, teff, logg, log)
    phases = {}
    for phase, value in re.findall(r"EMBER CPU (\w+):\s*([\d.EeDd+-]+)", log):
        phases.setdefault(phase, []).append(float(value.replace("D", "E")))
    chemical = re.findall(r"EMBER KERNEL chemical_equilibrium:\s*([\d.EeDd+-]+)\s+(\d+)", log)
    if len(chemical) != 1:
        raise ValueError("missing final chemical kernel timing")
    report = {"description": "Isolated replay from an accepted final atmosphere; process CPU seconds, not concurrent wall time",
              "source_sha256": digest(original), "instrumented_source_sha256": digest(source),
              "executable_sha256": digest(executable), "driver_sha256": digest(Path(__file__)),
              "compiler_flags": flags, "model": {k: record[k] for k in ["XH", "X3", "teff_K", "log_g"]},
              "depths": spec["depths"], "frequencies": spec["atmosphere_frequencies"],
              "cpu_seconds_by_call": phases,
              "cpu_seconds_by_phase": {k: sum(v) for k, v in phases.items()},
              "chemical_kernel_cpu_seconds": float(chemical[0][0].replace("D", "E")),
              "chemical_kernel_calls": int(chemical[0][1]),
              "timing_note": "Chemical kernel time is included within the phase times; do not add it again. Final formal-solution time includes validation diagnostics.",
              "wall_seconds": json.loads((model / "completed.json").read_text())["seconds"],
              "diagnostics": state,
              "relative_matching_state_change": {k: state[k]/record["diagnostics"][k]-1 for k in ["T", "Pgas"]}}
    (args.work / "report.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
