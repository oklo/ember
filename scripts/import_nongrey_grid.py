#!/usr/bin/env python3
"""Validate archived TLUSTY structures and import a complete non-grey family.

Requires a manifest produced by generate_nongrey_grid.py and the original
source outputs. Never fills missing cells or extrapolates to matching depth.
"""
import argparse
from bisect import bisect_right
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import re

GAS_CALCULATION = 'TLUSTY208/SYNSPEC54 non-grey LTE convection'
CONDENSATE_CALCULATION = GAS_CALCULATION+' with FastChem4 equilibrium depletion'


def source_failure(log):
    """Recognize fatal source messages even when Fortran STOP exits zero."""
    match=re.search(r'runtime error|DOES NOT CONVERGE|NOT CONVERGE IN|\bNaN\b|\bInfinity\b|SIGBUS|SIGSEGV|'
                    r'^[ \t]*STOP[ \t]+[^\r\n]+|\*{4}[ \t]+STOP in[^\r\n]+|h2minus:Stop',log,re.I|re.M)
    return match.group(0).strip() if match else None


def read_text(path):
    path = Path(path)
    return gzip.decompress(path.read_bytes()).decode() if path.suffix == ".gz" else path.read_text()


def source_inputs(inputs, spec, x, y, teff, logg, log):
    """Check archived input physics independently of the model's grid label."""
    from generate_nongrey_grid import composition
    abundance, masses = composition(x, y, spec["metals"])
    rows = inputs["atmosphere_input"].splitlines()
    if len(rows) != 98 or [float(v) for v in rows[0].split()] != [teff, logg]:
        raise ValueError("source input Teff/gravity mismatch")
    expected_start = ["T","F"] if "initial_structure" in inputs else ["T","T"]
    if rows[1].split() != expected_start or rows[3].strip() != "0" or rows[4].strip() != "92":
        raise ValueError("source must use LTE full-opacity atmosphere mode")
    if "initial_structure" in inputs:
        from generate_nongrey_grid import resample_initial_structure
        if int(inputs["initial_structure"].split()[0]) != spec["depths"]:
            raise ValueError("initial structure depth count mismatch")
        resample_initial_structure(inputs["initial_structure"],spec["depths"])
    for i,(row, expected) in enumerate(zip(rows[5:97], abundance)):
        mode, value, modifier = row.split()
        if mode != ("1" if i < 30 else "0") or modifier != "0" or not math.isclose(float(value), expected, rel_tol=1e-13):
            raise ValueError("source atmosphere abundance mismatch")
    actual = [float(v) for v in inputs["element_masses"].split()]
    if len(actual) != len(masses) or any(not math.isclose(a,b,rel_tol=1e-13) for a,b in zip(actual,masses)):
        raise ValueError("source baryonic mass mapping mismatch")
    marker = re.findall(r"EMBER ELEMENT MASSES:\s*([0-9.eEdD+-]+)\s+([0-9.eEdD+-]+)",log)
    if not marker or any(not math.isclose(float(a.replace("D","E")),b,rel_tol=1e-12) for a,b in zip(marker[-1],masses)):
        raise ValueError("source did not load the baryonic element masses")
    chemistry=re.findall(r"EMBER MOLECULAR EQUILIBRIUM TOLERANCE:\s*([0-9.eEdD+-]+)",log)
    if not chemistry or float(chemistry[-1].replace("D","E")) != 1e-8:
        raise ValueError("source molecular equilibrium tolerance mismatch")
    settings = {k:float(v.replace("D","E")) for k,v in re.findall(
        r"([A-Z][A-Z0-9]*)\s*=\s*([0-9.eEdD+-]+)",inputs["parameters"].upper())}
    required = {"IOPTAB":-1,"IFRYB":1,"IFMOL":1,"TMOLIM":10000,"IDLST":0,"IFRAYL":1,
                "HMIX0":spec["alpha"],"IFRSET":spec["atmosphere_frequencies"],"ND":spec["depths"],
                "TAUFIR":spec["tau_top"],"TAULAS":spec["tau_bottom"],"TAUDIV":.01,
                "ILGDER":1,"NITER":200,"DPSILT":spec.get("temperature_step_limit",1.03),
                "DERT":spec.get("convection_derivative_step",.001)}
    if any(settings.get(k) != v for k,v in required.items()) or not 0 < settings.get("CHMAX",1) <= 1e-6:
        raise ValueError("source atmosphere physical/numerical settings mismatch")


def source_state(log, convergence, teff, logg, opacity, tau=100., max_flux_error=.002,
                 max_correction=1e-5, require_chemical_closure=True):
    """Return matching state and diagnostics, checking independent outputs."""
    failure=source_failure(log)
    if failure:
        raise ValueError(f"failed or nonfinite source calculation: {failure}")
    if "FINAL MODEL ATMOSPHERE" not in log:
        raise ValueError("source has no final atmosphere")
    rows = []
    for line in log.rsplit("FINAL MODEL ATMOSPHERE", 1)[1].splitlines():
        tokens = line.replace("D", "E").split()
        if len(tokens) != 11 or not tokens[0].isdigit():
            continue
        row = [float(v) for v in tokens]
        if int(row[0]) != len(rows) + 1 or not all(math.isfinite(v) for v in row):
            raise ValueError("invalid source profile ordering or values")
        if any(row[j] <= 0 for j in [1, 2, 3, 4, 5, 6]):
            raise ValueError("nonpositive source thermodynamic state")
        if rows and (row[1] <= rows[-1][1] or row[2] <= rows[-1][2]):
            raise ValueError("unordered source column or optical depth")
        rows.append(row)
    if len(rows) < 20:
        raise ValueError("missing source depth points")
    chemical_error = None
    if require_chemical_closure:
        chemical = [line.split(":",1)[1].replace("D","E").split()
                    for line in log.splitlines() if line.strip().startswith("EMBER CHEMICAL DENSITY:")]
        if len(chemical) != len(rows):
            raise ValueError("missing independent chemical density check")
        errors = []
        for expected, words in zip(rows,chemical):
            if len(words) != 5:
                raise ValueError("invalid chemical density record")
            i,t,p,rho,actual = map(float,words)
            if i != expected[0] or any(not math.isfinite(v) or v <= 0 for v in [t,p,rho,actual]):
                raise ValueError("invalid chemical density state")
            if any(not math.isclose(a,b,rel_tol=1e-8) for a,b in zip([t,p,rho],[expected[3],expected[6],expected[5]])):
                raise ValueError("chemical density check does not match final structure")
            errors.append(abs(rho/actual-1))
        chemical_error = max(errors)
        if chemical_error > .002:
            raise ValueError(f"molecular density closure error: {chemical_error:g}")
    changes = []
    for line in convergence.splitlines():
        words = line.replace("D", "E").split()
        if len(words) == 9 and words[0].isdigit() and words[1].isdigit():
            values = [float(v) for v in words]
            if not all(math.isfinite(v) for v in values):
                raise ValueError("nonfinite convergence diagnostic")
            changes.append(values)
    if not changes:
        raise ValueError("missing independent convergence output")
    last = int(changes[-1][0]); final = [r for r in changes if r[0] == last]
    if len(final) != len(rows) or {int(r[1]) for r in final} != set(range(1, len(rows) + 1)):
        raise ValueError("incomplete final correction vector")
    correction = max(abs(r[6]) for r in final)
    flux_error = max(abs(r[10] - 1) for r in rows)
    if correction > max_correction or flux_error > max_flux_error:
        raise ValueError(f"unconverged source: correction={correction:g}, flux error={flux_error:g}")
    # This table column and reported surface H are independent of the
    # requested Teff. Source uses H=F/(4*pi), not flux F.
    surface = re.findall(r"TOTAL SURFACE FLUX\s*([0-9.+EDed-]+)", log)
    if not surface:
        raise ValueError("missing emergent bolometric flux")
    flux = 4 * math.pi * float(surface[-1].replace("D", "E"))
    effective_error = abs(flux / (5.670374419e-5 * teff**4) - 1)
    if effective_error > max_flux_error:
        raise ValueError("emergent flux does not match requested Teff")
    # Bound the converged physical profile by the original T/rho table.
    # Generated tables have common logarithmic density axes at every T.
    for r in rows:
        if not (opacity["temperature_K"][0] <= r[3] <= opacity["temperature_K"][-1]
                and opacity["density_g_cm3"][0] <= r[5] <= opacity["density_g_cm3"][-1]):
            raise ValueError("converged atmosphere uses opacity outside source support")
    # Radiation acceleration is very small on this branch. Check the
    # hydrostatic gas-pressure column independently; this also catches a
    # wrong requested gravity or pressure convention.
    hydro_error = max(abs(r[6] / (r[1] * 10**logg) - 1) for r in rows)
    if hydro_error > .01:
        raise ValueError("source gas pressure and requested gravity disagree")
    depths = [r[2] for r in rows]
    if not depths[0] < tau < depths[-1]:
        raise ValueError("matching optical depth is not bracketed by source")
    j = bisect_right(depths, tau) - 1
    a, b = rows[j:j+2]
    w = math.log(tau/a[2]) / math.log(b[2]/a[2])
    t = math.exp((1-w)*math.log(a[3]) + w*math.log(b[3]))
    pg = math.exp((1-w)*math.log(a[6]) + w*math.log(b[6]))
    density = math.exp((1-w)*math.log(a[5]) + w*math.log(b[5]))
    column = math.exp((1-w)*math.log(a[1]) + w*math.log(b[1]))
    return {"T": t, "Pgas": pg, "source_density": density, "column_mass": column,
            "flux_error": flux_error, "correction": correction,
            "effective_flux_error": effective_error, "hydrostatic_error": hydro_error,
            "chemical_density_error": chemical_error,
            "iterations": last, "depths": len(rows), "tau_bracket": [a[2], b[2]],
            "temperature_range": [min(r[3] for r in rows), max(r[3] for r in rows)],
            "density_range": [min(r[5] for r in rows), max(r[5] for r in rows)]}


def import_grid(manifest, output):
    manifest = Path(manifest); spec = json.loads(manifest.read_text())
    if spec.get("format") != 1 or spec.get("calculation") not in [GAS_CALCULATION,CONDENSATE_CALCULATION]:
        raise ValueError("unrecognized atmosphere source manifest")
    if not spec.get("source") or not spec.get("approximation") or not spec.get("provenance"):
        raise ValueError("source provenance/approximations required")
    axes = [spec[k] for k in ["hydrogen", "helium3", "teff_K", "log_g"]]
    for a in axes:
        if len(a) < 2 or any(not math.isfinite(v) for v in a) or any(v >= w for v,w in zip(a,a[1:])):
            raise ValueError("invalid source axis")
    if len(spec["metals"]) != 5 or any(v < 0 or not math.isfinite(v) for v in spec["metals"]):
        raise ValueError("invalid metal inventory")
    if axes[0][0] < 0 or axes[1][0] < 0 or axes[0][-1]+axes[1][-1]+sum(spec["metals"]) > 1+1e-12:
        raise ValueError("unphysical composition grid")
    expected = set(itertools.product(*axes)); records = {}
    depleted=spec['calculation']==CONDENSATE_CALCULATION
    if depleted:
        from validate_condensate_model import validate
        if spec.get('condensates')!={'mode':'equilibrium','grain_opacity':0,
                'grain_enthalpy':'reject condensing layers with convective flux fraction above 1e-8'}:
            raise ValueError('unrecognized condensate physical approximation')
        if not spec.get('archive_files_sha256') or not spec.get('opacity_planes'):
            raise ValueError('condensate source archive and opacity receipts required')
        for name,expected_hash in spec['archive_files_sha256'].items():
            path=manifest.parent/name
            if not path.resolve().is_relative_to(manifest.parent.resolve()):
                raise ValueError('condensate source path leaves archive')
            if hashlib.sha256(path.read_bytes()).hexdigest()!=expected_hash:
                raise ValueError('condensate archive checksum mismatch')
    for record in spec["models"]:
        key = tuple(record[k] for k in ["XH", "X3", "teff_K", "log_g"])
        if key not in expected or key in records:
            raise ValueError("duplicate or unexpected source model")
        if depleted:
            directory=manifest.parent/record['condensate_model']
            if not directory.resolve().is_relative_to(manifest.parent.resolve()):
                raise ValueError('condensate model path leaves archive')
            planes=[p for p in spec['opacity_planes'] if (p['XH'],p['X3'])==key[:2]]
            if len(planes)!=1:raise ValueError('missing or duplicate condensate opacity plane')
            result=validate(directory,spec,spec['provenance'],key,planes[0]['sha256'])
            records[key]=result['diagnostics']
            continue
        contents = []
        for kind in ["log", "convergence"]:
            path = manifest.parent / record[kind]
            if hashlib.sha256(path.read_bytes()).hexdigest() != record[kind+"_sha256"]:
                raise ValueError("source output checksum mismatch")
            contents.append(read_text(path))
        inputs = {}
        kinds = ["atmosphere_input", "element_masses", "parameters"]
        if "initial_structure" in record: kinds.append("initial_structure")
        for kind in kinds:
            path = manifest.parent / record[kind]
            if hashlib.sha256(path.read_bytes()).hexdigest() != record[kind+"_sha256"]:
                raise ValueError("source input checksum mismatch")
            inputs[kind] = read_text(path)
        source_inputs(inputs, spec, *key, contents[0])
        records[key] = source_state(*contents, key[2], key[3], spec["opacity"], tau=spec["tau"])
        if records[key]["depths"] != spec["depths"]:
            raise ValueError("source profile depth count mismatch")
    if set(records) != expected:
        raise ValueError(f"incomplete physical source grid: {len(records)}/{len(expected)} models")
    lines = ["EMBER_COMPOSITION_ATMOSPHERE 1", "source "+json.dumps(spec["source"]),
             "approximation "+json.dumps(spec["approximation"]), "basis baryon_mass",
             f'tau {spec["tau"]:.17g}', "metals "+" ".join(f"{v:.17g}" for v in spec["metals"])]
    for name, a in zip(["hydrogen", "helium3", "log_teff", "log_g"], axes):
        values = [math.log10(v) for v in a] if name == "log_teff" else a
        lines.append(f"{name} {len(a)} "+" ".join(f"{v:.17g}" for v in values))
    lines.append("data")
    for key in itertools.product(*axes):
        r = records[key]
        lines.append(f'{math.log10(r["T"]):.17g} {math.log10(r["Pgas"]):.17g}')
    Path(output).write_text("\n".join(lines)+"\n")
    return records


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("manifest",type=Path); p.add_argument("output",type=Path)
    a=p.parse_args(); records=import_grid(a.manifest,a.output)
    print(f"Imported {len(records)} converged non-grey source models into {a.output}")


if __name__ == "__main__":
    main()
