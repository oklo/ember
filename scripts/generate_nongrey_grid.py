#!/usr/bin/env python3
"""Calculate a helium-rich non-grey atmosphere family with pinned sources.

Run prepare_nongrey_sources.py first. A JSON specification controls the
physical composition, opacity sampling and atmosphere mesh. Expensive source
runs are restartable only when their complete input fingerprints agree.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import time
from import_nongrey_grid import source_state, source_inputs, import_grid, source_failure
from prepare_nongrey_sources import digest, data_digest, SOURCES, ROOT
from nongrey_opacity import validate_table, merge_isotherms

MU = 1.66053906660e-24
SOURCE_HMASS = 1.67333e-24
EXPLICIT = {1, 2, 6, 7, 8, 11, 12, 13, 14, 20, 26}
CALCULATION = "TLUSTY208/SYNSPEC54 non-grey LTE convection"


def sequence(triad):
    n, lo, hi = triad
    if int(n) != n or n < 2 or n > 100000 or not math.isfinite(lo+hi) or lo >= hi:
        raise ValueError("invalid logarithmic source axis")
    return [10**(lo+(hi-lo)*i/(n-1)) for i in range(n)]


def temperatures(spec):
    reference = sequence(spec["log_temperature"])
    values = spec.get("temperature_K",reference)
    # Explicit values preserve existing isotherms bit for bit when extending
    # a uniform logarithmic grid. They must still describe that same grid.
    if len(values) != len(reference) or any(not math.isfinite(v) or v <= 0 or
            not math.isclose(v,r,rel_tol=1e-12) for v,r in zip(values,reference)):
        raise ValueError("explicit temperatures disagree with logarithmic grid")
    return values


def composition(x, y, metals):
    """GS98 metal proxy with explicit H/He baryonic number counts.

    Non-He elements use the nearest integer mass to the source atomic weight;
    this is a declared representative isotope, not an isotope-resolved line
    library. The mean helium mass gives NHe/rho=(X3/3+X4/4)/mu exactly.
    """
    if len(metals) != 5 or any(v < 0 or not math.isfinite(v) for v in metals):
        raise ValueError("invalid fixed metal abundances")
    z = sum(metals); he4 = 1-x-y-z
    if not x > 0 or not y >= 0 or not he4 >= 0 or not y+he4 > 0:
        raise ValueError("unphysical atmosphere composition")
    elements = json.loads((SOURCES/"synple-elements.json").read_text())
    masses = [float(round(v)) for v in elements["mass"]]
    masses[0] = 1.; masses[1] = (y+he4)/(y/3+he4/4)
    fractions = [0.] * 99
    fractions[0] = x; fractions[1] = y+he4
    request = json.loads((ROOT/"data/opacity/sources/tops_gs98_x070_z020.request.json").read_text())
    tokens = request["mixture"].split()[4:]
    raw = {tokens[i+1].lower(): float(tokens[i]) for i in range(0,len(tokens),2)}
    total = sum(raw.values())
    for i, symbol in enumerate(elements["symbol"][2:],2):
        fractions[i] = z * raw.get(symbol.lower(),0.) / total
    # Positive tiny abundances avoid the source convention ABN=0 => solar.
    abundance = [max(v / (m*x),1e-99) for v,m in zip(fractions,masses)]
    weights = [m*MU/SOURCE_HMASS for m in masses]
    if abs(sum(a*w for a,w in zip(abundance,weights))*x*SOURCE_HMASS/MU-1) > 1e-12:
        raise ValueError("baryonic source normalization failed")
    return abundance, weights


def link(path, target):
    path=Path(path); target=Path(target).resolve()
    if path.is_symlink():
        if path.resolve() == target: return
        path.unlink()
    if path.exists(): raise FileExistsError(path)
    path.symlink_to(target, target_is_directory=target.is_dir())


def input_fingerprint(executable_sha256, directory, archived=False):
    """Hash executable identity and exact physical/numerical source inputs."""
    names = ["fort.5","fort.55","fort.2","fort.15","tas","ember-masses.dat","opacity.sha256","physics.json","fort.8",
             "ember-condensates.cfg","condensate-abundances.dat","condensates.sha256"]
    fingerprint=hashlib.sha256()
    fingerprint.update(executable_sha256.encode())
    for name in names:
        p=directory/name
        compressed=False
        if archived and not p.exists() and p.with_name(name+'.gz').exists():
            p=p.with_name(name+'.gz');compressed=True
        # TLUSTY creates an empty unit 2 while starting a full-table run.
        # It is not an input and must not invalidate a completed run's key.
        if p.exists() and (name != "fort.2" or p.stat().st_size):
            content=gzip.decompress(p.read_bytes()) if compressed else p.read_bytes()
            if name!='fort.2' or content:
                fingerprint.update(name.encode()+content)
    return fingerprint.hexdigest()


def execute(executable, directory, required):
    """Never treat a zero exit status alone as Fortran STOP success."""
    key=input_fingerprint(digest(executable),directory); receipt=directory/"completed.json"
    if receipt.exists():
        old=json.loads(receipt.read_text())
        if old["input_sha256"] == key and all((directory/f).exists() and digest(directory/f)==h for f,h in old["outputs"].items()):
            failure=source_failure((directory/'run.log').read_text())
            if failure:raise RuntimeError(f"cached source calculation failed ({failure}); see {directory/'run.log'}")
            return
    # Restart fresh; old generated files cannot masquerade as completion.
    receipt.unlink(missing_ok=True)
    for name in required:
        (directory/name).unlink(missing_ok=True)
    with (directory/"fort.5").open("rb") as source, (directory/"run.log").open("wb") as log:
        start=time.time()
        proc=subprocess.Popen([str(executable)],cwd=directory,stdin=source,stdout=log,stderr=subprocess.STDOUT)
        (directory/"running.json").write_text(json.dumps({"pid":proc.pid,"started_unix":start,"input_sha256":key})+"\n")
        code=proc.wait()
    (directory/"running.json").unlink(missing_ok=True)
    text=(directory/"run.log").read_text()
    failure=source_failure(text)
    if code or failure:
        raise RuntimeError(f"source calculation failed ({failure or code}); see {directory/'run.log'}")
    if not all((directory/f).exists() and (directory/f).stat().st_size > 0 for f in required):
        raise RuntimeError(f"source calculation stopped without outputs; see {directory/'run.log'}")
    receipt.write_text(json.dumps({"input_sha256":key,"seconds":time.time()-start,
        "outputs":{f:digest(directory/f) for f in required+["run.log"]}},indent=2)+"\n")


def opacity_inputs(directory, prepared, spec, x, y, temperature=None):
    directory.mkdir(parents=True,exist_ok=True)
    abundance, masses=composition(x,y,spec["metals"])
    text="9999 9.9\nT F\n 'tas'\n50\n99\n"
    for i,a in enumerate(abundance,1): text+=f'{2 if i in EXPLICIT else 1} {a:.17e} 0\n'
    text+=(SOURCES/"synspec-ap18-ions.dat").read_text()
    (directory/"fort.5").write_text(text)
    (directory/"ember-masses.dat").write_text("\n".join(f"{m:.17e}" for m in masses)+"\n")
    (directory/"tas").write_text("ND=1,IFMOL=1,TMOLIM=10000.\nIOH2H2=1,IOH2HE=1,IOH2H1=1,IOHHE=1\n")
    lo,hi=spec["wavelength_A"]; strength=spec["line_threshold"]
    (directory/"fort.55").write_text(f"-3 0 0\n1 0 0 0\n0 0 0 0 0\n1 0 0 1 0\n2 1 1\n"
        f"{lo:.17g} {-hi:.17g} 200 0 {strength:.17g} {spec['synthesis_spacing_A']:.17g}\n3 20 21 22\n{spec['microturbulence_km_s']:.17g}\n")
    t=temperatures(spec) if temperature is None else [temperature]
    r=sequence(spec["log_density"])
    (directory/"fort.2").write_text(f"{len(t)} {t[0]:.17e} {t[-1]:.17e}\n1\n"
        f"{len(r)} {r[0]:.17e} {r[-1]:.17e}\n{spec['opacity_frequencies']} 1 {lo:.17e} {hi:.17e}\n'opacity.bin' 1\n")
    (directory/"physics.json").write_text(json.dumps({"data":prepared["data_sha256"],
        "lines":prepared["line_list_sha256"]},sort_keys=True)+"\n")
    data=Path(prepared["synple"])/"data"; link(directory/"data",data)
    for unit,path in zip([19,20,21,22],prepared["line_lists"]): link(directory/f"fort.{unit}",path)
    for name in ["CIA_H2H2.dat","CIA_H2H.dat","CIA_H2He.dat","CIA_HHe.dat","irwin_bc.dat","tremblay.dat","tsuji.molec_bc2"]:
        link(directory/name,data/name)
    return abundance,masses


def resample_initial_structure(text, depths):
    """Refine a positive LTE molecular seed, preserving its column endpoints.

    This interpolates an initial guess only. Transfer, hydrostatic and energy
    balance must subsequently converge on the new mesh at full resolution.
    """
    words = text.replace("D", "E").split()
    if len(words) < 2:
        raise ValueError("missing initial atmosphere")
    n, parameters = int(words[0]), int(words[1])
    values = [float(v) for v in words[2:]]
    if n < 20 or parameters != -4 or len(values) != 5*n or depths < n:
        raise ValueError("invalid initial molecular LTE structure")
    if any(not math.isfinite(v) or v <= 0 for v in values):
        raise ValueError("nonpositive or nonfinite initial structure")
    mass = values[:n]
    if any(a >= b for a,b in zip(mass,mass[1:])):
        raise ValueError("unordered initial column mass")
    columns = [mass] + list(zip(*(values[n+4*i:n+4*i+4] for i in range(n))))
    refined = []
    for column in columns:
        result = []
        for i in range(depths):
            f = i*(n-1)/(depths-1); j = min(int(f),n-2); w = f-j
            result.append(math.exp((1-w)*math.log(column[j])+w*math.log(column[j+1])))
        refined.append(result)
    rows = [f"{depths} -4"] + [f"{v:.17e}" for v in refined[0]]
    rows += [" ".join(f"{c[i]:.17e}" for c in refined[1:]) for i in range(depths)]
    return "\n".join(rows)+"\n"


def atmosphere_inputs(directory, prepared, spec, table, abundance, masses, teff, logg, initial=None):
    directory.mkdir(parents=True,exist_ok=True)
    text=f"{teff:.17g} {logg:.17g}\nT T\n 'tas'\n0\n92\n"
    # TLUSTY's atomic partition functions support elements through Zn.
    # The GS98 proxy has no nonzero element beyond that range. Preserve
    # the common trace placeholders in molecular chemistry, without
    # activating unsupported atomic partition-function queries.
    if any(v > 1e-90 for v in abundance[30:]):
        raise ValueError("nonzero element beyond source atomic EOS support")
    text+="".join(f"{1 if i < 30 else 0} {v:.17e} 0\n" for i,v in enumerate(abundance[:92]))
    text+="0 0 0 -1 0 0 ' ' ' '\n"
    (directory/"fort.5").write_text(text)
    if initial is not None:
        (directory/"fort.8").write_text(resample_initial_structure(initial, spec["depths"]))
        (directory/"fort.5").write_text(text.replace("T T\n", "T F\n"))
    else:
        (directory/"fort.8").unlink(missing_ok=True)
    (directory/"ember-masses.dat").write_text("\n".join(f"{m:.17e}" for m in masses)+"\n")
    (directory/"tas").write_text("IOPTAB=-1,IFRYB=1,IFMOL=1,TMOLIM=10000.,IDLST=0,IFRAYL=1\n"
        f"HMIX0={spec['alpha']},ITEK=200,IACC=200,IFRSET={spec['atmosphere_frequencies']}\n"
        # Rybicki convection uses logarithmic gradients. CONOUT and the
        # initial model must use the same discretization (ILGDER=1).
        # Save a full initial-structure checkpoint after each iteration.
        # Only the final independent diagnostics can validate a grid cell.
        f"ND={spec['depths']},NITER=200,CHMAX=1.e-6,ILGDER=1,IPRIND=2\n"
        f"DPSILT={spec.get('temperature_step_limit', 1.03)},DERT={spec.get('convection_derivative_step',.001)}\n"
        f"TAUFIR={spec['tau_top']},TAULAS={spec['tau_bottom']},TAUDIV=0.01\n")
    # The table is a dependency, not merely its path.
    (directory/"fort.15").write_text("'opacity.bin' 1\n")
    link(directory/"opacity.bin",table)
    (directory/"opacity.sha256").write_text(digest(table)+"\n")
    (directory/"physics.json").write_text(json.dumps({"data":prepared["data_sha256"]},sort_keys=True)+"\n")
    link(directory/"data",Path(prepared["synple"])/"data")


def convective_tail(initial, depths):
    """Steepen only the tail of a trial structure to seed a convective zone.

    This is a Newton starting guess, not a physical atmosphere prescription.
    The target chemistry, hydrostatic structure and flux must all solve anew.
    """
    words=initial.split();n=int(words[0])
    resample_initial_structure(initial,n)
    if not 2<=depths<n-2:raise ValueError("invalid convective tail depth count")
    mass=list(map(float,words[2:2+n]))
    rows=[list(map(float,words[2+n+4*i:2+n+4*i+4])) for i in range(n)]
    for i in range(n-depths,n):rows[i][0]=max(rows[i][0],rows[i-1][0]*(mass[i]/mass[i-1])**.5)
    return f"{n} -4\n"+'\n'.join(format(v,'.17g') for v in mass)+'\n'+'\n'.join(
        ' '.join(format(v,'.17g') for v in row) for row in rows)+'\n'


def continuation_structure(initial, source_teff, source_logg, teff, logg):
    """Scale a trial LTE profile while preserving its hydrostatic gas pressure.

    T follows the effective-temperature ratio, all number/mass densities
    scale inversely, and column mass scales inversely with gravity. Thus
    both n*k*T and g*m remain fixed. This is an initial guess only: chemistry,
    optical depths and flux must be solved anew at the target coordinates.
    """
    words=initial.replace('D','E').split();n=int(words[0])
    resample_initial_structure(initial,n)
    if min(source_teff,teff)<=0 or not all(math.isfinite(v) for v in
            [source_teff,source_logg,teff,logg]):
        raise ValueError('invalid continuation coordinates')
    ratio=teff/source_teff;gravity=10**(source_logg-logg)
    mass=[float(v)*gravity for v in words[2:2+n]]
    rows=[list(map(float,words[2+n+4*i:2+n+4*i+4])) for i in range(n)]
    for row in rows:
        row[0]*=ratio
        for j in range(1,4):row[j]/=ratio
    result=f'{n} -4\n'+'\n'.join(format(v,'.17g') for v in mass)+'\n'+'\n'.join(
        ' '.join(format(v,'.17g') for v in row) for row in rows)+'\n'
    resample_initial_structure(result,n)
    return result


def archive(path, directory):
    target=directory/(path.name+".gz")
    target.write_bytes(gzip.compress(path.read_bytes(),compresslevel=9,mtime=0))
    return target


def checkpoint_initial_structure(directory, x, y, teff, logg):
    """Read an attested interrupted structure strictly as an initial guess.

    A checkpoint has no final flux/convergence claim and can never supply
    an imported grid cell. Its original input physics and saved bytes are
    checked before a new full atmosphere solve is allowed to use it.
    """
    directory=Path(directory)
    saved=json.loads((directory/"checkpoint.json").read_text())
    if saved.get("kind") != "unconverged_initial_guess":
        raise ValueError("unrecognized atmosphere checkpoint")
    if any(saved.get(k) != v for k,v in zip(["XH","X3","teff_K","log_g"],[x,y,teff,logg])):
        raise ValueError("initial checkpoint grid label mismatch")
    if any(not (directory/f).is_file() or digest(directory/f) != expected
           for f,expected in saved["files"].items()):
        raise ValueError("initial checkpoint checksum mismatch")
    required={"fort.5","fort.7","fort.9","tas","ember-masses.dat","run.log",
              "fort.15","opacity.sha256","physics.json"}
    if not required.issubset(saved["files"]):
        raise ValueError("incomplete initial checkpoint provenance")
    if input_fingerprint(saved["source_executable_sha256"],directory) != saved["input_sha256"]:
        raise ValueError("initial checkpoint input fingerprint mismatch")
    inputs={k:(directory/f).read_text() for k,f in
            [("atmosphere_input","fort.5"),("element_masses","ember-masses.dat"),("parameters","tas")]}
    if (directory/"fort.8").exists():
        if "fort.8" not in saved["files"]:
            raise ValueError("unattested initial checkpoint input structure")
        inputs["initial_structure"]=(directory/"fort.8").read_text()
    source_inputs(inputs,saved["specification"],saved["XH"],saved["X3"],teff,logg,
                  (directory/"run.log").read_text())
    structure=(directory/"fort.7").read_text()
    resample_initial_structure(structure,saved["specification"]["depths"])
    return structure


def completed_initial_structure(directory, spec, teff, logg, depths):
    """Read a completed starting model from working files or the offline archive."""
    directory=Path(directory)

    def contents(name):
        path=directory/name
        return path.read_bytes() if path.exists() else gzip.decompress(
            path.with_name(path.name+".gz").read_bytes())

    saved=json.loads(contents("completed.json"))
    if not {"fort.7","fort.9","run.log"}.issubset(saved["outputs"]):
        raise ValueError("incomplete initial model receipt")
    for filename,expected in saved["outputs"].items():
        if hashlib.sha256(contents(filename)).hexdigest()!=expected:
            raise ValueError("initial model receipt checksum mismatch")
    opacity={"temperature_K":temperatures(spec),"density_g_cm3":sequence(spec["log_density"])}
    state=source_state(contents("run.log").decode(),contents("fort.9").decode(),
                       teff,logg,opacity,spec["tau"],max_flux_error=.05)
    if state["depths"]>depths:
        return None
    structure=contents("fort.7").decode()
    resample_initial_structure(structure,state["depths"])
    return structure


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("prepared",type=Path); p.add_argument("specification",type=Path); p.add_argument("work",type=Path)
    p.add_argument("--jobs",type=int,default=1)
    p.add_argument("--opacity-only",action="store_true")
    p.add_argument("--reuse-opacity",type=Path,help="reuse matching verified isotherm receipts from this work directory")
    p.add_argument("--initial-models",type=Path,help="use completed structures or attested interrupted checkpoints as initial guesses; every final model is solved again")
    p.add_argument("--continuation-models",type=Path,
                   help="prefer nearby accepted atmospheres at the same composition, preserving trial gas pressure while changing Teff/gravity")
    p.add_argument("--plane",type=int,help="compute only this zero-based composition plane")
    p.add_argument("--convective-initialization",type=int,nargs='+',default=[],
                   help="composition planes to precondition with native CONREF before a separate canonical solve")
    p.add_argument("--convective-iterations",type=int,default=3,
                   help="number of initializer iterations with native CONREF (1..200); canonical replay remains unchanged")
    p.add_argument("--initializer",type=Path,help="receipt for a separately built, bounds-checked convective initializer; final source stays canonical")
    p.add_argument("--convective-tail-depths",type=int,default=0,help="seed a convective gradient in this many bottom trial layers before native refinement")
    a=p.parse_args(); prepared=json.loads(a.prepared.read_text()); spec=json.loads(a.specification.read_text())
    initializer=None
    if a.initializer:
        initializer=json.loads(a.initializer.read_text())
        if digest(initializer['executable'])!=initializer['executable_sha256']:
            raise ValueError("initializer executable changed")
        canonical_source=Path(prepared['tlusty_source'])/'tlusty/tlusty208.f'
        if digest(canonical_source)!=initializer['base_source_sha256']:
            raise ValueError("initializer does not derive from the selected canonical source")
    if a.convective_tail_depths and not 2<=a.convective_tail_depths<spec['depths']-2:
        raise ValueError("invalid convective tail size")
    if not 1<=a.convective_iterations<=200:
        raise ValueError("invalid convective refinement iteration count")
    initial_spec = None
    continuation_spec = None
    if a.continuation_models:
        continuation_spec=json.loads((a.continuation_models/'specification.json').read_text())
        if any(continuation_spec[k]!=spec[k] for k in ['hydrogen','helium3','metals']):
            raise ValueError('continuation composition axes differ')
    if a.initial_models:
        initial_spec=json.loads((a.initial_models/"specification.json").read_text())
        if any(initial_spec[k] != spec[k] for k in ["hydrogen","helium3","teff_K","log_g","metals"]):
            raise ValueError("initial family axes/composition do not match")
    for name in ["synspec","tlusty"]:
        if digest(prepared[name]) != prepared["executables"][name]: raise ValueError("source executable changed")
    if prepared.get("opacity_method") != spec.get("opacity_method"):
        raise ValueError("opacity method does not match source preparation")
    if len(prepared["line_lists"]) != 4 or len(prepared["line_list_sha256"]) != 4:
        raise ValueError("incomplete source line lists")
    for path, expected in zip(prepared["line_lists"],prepared["line_list_sha256"]):
        if digest(path) != expected: raise ValueError("source line list changed")
    if data_digest(Path(prepared["synple"])/"data") != prepared["data_sha256"]:
        raise ValueError("source continuum/chemistry data changed")
    if not 1 <= a.jobs <= 8: raise ValueError("jobs must be 1..8")
    if not 2 <= spec["log_temperature"][0] <= 21 or not 2 <= spec["log_density"][0] <= 19:
        raise ValueError("opacity axes exceed compiled source capacity")
    if not 20 <= spec["depths"] <= 400 or not 2 <= spec["atmosphere_frequencies"] <= 32000:
        raise ValueError("atmosphere exceeds compiled source capacity")
    if not 1 < spec.get("temperature_step_limit", 1.03) <= 1.25:
        raise ValueError("invalid temperature correction limiter")
    if not 0 < spec.get("convection_derivative_step",.001) <= .01:
        raise ValueError("invalid convection derivative step")
    if "initial_depths" in spec and not (20 <= spec["initial_depths"] <= spec["depths"]
            and 2 <= spec["initial_frequencies"] <= spec["atmosphere_frequencies"]):
        raise ValueError("invalid initial atmosphere resolution")
    if not spec["atmosphere_frequencies"] <= spec["opacity_frequencies"] <= 64000:
        raise ValueError("invalid spectral resolution")
    if not 0 < spec["tau_top"] < .01 < spec["tau"] < spec["tau_bottom"]:
        raise ValueError("invalid atmosphere matching depths")
    if not 0 < spec["alpha"] < 10 or spec["wavelength_A"][0] <= 0 or spec["wavelength_A"][1] <= spec["wavelength_A"][0]:
        raise ValueError("invalid convection or wavelength parameters")
    root=a.work.resolve(); root.mkdir(parents=True,exist_ok=True)
    for filename,value in [("specification.json",spec),("provenance.json",prepared)]:
        existing=root/filename
        if existing.exists() and json.loads(existing.read_text()) != value:
            raise ValueError("source/settings changed: use a new work directory")
        existing.write_text(json.dumps(value,indent=2)+"\n")
    opacity={"temperature_K":temperatures(spec),"density_g_cm3":sequence(spec["log_density"])}
    planes=list(itertools.product(spec["hydrogen"],spec["helium3"]))
    if a.plane is not None and not 0 <= a.plane < len(planes): raise ValueError("invalid composition plane")
    if any(not 0 <= n < len(planes) for n in a.convective_initialization):
        raise ValueError("invalid convective initialization plane")

    def calculate_opacity(item):
        n,(x,y)=item; directory=root/f"plane-{n:03d}"
        tabledir=directory/"opacity"
        abundance,masses=composition(x,y,spec["metals"])
        print(f"computing opacity XH={x:g} X3={y:g}",flush=True)
        # SYNSPEC's binary output is unit 63, regardless of the text-table
        # filename in fort.2. Retain the original sequential record format.
        # Independent isotherms retain completed work if another source
        # state fails. They also avoid hidden cross-temperature source caches.
        isotherms=[]
        for i,t in enumerate(opacity["temperature_K"]):
            isotherms.append(tabledir/f"temperature-{i:03d}")
        # Exercise hot line-profile branches before the long molecular rows.
        for i in reversed(range(len(isotherms))):
            d=isotherms[i]; t=opacity["temperature_K"][i]
            opacity_inputs(d,prepared,spec,x,y,t)
            if a.reuse_opacity and not (d/"completed.json").exists():
                candidates=(a.reuse_opacity/f"plane-{n:03d}"/"opacity").glob("temperature-*")
                old=next((v for v in candidates if (v/"completed.json").exists() and
                          (v/"fort.2").read_bytes()==(d/"fort.2").read_bytes()),None)
                if old is not None:
                    for filename in ["fort.63","fort.29","run.log","completed.json"]:
                        shutil.copy2(old/filename,d/filename)
            # execute rechecks both the new input fingerprint and every
            # copied output hash; a mismatched isotherm is recomputed.
            execute(prepared["synspec"],d,["fort.63","fort.29"])
            validate_table(d/"fort.63",abundance,[t],opacity["density_g_cm3"])
            print(f"opacity XH={x:g} X3={y:g}: T={t:g} complete",flush=True)
        merge_isotherms([d/"fort.63" for d in isotherms],tabledir/"fort.63")
        validate_table(tabledir/"fort.63", abundance, opacity["temperature_K"], opacity["density_g_cm3"])
        return n,x,y,directory,tabledir,abundance,masses

    def calculate_model(job):
        plane,it,teff,ig,logg=job
        n,x,y,directory,tabledir,abundance,masses=plane
        d=directory/f"model-{it:03d}-{ig:03d}"
        if (d/"validated.json").exists() and (d/"completed.json").exists():
            saved=json.loads((d/"completed.json").read_text())
            record=json.loads((d/"validated.json").read_text())
            kinds=["log","convergence","atmosphere_input","element_masses","parameters"]
            if "initial_structure" in record: kinds.append("initial_structure")
            intact=(input_fingerprint(prepared["executables"]["tlusty"],d)==saved["input_sha256"]
                    and all((d/f).exists() and digest(d/f)==h for f,h in saved["outputs"].items())
                    and all((root/record[k]).exists() and digest(root/record[k])==record[k+"_sha256"] for k in kinds)
                    and (d/"opacity.sha256").read_text().strip()==digest(tabledir/"fort.63"))
            if intact:
                # A changed starting-guess choice must not discard a final
                # solution already validated with these source settings.
                from import_nongrey_grid import read_text
                inputs={k:read_text(root/record[k]) for k in kinds if k not in ["log","convergence"]}
                log=read_text(root/record["log"])
                source_inputs(inputs,spec,x,y,teff,logg,log)
                state=source_state(log,read_text(root/record["convergence"]),teff,logg,opacity,spec["tau"])
                if state["depths"] != spec["depths"] or any(record[k] != v for k,v in
                        zip(["XH","X3","teff_K","log_g"],[x,y,teff,logg])):
                    raise ValueError("cached atmosphere label or depth mismatch")
                (d/"failure.json").unlink(missing_ok=True)
                print(f"reused validated XH={x:g} X3={y:g} Teff={teff:g} logg={logg:g}",flush=True)
                return record
        (d/"failure.json").unlink(missing_ok=True)
        (d/"validated.json").unlink(missing_ok=True)
        initial = None
        continuation = None
        if a.continuation_models:
            candidates=[]
            for file in (a.continuation_models/f'plane-{n:03d}').glob('model-*/validated.json'):
                candidate=json.loads(file.read_text())
                if candidate['XH']!=x or candidate['X3']!=y:
                    raise ValueError('continuation composition label mismatch')
                distance=abs(math.log(teff/candidate['teff_K']))+.15*abs(logg-candidate['log_g'])*math.log(10)
                candidates.append((distance,str(file),candidate))
            if candidates:
                _,file,candidate=min(candidates,key=lambda v:(v[0],v[1]))
                old=Path(file).parent
                text=completed_initial_structure(old,continuation_spec,candidate['teff_K'],candidate['log_g'],spec['depths'])
                if text is not None:
                    initial=continuation_structure(text,candidate['teff_K'],candidate['log_g'],teff,logg)
                    continuation={'source':str(old.resolve()),'source_receipt_sha256':digest(old/'completed.json'),
                                  'teff_K':candidate['teff_K'],'log_g':candidate['log_g'],
                                  'method':'pressure-preserving temperature/gravity continuation; starting guess only'}
        if initial is None and a.initial_models:
            old=a.initial_models/f"plane-{n:03d}"/d.name
            if (old/"completed.json").exists() or (old/"completed.json.gz").exists():
                initial=completed_initial_structure(old,initial_spec,teff,logg,spec["depths"])
            elif (old/"checkpoint.json").exists():
                initial=checkpoint_initial_structure(old,x,y,teff,logg)
        if initial is None and "initial_depths" in spec:
            coarse = {**spec,"depths":spec["initial_depths"],
                      "atmosphere_frequencies":spec["initial_frequencies"]}
            seed = directory/f"seed-{it:03d}-{ig:03d}"
            atmosphere_inputs(seed,prepared,coarse,tabledir/"fort.63",abundance,masses,teff,logg)
            execute(prepared["tlusty"],seed,["fort.7","fort.9"])
            # A coarser solution is a starting guess, never a grid
            # cell. The final model still has the strict 0.2% test.
            source_state((seed/"run.log").read_text(),(seed/"fort.9").read_text(),
                         teff,logg,opacity,spec["tau"],max_flux_error=.05)
            initial = (seed/"fort.7").read_text()
        if n in a.convective_initialization:
            if initial is None:raise ValueError("convective initialization requires a completed starting structure")
            if a.convective_tail_depths:initial=convective_tail(initial,a.convective_tail_depths)
            precondition=directory/f"convective-initial-{it:03d}-{ig:03d}"
            atmosphere_inputs(precondition,prepared,spec,tabledir/"fort.63",abundance,masses,teff,logg,initial)
            tas=precondition/"tas"
            tas.write_text(tas.read_text()+f"ICONRE={a.convective_iterations},ICONRS=1,IMUCON=200,IDEEPC=3,CRFLIM=-1\n")
            executable=initializer['executable'] if initializer else prepared['tlusty']
            execute(executable,precondition,["fort.7","fort.9"])
            source_state((precondition/"run.log").read_text(),(precondition/"fort.9").read_text(),
                         teff,logg,opacity,spec["tau"])
            initial=(precondition/"fort.7").read_text()
        atmosphere_inputs(d,prepared,spec,tabledir/"fort.63",abundance,masses,teff,logg,initial)
        execute(prepared["tlusty"],d,["fort.7","fort.9"])
        state=source_state((d/"run.log").read_text(),(d/"fort.9").read_text(),teff,logg,opacity,spec["tau"])
        if state["depths"] != spec["depths"]:
            raise ValueError("final atmosphere depth count mismatch")
        record={"XH":x,"X3":y,"teff_K":teff,"log_g":logg,"diagnostics":state}
        if continuation is not None:record['continuation']=continuation
        if n in a.convective_initialization:
            files={}
            for filename in ['fort.5','fort.7','fort.8','fort.9','fort.15','tas','ember-masses.dat','opacity.sha256','physics.json','completed.json','run.log']:
                path=archive(precondition/filename,precondition)
                files[str(path.relative_to(root))]=digest(path)
            record['initialization']={'method':f'native CONREF first {a.convective_iterations} iterations, followed by independent canonical solve',
                'convective_tail_depths':a.convective_tail_depths,'source':initializer,
                'executable_sha256':digest(executable),'files_sha256':files}
        files = [("log","run.log"),("convergence","fort.9"),
                              ("atmosphere_input","fort.5"),("element_masses","ember-masses.dat"),
                              ("parameters","tas")]
        if initial is not None: files.append(("initial_structure","fort.8"))
        source_inputs({kind:(d/filename).read_text() for kind,filename in files if kind not in ["log","convergence"]},
                      spec,x,y,teff,logg,(d/"run.log").read_text())
        for name,filename in files:
            path=archive(d/filename,d); record[name]=str(path.relative_to(root));record[name+"_sha256"]=digest(path)
        (d/"validated.json").write_text(json.dumps(record,indent=2)+"\n")
        print(f"validated XH={x:g} X3={y:g} Teff={teff:g} logg={logg:g}",flush=True)
        return record

    def run_model(job):
        try:
            calculate_model(job)
            return True
        except Exception as error:
            plane,it,teff,ig,logg=job
            n,x,y,directory,*_=plane
            d=directory/f"model-{it:03d}-{ig:03d}"
            d.mkdir(parents=True,exist_ok=True)
            failure={"XH":x,"X3":y,"teff_K":teff,"log_g":logg,"error":str(error)}
            (d/"failure.json").write_text(json.dumps(failure,indent=2)+"\n")
            print(f"REJECTED XH={x:g} X3={y:g} Teff={teff:g} logg={logg:g}: {error}",flush=True)
            return False

    selected=list(enumerate(planes)) if a.plane is None else [(a.plane,planes[a.plane])]
    # Opacity jobs retain large molecular line lists; limit those to four.
    with ThreadPoolExecutor(max_workers=min(a.jobs,4)) as pool:
        ready=list(pool.map(calculate_opacity,selected))
    if a.opacity_only: return
    jobs=[(plane,it,teff,ig,logg) for it,teff in enumerate(spec["teff_K"])
          for ig,logg in enumerate(spec["log_g"]) for plane in ready]
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        accepted=list(pool.map(run_model,jobs))
    models=[json.loads(f.read_text()) for f in sorted(root.glob("plane-*/model-*/validated.json"))]
    manifest={**spec,"format":1,"calculation":CALCULATION,"provenance":prepared,"opacity":opacity,"models":models}
    (root/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    if not all(accepted):
        raise RuntimeError(f"{accepted.count(False)} rejected atmosphere models; see model-*/failure.json")
    if a.plane is None: import_grid(root/"manifest.json",root/"atmosphere.dat")


if __name__ == "__main__":
    main()
