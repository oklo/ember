#!/usr/bin/env python3
"""Archive a complete, accepted family without its large opacity binaries.

The archive retains original atmosphere inputs/outputs, source receipts,
opacity metadata and hashes, and the pinned generation recipe. It can be
reimported offline; recomputing the spectra requires the source preparation.
"""
import argparse
import gzip
import itertools
import json
from pathlib import Path
import shutil
import tempfile

from generate_nongrey_grid import composition, temperatures, sequence, input_fingerprint
from import_nongrey_grid import import_grid
from nongrey_opacity import validate_table
from prepare_nongrey_sources import digest


def archive_family(work, destination):
    work = Path(work).resolve(); destination = Path(destination).resolve()
    if destination.exists():
        raise FileExistsError(destination)
    manifest = json.loads((work / "manifest.json").read_text())
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent) as temporary:
        root = Path(temporary) / "archive"; root.mkdir()
        import_grid(work / "manifest.json", root / "atmosphere.dat")
        files = {}

        def copy(relative, compress=False):
            source = work / relative
            target = root / (str(relative) + (".gz" if compress else ""))
            target.parent.mkdir(parents=True, exist_ok=True)
            if compress:
                target.write_bytes(gzip.compress(source.read_bytes(), compresslevel=9, mtime=0))
            else:
                shutil.copyfile(source, target)
            files[str(target.relative_to(root))] = digest(target)

        for name in ["specification.json", "provenance.json"]:
            copy(name)
        opacity_planes = []
        axes = itertools.product(manifest["hydrogen"], manifest["helium3"])
        for i, (x, y) in enumerate(axes):
            relative = Path(f"plane-{i:03d}/opacity")
            source = work / relative / "fort.63"
            abundance, _ = composition(x, y, manifest["metals"])
            table = validate_table(source, abundance, temperatures(manifest),
                                   sequence(manifest["log_density"]))
            opacity_planes.append({"XH": x, "X3": y, "sha256": digest(source),
                                   "shape_frequency_density_temperature": table["shape"]})
            for isotherm in sorted((work / relative).glob("temperature-*")):
                receipt = json.loads((isotherm / "completed.json").read_text())
                if input_fingerprint(manifest["provenance"]["executables"]["synspec"], isotherm) != receipt["input_sha256"]:
                    raise ValueError("opacity input fingerprint mismatch")
                for name, expected in receipt["outputs"].items():
                    if digest(isotherm / name) != expected:
                        raise ValueError("opacity receipt checksum mismatch")
                for name in ["fort.5", "fort.55", "fort.2", "tas", "ember-masses.dat",
                             "physics.json", "fort.29", "run.log", "completed.json"]:
                    copy((isotherm / name).relative_to(work), compress=True)
        for model in manifest["models"]:
            directory = Path(model["log"]).parent
            if 'initialization' in model:
                for name,expected in model['initialization']['files_sha256'].items():
                    if digest(work/name)!=expected:raise ValueError('initialization archive checksum mismatch')
                    copy(name)
            loaded_opacity = (work / directory / "opacity.sha256").read_text().strip()
            plane = next(p for p in opacity_planes if (p["XH"], p["X3"]) ==
                         (model["XH"], model["X3"]))
            if loaded_opacity != plane["sha256"]:
                raise ValueError("atmosphere used a different opacity table")
            receipt = json.loads((work / directory / "completed.json").read_text())
            if input_fingerprint(manifest["provenance"]["executables"]["tlusty"], work / directory) != receipt["input_sha256"]:
                raise ValueError("atmosphere input fingerprint mismatch")
            for name, expected in receipt["outputs"].items():
                if digest(work / directory / name) != expected:
                    raise ValueError("atmosphere receipt checksum mismatch")
            for kind in ["log", "convergence", "atmosphere_input", "element_masses",
                         "parameters", "initial_structure"]:
                if kind in model:
                    copy(model[kind])
            for name in ["completed.json", "fort.7", "fort.15", "opacity.sha256", "physics.json"]:
                copy(directory / name, compress=True)
        manifest["opacity_planes"] = opacity_planes
        manifest["archive_files_sha256"] = files
        (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        import_grid(root / "manifest.json", root / "reimport.dat")
        if (root / "atmosphere.dat").read_bytes() != (root / "reimport.dat").read_bytes():
            raise ValueError("archive does not reproduce the accepted grid")
        (root / "reimport.dat").unlink()
        root.rename(destination)
    return destination / "atmosphere.dat"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    print(archive_family(args.work, args.destination))


if __name__ == "__main__":
    main()
