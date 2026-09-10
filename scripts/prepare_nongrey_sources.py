#!/usr/bin/env python3
"""Fetch pinned TLUSTY/SYNSPEC inputs and build optional atmosphere tools.

Normal Ember builds do not run this script or need Fortran/network. All work
goes below the explicitly supplied cache directory. Original archives are
retained. The public FTP login is published in Synple's linelists/makefile.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from ftplib import FTP
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ROOT / "data/atmosphere/sources"
ASSETS = {
    "tl208-s54.tar.gz": (
        "https://www.as.arizona.edu/~hubeny/tlusty208-package/tl208-s54.tar.gz",
        "ec9febdc1795f2c1bbe948ea1736bed663aac9238b290483b517201132f85e4a"),
    "synple.tar.gz": (
        "https://codeload.github.com/callendeprieto/synple/tar.gz/703549b60031c5f6dd2446ad933b129c070a6560",
        "7e98b2f20f5c2c8ee399769f63a2d511cf69fcf2357975d94dcb798d1d5f7967"),
    "gfATOc.19.gz": ("ftp:gfATOc.19.gz", "0fb5f9945bce3defc151da7331d0daf96136e360eda1703b014a044b2839e737"),
    "gfMOLsun.20.gz": ("ftp:gfMOLsun.20.gz", "68491ac8068cfbf86491b639cb4a9a2461a8807c865b8b15b2c08e0b297a2586"),
    "gfTiO.20.gz": ("ftp:gfTiO.20.gz", "cee8cf61473303f5aca8f93f90cd0c404fb8fbcc4c40df9283c54edcb193a456"),
    "H2O-8.20.gz": ("ftp:H2O-8.20.gz", "8ca3716c74c8cd5e00905b3134f02bb164955cb1be97379174897766075aaf29"),
    "optables.tar.gz": (
        "https://www.as.arizona.edu/~hubeny/tlusty208-package/optables.tar.gz",
        "f0fa3faa6c9f5ed45a75f235629e91cbb71b4a067feaf4abd455eca0ffdef57a"),
}
LINES = ["gfATOc.19", "gfMOLsun.20", "gfTiO.20", "H2O-8.20"]


def digest(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def data_digest(directory):
    root = Path(directory)
    result = hashlib.sha256()
    for path in sorted(root.rglob("*")):
        if path.is_file():
            result.update(str(path.relative_to(root)).encode()+b"\0"+digest(path).encode()+b"\0")
    return result.hexdigest()


def fetch(directory, name):
    path = directory / name
    url, expected = ASSETS[name]
    if not path.exists():
        temporary = path.with_suffix(path.suffix + ".part")
        with temporary.open("wb") as out:
            if url.startswith("ftp:"):
                with FTP("ftp.ll.iac.es", timeout=120) as ftp:
                    ftp.login("carlos", "allende")
                    ftp.cwd("linelists")
                    ftp.retrbinary("RETR " + url[4:], out.write, 1024 * 1024)
            else:
                with urllib.request.urlopen(url, timeout=180) as src:
                    shutil.copyfileobj(src, out)
        if digest(temporary) != expected:
            raise ValueError(f"download checksum mismatch: {name}")
        temporary.replace(path)
    if digest(path) != expected:
        raise ValueError(f"cached source checksum mismatch: {name}")
    return path


def extract(archive, directory):
    directory.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive) as source:
        # Upstream contains absolute symlinks to the author's external disk.
        # Install only regular files/directories; link verified data explicitly.
        members = [m for m in source if m.isfile() or m.isdir()]
        source.extractall(directory, members=members, filter="data")
        return directory / members[0].name.split("/")[0]


def run(args, cwd, log, stdin=None):
    with Path(log).open("wb") as out:
        result = subprocess.run(args, cwd=cwd, stdin=stdin, stdout=out, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"source command failed ({result.returncode}); see {log}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cache", type=Path)
    p.add_argument("--compiler", default="gfortran")
    p.add_argument("--offline", action="store_true")
    a = p.parse_args()
    root = a.cache.resolve(); downloads = root / "downloads"
    downloads.mkdir(parents=True, exist_ok=True)
    for name in ASSETS:
        if a.offline and not (downloads / name).exists():
            raise FileNotFoundError(f"offline source missing: {downloads / name}")
        fetch(downloads, name)
    tl = extract(downloads / "tl208-s54.tar.gz", root / "tlusty-source")
    sy = extract(downloads / "synple.tar.gz", root / "synple-source")
    tables = extract(downloads / "optables.tar.gz", root / "benchmark-source")
    for name, directory in [("tlusty208", tl), ("synspec54", sy)]:
        with (SOURCES / (name + "-ember.patch")).open("rb") as patch:
            run(["patch", "-p1", "--batch"], directory, root / (name + "-patch.log"), patch)
    flags = [a.compiler, "-O2", "-g", "-fno-automatic", "-std=legacy", "-fallow-argument-mismatch"]
    run(flags + ["-fcheck=bounds", "-fbacktrace", "-o", "tlusty.exe", "tlusty208.f"], tl / "tlusty", root / "tlusty-build.log")
    run(flags + ["-fcheck=bounds", "-fbacktrace", "-o", "synspec54", "synspec54.f"],
        sy / "synspec", root / "synspec-build.log")
    sampled = root / "sampling-source"
    shutil.copytree(sy / "synspec", sampled / "synspec", dirs_exist_ok=True)
    with (SOURCES / "synspec54-sampling.patch").open("rb") as patch:
        run(["patch", "-p1", "--batch"], sampled, root / "sampling-patch.log", patch)
    run(flags + ["-fcheck=bounds", "-fbacktrace", "-o", "synspec54", "synspec54.f"],
        sampled / "synspec", root / "sampling-build.log")
    converter = root / "list2bin"
    run([a.compiler, "-O2", "-std=legacy", "-o", str(converter), "list2bin.f"],
        sy / "synspec", root / "list2bin-build.log")
    lines = root / "lines"; lines.mkdir(exist_ok=True)

    def convert(name):
        directory = lines / name; directory.mkdir(exist_ok=True)
        with gzip.open(downloads / (name + ".gz"), "rb") as f, (directory / "ascii").open("wb") as out:
            shutil.copyfileobj(f, out)
        with (directory / "ascii").open("rb") as f:
            run([str(converter)], directory, directory / "conversion.log", f)
        if "lines included" not in (directory / "conversion.log").read_text():
            raise ValueError("line conversion did not finish")
        return str(directory / "fort.12")

    with ThreadPoolExecutor(max_workers=4) as pool:
        binary = list(pool.map(convert, LINES))
    receipt = {
        "format": 1, "tlusty": str(tl / "tlusty/tlusty.exe"),
        "synspec": str(sampled / "synspec/synspec54"), "dense_synspec": str(sy / "synspec/synspec54"),
        "opacity_method": "sampling", "synple": str(sy),
        "tlusty_source": str(tl), "benchmark_tables": str(tables),
        "line_lists": binary,
        "inputs": {name: {"url": url, "sha256": sha} for name, (url, sha) in ASSETS.items()},
        "patches": {name: digest(SOURCES / name) for name in ["tlusty208-ember.patch", "synspec54-ember.patch", "synspec54-sampling.patch"]},
        "line_list_sha256": [digest(path) for path in binary],
        "data_sha256": data_digest(sy / "data"),
        "compiler": subprocess.check_output([a.compiler, "--version"], text=True).splitlines()[0],
    }
    receipt["executables"] = {name: digest(receipt[name]) for name in ["tlusty", "synspec", "dense_synspec"]}
    (root / "prepared.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(root / "prepared.json")


if __name__ == "__main__":
    main()
