#!/usr/bin/env python3
"""Copy a lifetime configuration's complete table set into a portable directory.

Numerical tables are copied byte for byte. Only filenames in configurations
and family manifests change. The original inputs are never modified.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil


INPUT_KEYS = ("eos", "opacity_low", "opacity_warm", "opacity_bridge",
              "opacity_hot", "conduction", "atmosphere", "collisions", "composition")
OPTIONAL_INPUT_KEYS = ("atmosphere_main_sequence",)


def input_keys(configuration):
    return INPUT_KEYS + tuple(key for key in OPTIONAL_INPUT_KEYS if key in configuration)


def is_family(path, role):
    if role == "eos":
        with path.open("rb") as stream:
            header = stream.readline(128)
        return not header.startswith(b"EMBER_VARIABLE_METAL_HELMHOLTZ_BINARY ")
    return role.startswith("opacity_")


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def settings(path):
    result = {}
    for line in path.read_text().splitlines():
        row = shlex.split(line, comments=True)
        if not row:
            continue
        if len(row) != 2 or row[0] in result:
            raise ValueError(f"invalid or duplicate configuration entry: {line}")
        result[row[0]] = row[1]
    return result


def family(path, role):
    """Return header lines and (numeric prefix, child filename) entries."""
    lines = path.read_text().splitlines()
    if role == "eos":
        if lines[0] != "EMBER_VARIABLE_METAL_HELMHOLTZ 2":
            raise ValueError(f"unsupported EOS manifest: {path}")
        count = 1
        for line, name in zip(lines[1:4], ("metals", "hydrogen_share", "helium3_share")):
            row = shlex.split(line)
            if row[0] != name or len(row) != int(row[1]) + 2:
                raise ValueError(f"invalid EOS axis: {path}")
            count *= int(row[1])
        header, rows = lines[:4], [shlex.split(line) for line in lines[4:]]
        width = 1
    else:
        row = shlex.split(lines[0])
        if len(row) != 5 or row[:2] != ["EMBER_OPACITY_MIXTURE", "1"]:
            raise ValueError(f"unsupported opacity manifest: {path}")
        count = int(row[2])
        header, rows = lines[:1], [shlex.split(line) for line in lines[1:]]
        width = 2
    rows = [row for row in rows if row]
    if len(rows) != count or any(len(row) != width for row in rows):
        raise ValueError(f"wrong family entry count: {path}")
    return header, [(row[:-1], (path.parent / row[-1]).resolve(strict=True)) for row in rows]


def write_checked(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != content:
            raise ValueError(f"refusing to replace different packaged input: {path}")
    else:
        with path.open("xb") as stream:
            stream.write(content)


def verify(manifest):
    record = json.loads(manifest.read_text())
    root = manifest.parent.resolve()
    for entry in record["files"]:
        path = (root / entry["path"]).resolve(strict=True)
        if not path.is_relative_to(root):
            raise ValueError(f"packaged path escapes data directory: {path}")
        if path.stat().st_size != entry["size"] or digest(path) != entry["sha256"]:
            raise ValueError(f"packaged input differs: {path}")
    config = (root / record["configuration"]).resolve(strict=True)
    cfg = settings(config)
    discovered = {config}
    for role in input_keys(cfg):
        path = (config.parent / cfg[role]).resolve(strict=True)
        discovered.add(path)
        if is_family(path, role):
            _, children = family(path, role)
            discovered.update(child for _, child in children)
    declared = {(root / entry["path"]).resolve() for entry in record["files"]}
    if discovered != declared:
        raise ValueError("manifest does not match the complete runtime input set")
    print(f"Verified {len(declared)} packaged files; no external runtime table paths.")


def package(config, destination, metadata, manifest, source_root):
    config = config.resolve(strict=True)
    destination = destination.resolve()
    manifest = manifest.resolve()
    source_root = source_root.resolve()
    if not destination.is_relative_to(manifest.parent):
        raise ValueError("package must be inside the manifest's data directory")
    provenance = json.loads(metadata.read_text())
    cfg = settings(config)
    copied = {}
    entries = []

    def record(source, target, role, original_hash):
        info = provenance[role]
        for key in ("generator", "licence"):
            if not info.get(key):
                raise ValueError(f"missing {key} for {role}")
        entries.append(dict(path=target.relative_to(manifest.parent).as_posix(),
                            size=target.stat().st_size, sha256=digest(target),
                            source=Path(os.path.relpath(source, source_root)).as_posix(),
                            source_sha256=original_hash, role=role, **info))

    def leaf(source, role):
        if source in copied:
            return copied[source]
        checksum = digest(source)
        target = destination / "tables" / f"{checksum}.dat"
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            if digest(target) != checksum:
                raise ValueError(f"corrupt existing table: {target}")
        else:
            with source.open("rb") as src, target.open("xb") as dst:
                shutil.copyfileobj(src, dst)
        if digest(target) != checksum:
            raise ValueError(f"copy verification failed: {target}")
        copied[source] = target
        record(source, target, role, checksum)
        return target

    for role in input_keys(cfg):
        source = (config.parent / cfg[role]).resolve(strict=True)
        if is_family(source, role):
            header, children = family(source, role)
            target = destination / f"{role}.dat"
            rows = header[:]
            for prefix, child in children:
                local = leaf(child, role).relative_to(target.parent).as_posix()
                rows.append(" ".join(prefix + [json.dumps(local)]))
            write_checked(target, ("\n".join(rows) + "\n").encode())
            copied[source] = target
            record(source, target, role, digest(source))
        else:
            target = leaf(source, role)
        cfg[role] = target.relative_to(destination).as_posix()
    target = destination / "configuration.txt"
    write_checked(target, ("".join(f"{key} {json.dumps(value)}\n" for key, value in cfg.items())).encode())
    record(config, target, "configuration", digest(config))
    result = dict(format=1, configuration=target.relative_to(manifest.parent).as_posix(),
                  original_input_count=len(copied),
                  original_input_bytes=sum(path.stat().st_size for path in copied),
                  limitations=provenance.get("limitations", []),
                  files=sorted(entries, key=lambda row: (row["path"], row["source"])))
    write_checked(manifest, (json.dumps(result, indent=2, sort_keys=True) + "\n").encode())
    verify(manifest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--destination", type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--source-root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    if args.verify:
        verify(args.verify)
    elif all((args.config, args.destination, args.metadata, args.manifest)):
        package(args.config, args.destination, args.metadata, args.manifest, args.source_root)
    else:
        parser.error("supply --verify, or --config, --destination, --metadata and --manifest")


if __name__ == "__main__":
    main()
