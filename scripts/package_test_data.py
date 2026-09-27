#!/usr/bin/env python3
"""Verify, bundle or install the exact external inputs used by Ember's tests.

The committed manifest fixes every filename, size and SHA-256. Installation
verifies the entire archive before writing, refuses different existing files,
and never extracts links or filenames supplied only by the archive.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tarfile


def checksum(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--manifest', type=Path)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--verify', action='store_true')
    mode.add_argument('--pack', type=Path)
    mode.add_argument('--install', type=Path)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    manifest = args.manifest or root / 'data/TEST_DATA_MANIFEST.json'
    record = json.loads(manifest.read_text())
    if record['format'] != 1:
        raise ValueError('unsupported test-data manifest')
    files = {}
    for entry in record['files']:
        name = entry['path']
        p = Path(name)
        if p.is_absolute() or '..' in p.parts or p.parts[0] != 'data' or name in files:
            raise ValueError('invalid or duplicate test-data path')
        if not (root / p).resolve().is_relative_to(root):
            raise ValueError('test-data path leaves checkout')
        files[name] = entry

    def verify_file(name):
        p, entry = root / name, files[name]
        if p.stat().st_size != entry['size'] or checksum(p) != entry['sha256']:
            raise ValueError('test input differs: ' + name)

    if args.install:
        with tarfile.open(args.install, 'r:gz') as archive:
            members = archive.getmembers()
            if len(members) != len(files) or {m.name for m in members} != set(files):
                raise ValueError('archive does not match the declared test inputs')
            for member in members:
                entry = files[member.name]
                if not member.isfile() or member.size != entry['size']:
                    raise ValueError('invalid archive entry: ' + member.name)
                with archive.extractfile(member) as stream:
                    if hashlib.file_digest(stream, 'sha256').hexdigest() != entry['sha256']:
                        raise ValueError('corrupt archive entry: ' + member.name)
                p = root / member.name
                if p.exists():
                    verify_file(member.name)
            for member in members:
                p = root / member.name
                if p.exists():
                    continue
                p.parent.mkdir(parents=True, exist_ok=True)
                with archive.extractfile(member) as source, p.open('xb') as target:
                    shutil.copyfileobj(source, target)
                verify_file(member.name)
    else:
        for name in files:
            verify_file(name)
        if args.pack:
            with args.pack.open('xb') as output, tarfile.open(fileobj=output, mode='w:gz', compresslevel=3) as archive:
                for name in files:
                    p = root / name
                    info = tarfile.TarInfo(name)
                    info.size = files[name]['size']
                    info.mode = 0o644
                    with p.open('rb') as stream:
                        archive.addfile(info, stream)
    print(f"Verified {len(files)} test inputs ({sum(v['size'] for v in files.values())} bytes).")


if __name__ == '__main__':
    main()
