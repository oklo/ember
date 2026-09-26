#!/usr/bin/env python3
"""Retain immutable research inputs outside temporary directories.

Verified content is hard-linked into the workspace on the same filesystem,
so temporary-path removal cannot discard the last copy and large source
tables are not duplicated. Use these links only for immutable inputs.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
from write_scientific_result import write_result


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('plan', type=Path)
    ap.add_argument('store', type=Path)
    ap.add_argument('report', type=Path)
    a = ap.parse_args()
    if a.report.exists():
        raise FileExistsError('preserve the existing retention report')
    plan = json.loads(a.plan.read_text())
    a.store.mkdir(parents=True, exist_ok=True)
    entries = []; failures = []; linked_bytes = 0; copied_bytes = 0
    for name, expected in plan['input_sha256'].items():
        try:
            source = Path(name).resolve(strict=True)
            actual = digest(source)
            if actual != expected:
                raise ValueError('input differs from the recorded checksum')
            destination = a.store/actual[:2]/actual
            destination.parent.mkdir(exist_ok=True)
            if destination.exists():
                if digest(destination) != actual:
                    raise ValueError('retained content is corrupt')
                method = 'existing'
            elif source.stat().st_dev == destination.parent.stat().st_dev:
                os.link(source, destination); method = 'hard_link'
                linked_bytes += source.stat().st_size
            else:
                with source.open('rb') as src, destination.open('xb') as dst:
                    shutil.copyfileobj(src, dst)
                method = 'copy'; copied_bytes += source.stat().st_size
            if digest(destination) != actual:
                raise ValueError('retained content failed verification')
            entries.append(dict(original=str(source), retained=str(destination.resolve()),
                                sha256=actual, bytes=source.stat().st_size, method=method))
        except (OSError, ValueError) as e:
            failures.append(dict(path=name, error=str(e)))
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  outcome='retained_verified_inputs' if not failures else 'incomplete_retention',
                  plan_sha256=digest(a.plan), entries=entries, failures=failures,
                  newly_linked_bytes=linked_bytes, copied_bytes=copied_bytes,
                  limitations=['Hard links protect against deletion of a temporary pathname, not modification of the shared inode.',
                               'Only immutable source inputs belong in this store; active logs and checkpoints are separate artifacts.'])
    write_result(a.report, report)
    print(len(entries), 'inputs retained;', len(failures), 'failures;', copied_bytes, 'copied bytes')
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
