#!/usr/bin/env python3
"""Inventory local research artifacts without following links or changing files."""

import argparse
import csv
import gzip
import heapq
import json
import os
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import time


def inventory(root, destination):
    destination.mkdir(parents=True, exist_ok=False)
    groups = defaultdict(lambda: [0, 0, 0])
    names = defaultdict(lambda: [0, 0, 0])
    extensions = defaultdict(lambda: [0, 0, 0])
    hardlinks = set()
    largest = []
    errors = []
    totals = [0, 0, 0]
    last_update = time.monotonic()
    root_device = root.stat().st_dev
    excluded = {root / '.git', destination}

    def add(counter, size, allocated):
        counter[0] += 1
        counter[1] += size
        counter[2] += allocated

    with gzip.open(destination / 'files.csv.gz', 'wt', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['path', 'bytes', 'allocated_bytes', 'mtime_ns', 'device',
                         'inode', 'links', 'kind', 'symlink_target'])
        pending = [root]
        while pending:
            directory = pending.pop()
            try:
                entries = os.scandir(directory)
            except OSError as error:
                errors.append({'path': str(directory), 'error': str(error)})
                continue
            with entries:
                for entry in entries:
                    path = Path(entry.path)
                    if path in excluded:
                        continue
                    try:
                        stat = entry.stat(follow_symlinks=False)
                        if stat.st_dev != root_device:
                            errors.append({'path': str(path), 'error': 'different filesystem; skipped'})
                            continue
                        if entry.is_dir(follow_symlinks=False):
                            pending.append(path)
                            continue
                        relative = path.relative_to(root)
                        link = entry.is_symlink()
                        target = os.readlink(path) if link else ''
                        kind = 'symlink' if link else 'file' if entry.is_file(follow_symlinks=False) else 'special'
                        allocated = stat.st_blocks * 512
                        identity = (stat.st_dev, stat.st_ino)
                        unique_allocated = allocated
                        if stat.st_nlink > 1:
                            if identity in hardlinks:
                                unique_allocated = 0
                            hardlinks.add(identity)
                        writer.writerow([str(relative), stat.st_size, allocated, stat.st_mtime_ns,
                                         *identity, stat.st_nlink, kind, target])
                        add(totals, stat.st_size, unique_allocated)
                        for depth in [1, 2, 3]:
                            if len(relative.parts) > depth:
                                add(groups['/'.join(relative.parts[:depth])], stat.st_size, unique_allocated)
                        add(names[path.name], stat.st_size, unique_allocated)
                        add(extensions[''.join(path.suffixes[-2:]) or '(none)'], stat.st_size, unique_allocated)
                        item = (stat.st_size, str(relative), allocated)
                        if len(largest) < 200:
                            heapq.heappush(largest, item)
                        elif item > largest[0]:
                            heapq.heapreplace(largest, item)
                    except OSError as error:
                        errors.append({'path': str(path), 'error': str(error)})
            if time.monotonic() - last_update > 20:
                print(json.dumps({'files_scanned': totals[0], 'logical_GiB': totals[1] / 2**30}), flush=True)
                last_update = time.monotonic()

    def ranked(mapping, limit=None):
        items = sorted(mapping.items(), key=lambda item: item[1][2], reverse=True)
        return [{'name': name, 'files': values[0], 'bytes': values[1], 'allocated_bytes': values[2]}
                for name, values in items[:limit]]

    report = {
        'utc': datetime.now(timezone.utc).isoformat(),
        'root': str(root),
        'excluded': [str(path) for path in sorted(excluded)],
        'files': totals[0],
        'logical_bytes': totals[1],
        'unique_allocated_bytes': totals[2],
        'directories': ranked(groups),
        'basenames': ranked(names, 150),
        'extensions': ranked(extensions, 100),
        'largest_files': [{'path': path, 'bytes': size, 'allocated_bytes': allocated}
                          for size, path, allocated in sorted(largest, reverse=True)],
        'errors': errors,
        'measurement': 'lstat sizes; allocated bytes count each hardlinked inode once; symlinks are not followed; APFS shared clone extents are not measured separately',
    }
    (destination / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'destination': str(destination), 'files': totals[0],
                      'logical_GiB': totals[1] / 2**30,
                      'unique_allocated_GiB': totals[2] / 2**30,
                      'errors': len(errors)}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('.'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    inventory(args.root.resolve(), args.output.resolve())
