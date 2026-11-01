"""Consolidate checked, byte-identical physics inputs; preserve old paths.

The plan is written separately from execution. Only archived material files are
eligible. Logs, restart states, source code, environments and Git are excluded.
"""
import argparse
from collections import defaultdict
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
STORE = ROOT / 'data/retained/sha256'


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def eligible(name):
    p = Path(name)
    if p.parts[0] not in ('data', 'out'):
        return False
    if any(s in p.parts for s in ('.git', 'condaenv', 'site-packages', 'pkgs')):
        return False
    return (p.suffix in ('.dat', '.bin', '.hdf5') or
            p.name == 'source.json.gz' or p.name.endswith('_raw.json.gz'))


def plan(plan_path, archives=False):
    db = sqlite3.connect(HERE / 'inventory.sqlite')
    groups = defaultdict(list)
    for name, size, allocated, mtime, inode, checksum in db.execute('''
        select f.path,f.bytes,f.allocated,f.mtime,f.inode,h.sha256
        from files f join hashes h using(path) where h.sha256 is not null'''):
        path = ROOT/name
        additional = (not eligible(name) and Path(name).parts[0] == 'out' and
                      not any(s in Path(name).parts for s in ('.git','condaenv','site-packages','pkgs','metallicity-map-oct27')) and
                      (name.endswith(('.csv','.csv.gz','.json.gz','.npz')) or Path(name).name=='fort.63'))
        if (additional if archives else eligible(name)) and path.is_file() and not path.is_symlink():
            groups[checksum].append(dict(path=name, bytes=size,
                allocated=allocated, mtime_ns=mtime, inode=inode))
    actions = []
    for checksum, files in groups.items():
        # Already shared hard links consume one allocation; leave those alone.
        inodes = {f['inode']: f['allocated'] for f in files}
        if len(inodes) < 2:
            continue
        actions.append(dict(sha256=checksum, files=files,
            store=str((STORE/checksum[:2]/checksum).relative_to(ROOT)),
            estimated_reclaimed_bytes=sum(inodes.values())-max(inodes.values())))
    with plan_path.open('x') as stream:
        for action in actions:
            stream.write(json.dumps(action)+'\n')
    print({'groups': len(actions), 'paths': sum(len(a['files']) for a in actions),
           'estimated_GiB': sum(a['estimated_reclaimed_bytes'] for a in actions)/2**30})


def apply(plan_path, receipt_path):
    last = time.monotonic()
    with plan_path.open() as plan_file, receipt_path.open('a', buffering=1) as receipts:
        for index, line in enumerate(plan_file, 1):
            action = json.loads(line)
            target = ROOT/action['store']
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.exists():
                source = ROOT/action['files'][0]['path']
                if digest(source) != action['sha256']:
                    raise ValueError(f'changed source: {source}')
                os.link(source, target)
            if digest(target) != action['sha256']:
                raise ValueError(f'corrupt retained input: {target}')
            target.chmod(0o444)
            for entry in action['files']:
                path = ROOT/entry['path']
                if path.is_symlink() and path.resolve() == target:
                    continue  # A previously completed operation is resumable.
                stat = path.lstat()
                if (path.resolve() != path or stat.st_size != entry['bytes'] or
                        stat.st_mtime_ns != entry['mtime_ns'] or stat.st_ino != entry['inode']):
                    raise ValueError(f'changed inventory entry: {path}')
                if digest(path) != action['sha256']:
                    raise ValueError(f'changed contents: {path}')
                # Create the replacement first; os.replace never leaves a gap.
                temporary = path.with_name(path.name+'.consolidating')
                temporary.symlink_to(os.path.relpath(target, path.parent))
                os.replace(temporary, path)
                receipts.write(json.dumps(dict(path=entry['path'], sha256=action['sha256'],
                    retained=action['store'], bytes=entry['bytes'],
                    action='duplicate replaced by relative link to read-only retained input'))+'\n')
            if time.monotonic()-last > 20:
                print({'groups_completed': index}, flush=True)
                last = time.monotonic()
    print('Input consolidation complete.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=('plan', 'apply'))
    parser.add_argument('--archives', action='store_true', help='Include only older tabular/gzip archives and opacity scratch tables.')
    parser.add_argument('--record-prefix', default='input-consolidation')
    args = parser.parse_args()
    plan_path=HERE/(args.record_prefix+'-plan.jsonl')
    receipt_path=HERE/(args.record_prefix+'-receipts.jsonl')
    plan(plan_path,args.archives) if args.operation == 'plan' else apply(plan_path,receipt_path)
