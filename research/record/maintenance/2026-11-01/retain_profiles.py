"""Keep representative interiors and remove reproducible full profile dumps.

Surface histories and restart states remain intact. The 0.1-Msun study, current
paused campaign, long integrations and runs without verifiable launch inputs are
excluded. Plan and apply are separate operations; each removal has a receipt.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import sqlite3
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def plan():
    db = sqlite3.connect(HERE/'inventory.sqlite')
    entries, excluded, checksum_cache = [], [], {}

    def checksum(path):
        key = str(path.resolve())
        if key not in checksum_cache:
            checksum_cache[key] = digest(path)
        return checksum_cache[key]

    for name, size, mtime in db.execute("select path,bytes,mtime from files where path like '%/interior_full.csv.gz'"):
        path = ROOT/name
        if size < 8*2**20 or 'metallicity-map-oct27' in path.parts:
            continue
        work = path.parent
        try:
            report = json.loads((work/'report.json').read_text())
            if not 0 < report['cpu_seconds'] <= 10800:
                raise ValueError('integration exceeds three CPU hours')
            launch = json.loads((work.parent/'launch.json').read_text())
            command = launch['command']
            if command[1] != '--lifetime':
                raise ValueError('unknown replay command')
            executable, configuration = Path(command[0]), Path(command[2])
            settings = dict(shlex.split(line) for line in configuration.read_text().splitlines()
                            if line.strip() and not line.startswith('#'))
            if abs(float(settings['mass_Msun'])-.1) < 1e-8:
                raise ValueError('preserve full 0.1-Msun interiors')
            for key, value in settings.items():
                if value.startswith('/') and not Path(value).is_file():
                    raise ValueError(f'missing configuration input: {key}')
            if checksum(executable) != launch['executable_sha256']:
                raise ValueError('executable differs from launch identity')
            if checksum(configuration) != launch['configuration_sha256']:
                raise ValueError('configuration differs from launch identity')
            restart = Path(command[command.index('--restart')+1]) if '--restart' in command else None
            if restart and checksum(restart) != launch['restart_sha256']:
                raise ValueError('restart differs from launch identity')
            retained = [work/n for n in ['history.jsonl', 'seed.checkpoint', 'final.checkpoint', 'interior.csv']]
            if not all(p.is_file() and p.stat().st_size for p in retained):
                raise ValueError('missing compact history or restart state')
            rows = [json.loads(line) for line in retained[0].read_text().splitlines() if line.strip()]
            if len(rows) < 100:
                raise ValueError('short history already compact')
            indices = {0, len(rows)-1}
            # Uniform samples preserve late evolution; log-time samples resolve
            # the early contraction. Always retain scalar extrema and changes
            # in the number of radiative boundaries, with adjacent models.
            for i in range(33):
                indices.add(round(i*(len(rows)-1)/32))
            positive = [(i, math.log10(r['years'])) for i,r in enumerate(rows) if r['years']>0]
            if positive:
                lo, hi = positive[0][1], positive[-1][1]
                for i in range(33):
                    target = lo+(hi-lo)*i/32
                    indices.add(min(positive, key=lambda x: abs(x[1]-target))[0])
            for field in ('Teff_K', 'luminosity_Lsun', 'central_He3', 'nuclear_fraction'):
                available = [i for i,r in enumerate(rows) if field in r]
                if available:
                    indices.add(min(available, key=lambda i: rows[i][field]))
                    indices.add(max(available, key=lambda i: rows[i][field]))
            for i in range(1, len(rows)):
                if rows[i].get('radiative_boundaries') != rows[i-1].get('radiative_boundaries'):
                    indices.update(range(max(0,i-2), min(len(rows),i+3)))
            keep = sorted({rows[i]['accepted'] for i in indices})
            if len(keep) > .8*len(rows):
                raise ValueError('most models carry selected structural changes')
            entries.append(dict(path=name, bytes=size, mtime_ns=mtime,
                measured_segment_cpu_seconds=report['cpu_seconds'], selected_models=keep,
                history_models=len(rows), command=command,
                executable_sha256=launch['executable_sha256'],
                retained_sha256={str(p.relative_to(ROOT)): checksum(p) for p in retained},
                scope='Replay from the retained restart; total campaign time can exceed three hours.'))
        except (OSError, ValueError, KeyError, TypeError) as error:
            excluded.append(dict(path=name, reason=str(error)))
    (HERE/'profile-retention-plan.json').write_text(json.dumps(entries, indent=2)+'\n')
    (HERE/'profile-retention-exclusions.json').write_text(json.dumps(excluded, indent=2)+'\n')
    print({'archives': len(entries), 'GiB_before_sampling': sum(e['bytes'] for e in entries)/2**30,
           'excluded': len(excluded)})


def retain(entry):
    path = ROOT/entry['path']
    destination = path.with_name('interior_selected.csv.gz')
    receipt = path.with_name('profile_retention.json')
    if not path.exists() and receipt.is_file():
        result = json.loads(receipt.read_text())
        if digest(destination) != result['selected_sha256']:
            raise ValueError(f'changed retained archive: {destination}')
        return result
    stat = path.stat()
    if (stat.st_size,stat.st_mtime_ns) != (entry['bytes'],entry['mtime_ns']):
        raise ValueError(f'changed archive: {path}')
    for name, checksum in entry['retained_sha256'].items():
        if digest(ROOT/name) != checksum:
            raise ValueError(f'changed compact research record: {name}')
    original_hash = digest(path)
    temporary = destination.with_suffix('.tmp')
    keep = set(entry['selected_models'])
    found, selected_lines, total_lines = set(), 0, 0
    with gzip.open(path, 'rb') as source, gzip.open(temporary, 'wb', compresslevel=6) as output:
        header = source.readline()
        if not header.startswith(b'years,model,cell,'):
            raise ValueError(f'unknown profile format: {path}')
        output.write(header)
        for line in source:
            if not line.strip():
                continue
            model = int(line.split(b',', 3)[1])
            total_lines += 1
            if model in keep:
                found.add(model)
                output.write(line)
                selected_lines += 1
    if found != keep:
        temporary.unlink()
        return dict(path=entry['path'], status='kept original: history/profile model mismatch',
                    missing=sorted(keep-found))
    with gzip.open(temporary, 'rb') as check:
        if sum(1 for line in check if line.strip()) != selected_lines+1:
            raise ValueError('selected archive failed integrity check')
    if temporary.stat().st_size >= stat.st_size*.8:
        temporary.unlink()
        return dict(path=entry['path'], status='kept original: little space saving')
    os.replace(temporary,destination)
    result = dict(entry, original_sha256=original_hash,
        selected_path=str(destination.relative_to(ROOT)), selected_sha256=digest(destination),
        selected_bytes=destination.stat().st_size, original_rows=total_lines,
        selected_rows=selected_lines, status='full intermediate archive removed',
        selection='Endpoints; 33 index and 33 log-age samples; Teff/L/He3/Q extrema; radiative-boundary changes +/- two models. All cells at each selected model.',
        history_and_restart_bytes_unchanged=True,
        supersedes='profile_archive.json describes the original full archive; this receipt records its intentional removal.')
    receipt.write_text(json.dumps(result,indent=2)+'\n')
    path.unlink()
    return result


def apply():
    entries = json.loads((HERE/'profile-retention-plan.json').read_text())
    last = time.monotonic()
    with (HERE/'profile-retention-receipts.jsonl').open('a',buffering=1) as log, \
            ProcessPoolExecutor(max_workers=4) as pool:
        for i,result in enumerate(pool.map(retain, entries),1):
            log.write(json.dumps(result)+'\n')
            if time.monotonic()-last>20:
                print({'processed': i, 'total': len(entries)},flush=True)
                last=time.monotonic()
    print('Profile retention complete.',flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation',choices=('plan','apply'))
    args=parser.parse_args()
    plan() if args.operation=='plan' else apply()
