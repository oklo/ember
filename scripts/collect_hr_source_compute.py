"""Freeze a conservative, explicitly approximate account of local source work.

This reads completed-job receipts; it never starts source calculations. Shared
or copied receipts with the same native input key are counted once. The output
contains the evidence needed to redraw the figure without the scratch files.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import gzip
import hashlib
import itertools
import json
import math
import os
import re

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def collect(history, output):
    rows = [json.loads(line) for line in gzip.decompress(history.read_bytes()).decode().splitlines()]
    # Single-thread native source elapsed times are a proxy for CPU time, not
    # process CPU measurements. Remove documented system sleep where possible.
    power = ROOT / 'out/runtime-sept15-v1/power-events.txt'
    sleeps = []
    for line in power.read_text().splitlines():
        m = re.match(r'(.{25}) Sleep\s+.*? (\d+) secs', line)
        if m:
            start = datetime.strptime(m[1], '%Y-%m-%d %H:%M:%S %z').timestamp()
            sleeps.append((start, start + float(m[2])))
    roots = sorted(p for p in Path('/tmp').glob('ember*') if p.is_dir() and any(
        word in p.name for word in ('nongrey', 'atmosphere', 'sources')))
    seen, jobs, omitted = {}, [], []
    constants = dict(G_cgs=6.67430e-8, mass_g=0.1 * 1.988409870698051e33,
                     sigma_cgs=5.670374419e-5, solar_luminosity=3.828e33)

    def point(t, g):
        luminosity = (4 * math.pi * constants['sigma_cgs'] * constants['G_cgs']
                      * constants['mass_g'] * t**4 / 10**g)
        return [float(t), math.log10(luminosity / constants['solar_luminosity'])]

    def specification(directory):
        for parent in (directory, *directory.parents):
            if str(parent) == '/tmp':
                break
            path = parent / 'specification.json'
            if path.exists():
                d = read_json(path)
                if d.get('teff_K') and d.get('log_g'):
                    return path, d
        return None, None

    for root in roots:
        for directory, dirs, files in os.walk(root):
            dirs[:] = sorted(d for d in dirs if d not in {
                'site-packages', '.git', '__pycache__', 'venv', 'source',
                'synple', 'synspec', 'tlusty', 'include', 'lib', 'bin'})
            if 'completed.json' not in files:
                continue
            p = Path(directory) / 'completed.json'
            try:
                receipt = read_json(p)
            except (OSError, ValueError):
                continue  # A concurrently written receipt is not a completed record.
            if not {'seconds', 'outputs', 'input_sha256'} <= receipt.keys():
                continue
            channel = ('opacity' if 'fort.63' in receipt['outputs'] else
                       'atmosphere' if 'fort.7' in receipt['outputs'] else None)
            if channel is None:
                continue
            key = (channel, receipt['input_sha256'])
            if key in seen:
                seen[key]['duplicate_receipts'].append(str(p))
                continue
            if channel == 'atmosphere':
                deck = p.parent / 'fort.5'
                if deck.exists():
                    first = deck.open().readline()
                elif deck.with_suffix('.5.gz').exists():
                    first = gzip.decompress(deck.with_suffix('.5.gz').read_bytes()).decode().splitlines()[0]
                else:
                    omitted.append(dict(path=str(p), reason='No retained coordinate input'))
                    continue
                try:
                    t, g = map(float, first.split()[:2])
                    assert 1000 <= t <= 10000 and 3 <= g <= 7
                except (ValueError, AssertionError):
                    omitted.append(dict(path=str(p), reason='Not a stellar atmosphere coordinate'))
                    continue
                points, evidence = [point(t, g)], dict(teff_K=t, log_g=g, input_first_line=first.strip())
            else:
                spec_path, spec = specification(p.parent)
                if spec is None:
                    omitted.append(dict(path=str(p), reason='No retained target atmosphere grid'))
                    continue
                coords = sorted(set(itertools.product(spec['teff_K'], spec['log_g'])))
                points = [point(t, g) for t, g in coords]
                evidence = dict(specification=str(spec_path), specification_sha256=digest(spec_path),
                                target_teff_logg=coords)
            seconds = float(receipt['seconds'])
            end = p.stat().st_mtime
            removed = sum(max(0, min(end, b) - max(end-seconds, a)) for a, b in sleeps)
            seconds = max(0, seconds-removed)
            job = dict(channel=channel, receipt=str(p), receipt_sha256=digest(p),
                       native_input_sha256=receipt['input_sha256'], recorded_seconds=receipt['seconds'],
                       documented_sleep_removed_seconds=removed, estimated_worker_seconds=seconds,
                       points=points, allocation='equal shares over listed HR coordinates',
                       coordinate_evidence=evidence, duplicate_receipts=[])
            seen[key] = job
            jobs.append(job)

    # Source-plane sums are preferable to parallel batch wall time. Where only
    # batch timing survives, worker count times active elapsed time is an upper
    # estimate, including importer overhead and unused tail capacity.
    eos = [
        ('variable_metal_eos_global_grid_v1', 1, 8, 'batch'),
        ('variable_metal_eos_extension_v1', 1, None, 'planes'),
        ('variable_metal_low_density_v2', 4, None, 'planes'),
        ('high_metal_eos_source_v1', 8, 3, 'batch'),
    ]
    for name, segment, workers, mode in eos:
        source = ROOT / f'docs/results/{name}.json'
        d = read_json(source)
        run_path = ROOT / f'docs/results/metal_main_sequence_v{segment}.json'
        run = read_json(run_path)
        age = run['initial_age_seconds'] / 31557600
        r = min(rows, key=lambda r: abs(r['age_years']-age))
        assert abs(r['age_years']-age) < 0.01
        times = [p['elapsed_seconds'] for p in d['planes']] if mode == 'planes' else []
        seconds = sum(times) if times else workers*d['elapsed_seconds']
        jobs.append(dict(channel='eos', receipt=str(source), receipt_sha256=digest(source),
                         estimated_worker_seconds=seconds,
                         timing_method='sum of serial plane worker times' if times else
                         'active batch elapsed time times worker count; upper estimate',
                         worker_count=workers, recorded_plane_seconds=times,
                         recorded_batch_seconds=d.get('elapsed_seconds'),
                         points=[[r['Teff'], math.log10(r['luminosity']/constants['solar_luminosity'])]],
                         allocation='first accepted-track state using this EOS family',
                         first_use_age_years=r['age_years'], first_use_report=str(run_path),
                         first_use_report_sha256=digest(run_path)))

    assert all(math.isfinite(j['estimated_worker_seconds']) and j['estimated_worker_seconds'] >= 0 for j in jobs)
    data = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                scope='Retained local spectral-opacity and atmosphere source jobs, plus four selected variable-metal EOS construction batches; not total project compute.',
                history_sha256=digest(history), collector_sha256=digest(__file__), constants=constants,
                channels=['opacity', 'eos', 'atmosphere'], jobs=jobs, omitted_receipts=omitted,
                channel_worker_hours={c: sum(j['estimated_worker_seconds'] for j in jobs if j['channel']==c)/3600
                                      for c in ('opacity', 'eos', 'atmosphere')},
                documented_sleep_log=dict(path=str(power), sha256=digest(power), intervals_unix=sleeps),
                duplicate_receipts_excluded=sum(len(j.get('duplicate_receipts', [])) for j in jobs),
                limitations=[
                    'Serial native wall durations approximate worker CPU time; remaining scheduling, sleep and I/O effects are not measured. Documented sleep is subtracted using receipt modification times as finish-time estimates.',
                    'EOS batch-worker estimates include idle tail capacity; they are not process CPU measurements.',
                    'Shared spectral construction is divided once over its declared target atmosphere grid, not charged again at every stellar timestep.',
                    'Atmosphere coordinates are projected into the HR diagram at fixed mass using L=4 pi sigma G M Teff^4/g. They describe source-grid coverage, not additional evolved stellar models.',
                    'EOS construction is allocated at first track use, not where interior material temperatures lie on the HR axes.',
                    'Different compositions and completed numerical controls are summed. Identical native input keys are counted once, conservatively excluding both copies and any exact repeated runs.',
                    'Unfinished jobs, missing receipts, remote-service compute, much older EOS work, evolution solves, compilation and human or agent work are excluded.',
                    'Completed native stages can include atmosphere attempts subsequently rejected by independent validation; color is compute expenditure, not a validity mask.',
                    'Source work for comparison trajectories and advance coverage is included; some cells lie beyond the current primary track.',
                ])
    output.write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps(dict(jobs=len(jobs), omitted=len(omitted), hours=data['channel_worker_hours'],
                          duplicates=data['duplicate_receipts_excluded']), indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('history', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    collect(a.history, a.output)
