"""Account for native CPU work in the model sequences retained in a plot.

Exact saved age, temperature and luminosity identify a plotted state. Only
runs contributing those states are included. Required timestep comparisons
and adaptive retries within those runs remain part of their computation;
independent controls and superseded tracks are excluded.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import collections
import csv
import gzip
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads(Path(path).read_text())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def row_key(row):
    return row['age_years'], row['Teff'], row['luminosity']


def collect(directory, output):
    tracks = {name: [json.loads(line) for line in gzip.decompress(
        (directory / filename).read_bytes()).decode().splitlines()]
        for name, filename in [('metal', 'metal_history.jsonl.gz'),
                               ('fixed', 'consistent_cn_history.jsonl.gz')]}
    track_keys = {name: set(map(row_key, rows)) for name, rows in tracks.items()}
    wanted = set.union(*track_keys.values())
    contributors = collections.defaultdict(list)
    candidates = []
    for prefix in ['cn_main_sequence', 'metal_main_sequence']:
        candidates += [p for p in (ROOT/'docs/results').glob(prefix+'_v*.json')
                       if p.stem.rsplit('_v', 1)[-1].isdigit()]
    candidates += list((ROOT/'out/completed-metal-segments').glob('*/result.json'))
    candidates += [ROOT/'docs/results/cn_x045_crossing_v1.json']
    runs = []
    for path in sorted(candidates):
        data = read(path)
        rows = data.get('history', [])
        used = [r for r in rows if row_key(r) in wanted]
        if not used:
            continue
        if len(used) != len(rows):
            raise ValueError(f'Partly retained run needs interval accounting: {path}')
        cpu = data.get('native_total_cpu_seconds')
        if cpu is None or cpu < 0:
            raise ValueError(f'Missing native process CPU time: {path}')
        memberships = {name: sum(row_key(r) in keys for r in used)
                       for name, keys in track_keys.items()}
        # Each current segment belongs wholly to the common history or one
        # branch. Do not apportion a shared process time by arbitrary weights.
        if all(n == len(rows) for n in memberships.values()):
            group = 'shared_CN_prefix'
        elif memberships['metal'] == len(rows) and memberships['fixed'] == 0:
            group = 'metal_branch'
        elif memberships['fixed'] == len(rows) and memberships['metal'] == 0:
            group = 'fixed_metal_branch'
        else:
            raise ValueError(f'Mixed segment needs separate accounting: {path}')
        run = dict(report=str(path), report_sha256=digest(path),
                   states=len(rows), first_age_years=rows[0]['age_years'],
                   last_age_years=rows[-1]['age_years'], group=group,
                   native_cpu_seconds=cpu, cpu_field='native_total_cpu_seconds',
                   track_memberships=memberships)
        runs.append(run)
        for row in used:
            contributors[row_key(row)].append(str(path))
    missing = wanted - contributors.keys()
    duplicated = {str(k): paths for k, paths in contributors.items() if len(paths) != 1}
    if missing or duplicated:
        raise ValueError(f'Incomplete or duplicate provenance: {len(missing)} missing, {duplicated}')

    pp_provenance = ROOT/'docs/reports/2026-09-11/evolution_latest_provenance.json'
    pp = read(pp_provenance)
    if digest(directory/'pp_comparison_history.csv') != pp['history_csv_sha256']:
        raise ValueError('PP plotted history differs from its timing provenance')
    receipts = []
    for item in pp['segment_receipts']:
        path = Path(item['source'])
        if digest(path) != item['sha256']:
            raise ValueError(f'Changed PP receipt: {path}')
        data = read(path)
        receipts.append(dict(path=str(path), sha256=digest(path),
                             cpu_seconds=data['child_user_seconds']+data['child_system_seconds']))
    pp_cpu = sum(r['cpu_seconds'] for r in receipts)
    if abs(pp_cpu-pp['cpu_seconds']) > 1e-8:
        raise ValueError('PP summed CPU differs from frozen account')
    grouped = collections.defaultdict(float)
    for run in runs:
        grouped[run['group']] += run['native_cpu_seconds']
    grouped['pp_comparison'] = pp_cpu
    total = sum(grouped.values())
    rates = dict(central=2e9, low=.6e9, high=6e9)
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  scope='Native CPU work in only the retained plotted sequences, with shared history counted once',
                  sequence_states={name: len(rows) for name, rows in tracks.items()},
                  shared_states=len(set.intersection(*track_keys.values())),
                  unique_CN_and_metal_states=len(wanted),
                  all_plotted_CN_and_metal_states_matched_once=True,
                  runs=runs, pp_receipts=receipts,
                  cpu_seconds_by_group=dict(grouped),
                  cpu_hours_by_group={k: v/3600 for k, v in grouped.items()},
                  total_cpu_seconds=total, total_cpu_hours=total/3600,
                  estimated_flops_per_cpu_second=rates,
                  estimated_operations={k: total*v for k, v in rates.items()},
                  standalone_track_cpu_hours=dict(
                      metal=(grouped['shared_CN_prefix']+grouped['metal_branch'])/3600,
                      fixed=(grouped['shared_CN_prefix']+grouped['fixed_metal_branch'])/3600,
                      pp=pp_cpu/3600),
                  input_sha256={str(p): digest(p) for p in [
                      Path(__file__), pp_provenance, directory/'metal_history.jsonl.gz',
                      directory/'consistent_cn_history.jsonl.gz', directory/'pp_comparison_history.csv']},
                  limitations=[
                      'Counts native process user plus system CPU, including table loading, required full/two-half-step checks and adaptive retries in the retained runs. Python orchestration and compilation are excluded.',
                      'Only processes supplying the saved curves are counted. Independent physical controls, obsolete sequences and superseded portions of the track are excluded.',
                      'The 2 GFLOP/s central conversion, range 0.6–6, is a broad engineering estimate informed by CPU arithmetic-kernel measurements; it is not a measurement of whole-program Ember throughput. CPU time is the stronger comparison.',
                      'The LBA97 curve is read from the publication and has no local sequence-computation charge. Source-building time is reported separately.',
                      'The metal and fixed-metal standalone totals overlap; add the disjoint groups, not the standalone totals.'])
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(hours=result['cpu_hours_by_group'], total_hours=total/3600,
                          operations=result['estimated_operations'], shared_states=result['shared_states']), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    collect(args.directory, args.output)
