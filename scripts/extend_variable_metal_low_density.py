#!/usr/bin/env python3
"""Add low-density source rows to a variable-metal FreeEOS family.

Only new densities are calculated, plus one overlap check per temperature.
Existing raw responses and all valid potential nodes must remain exact.
The original hot/dense omissions, source options and composition are retained.
A limited pilot and a full run can share the same immutable work directory.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time

from eos_source_coverage import absent_source_rows
from write_scientific_result import write_result


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_raw(path, data):
    raw = (json.dumps(data, separators=(',', ':'), allow_nan=False)+'\n').encode()
    temporary = path.with_suffix(path.suffix+'.part')
    temporary.write_bytes(gzip.compress(raw, mtime=0))
    temporary.replace(path)


def overlap_difference(a, b):
    """Compare source responses on physical scales, including zero derivatives.

    Energy-density derivatives can vanish in neutral ideal gas. Comparing
    their absolute difference to 1 erg/g would demand meaningless precision;
    use the thermal energy and entropy scales of the same source state.
    """
    energy = max(a[4]/a[2], b[4]/b[2], abs(a[10]), abs(b[10]))
    entropy = max(energy/a[3], energy/b[3], abs(a[13]), abs(b[13]))
    differences = []
    for i in range(2, len(a)):
        scale = max(abs(a[i]), abs(b[i]))
        if i in [5, 9, 10, 16, 17]:
            scale = max(scale, energy)
        elif i in [6, 11, 12, 13]:
            scale = max(scale, entropy)
        elif i not in [2, 3, 4]:
            scale = max(scale, 1.)
        differences.append(abs(a[i]-b[i])/scale)
    return max(differences)


def extend(job):
    work = Path(job['work'])
    work.mkdir(parents=True, exist_ok=True)
    receipt = work/'low_density.json'
    if receipt.exists():
        saved = json.loads(receipt.read_text())
        if saved['job'] != job:
            raise ValueError('changed extension job')
        for path, expected in saved['artifacts_sha256'].items():
            if digest(path) != expected:
                raise ValueError('completed extension changed')
        return saved
    for path, expected in job['input_sha256'].items():
        if digest(path) != expected:
            raise ValueError('source input changed: '+path)
    old = json.loads(gzip.decompress(Path(job['raw']).read_bytes()))
    absent_source_rows(old)
    if (old['options'] != [3, 223, -2]
            or old['probe_sha256'] != digest(job['probe'])
            or old.get('precision_fallback') or old.get('source_consistency_exclusions')):
        raise ValueError('unsupported source physics or precision treatment')
    qs, ts = old['logQ'], old['logT']
    step = (qs[-1]-qs[0])/(len(qs)-1)
    count = round((qs[0]-job['minimum_logQ'])/step)
    if count < 2 or abs(qs[0]-count*step-job['minimum_logQ']) > 1e-12:
        raise ValueError('new lower density must extend the original grid')
    prefix = [qs[0]-(count-i)*step for i in range(count)]
    if prefix[1]-prefix[0] != qs[1]-qs[0]:
        raise ValueError('source derivative spacing changed in floating point')
    scale = old['source_mass_scale']
    rows, overlaps, caches = [], [], {}
    started = time.monotonic()
    for it, t in enumerate(ts):
        request = ' '.join(map(str, old['eps']))+'\n3 223 -2\n'+''.join(
            f'{math.log(scale)+math.log(10)*(q+1.5*(t-6)):.17g} {math.log(10)*t:.17g}\n'
            for q in prefix+[qs[0]])
        fingerprint = hashlib.sha256((old['probe_sha256']+request).encode()).hexdigest()
        cache = work/f'temperature-{it:03d}.json.gz'
        if cache.exists():
            saved = json.loads(gzip.decompress(cache.read_bytes()))
            if saved['request'] != request or saved['input_sha256'] != fingerprint:
                raise ValueError('cached source request changed')
            extra = saved['rows']
        else:
            run = subprocess.run([job['probe']], input=request, text=True,
                                 capture_output=True, timeout=600)
            try:
                extra = [list(map(float, line.split())) for line in run.stdout.splitlines()]
            except ValueError:
                extra = []
            if (run.returncode or len(extra) != count+1
                    or any(len(r) != 22 or not all(map(math.isfinite, r)) for r in extra)):
                save_raw(work/f'failed-temperature-{it:03d}.json.gz',
                         dict(request=request, returncode=run.returncode,
                              stdout=run.stdout, stderr=run.stderr))
                raise ValueError(f'source failed at logT={t:g}')
            for r in extra:
                r[2] /= scale
                for k in [5, 6, 9, 10, 11, 12, 13, 16, 17]:
                    r[k] *= scale
            save_raw(cache, dict(request=request, input_sha256=fingerprint,
                                 rows=extra, stderr=run.stderr))
        retained = old['data'][it*len(qs):(it+1)*len(qs)]
        original, overlap = retained[0], extra[-1]
        if original is None or original[0] != 0 or overlap[0] != 0:
            raise ValueError('overlap requires two converged source states')
        difference = overlap_difference(original, overlap)
        overlaps.append(difference)
        if difference > 1e-8:
            raise ValueError(f'independent density overlap differs at logT={t:g}: {difference:g}')
        rows.extend(extra[:-1])
        rows.extend(retained)
        caches[str(cache)] = digest(cache)
    raw = {**old, 'logQ': prefix+qs, 'data': rows,
           'low_density_extension': dict(parent=job['raw'],
                                        parent_sha256=digest(job['raw']),
                                        added_source_rows=len(ts)*count,
                                        overlap_max_relative_difference=max(overlaps))}
    absent_source_rows(raw)
    for it in range(len(ts)):
        if rows[it*(len(qs)+count)+count:(it+1)*(len(qs)+count)] != old['data'][it*len(qs):(it+1)*len(qs)]:
            raise ValueError('old source rows changed')
    target = work/'source.json.gz'
    save_raw(target, raw)
    potential = Path(job['output_potential'])
    potential.parent.mkdir(parents=True, exist_ok=True)
    with (work/'import.log').open('w') as log:
        subprocess.run([sys.executable, job['importer'], str(target), str(potential)],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    before = Path(job['potential']).read_text().splitlines()
    after = potential.read_text().splitlines()
    before, after = before[before.index('data')+1:], after[after.index('data')+1:]
    old_nq, new_nq = len(qs)-4, len(qs)+count-4
    if len(before) != (len(ts)-4)*old_nq or len(after) != (len(ts)-4)*new_nq:
        raise ValueError('potential dimensions changed unexpectedly')
    retained_nodes = 0
    for it in range(len(ts)-4):
        for iq in range(old_nq):
            a, b = before[it*old_nq+iq], after[it*new_nq+count+iq]
            if a.startswith('1 '):
                if a != b:
                    raise ValueError('valid old potential node changed')
                retained_nodes += 1
            elif b != a:
                raise ValueError('old missing-node mask changed')
    for path, expected in job['input_sha256'].items():
        if digest(path) != expected:
            raise ValueError('source changed during extension')
    result = dict(job=job, raw=str(target), potential=str(potential),
                  added_source_rows=len(ts)*count, reused_source_rows=len(old['data']),
                  overlap_queries=len(ts), maximum_overlap_difference=max(overlaps),
                  preserved_valid_potential_nodes=retained_nodes,
                  failed_added_source_rows=sum(r[0] != 0 for it in range(len(ts))
                                              for r in rows[it*(len(qs)+count):it*(len(qs)+count)+count]),
                  artifacts_sha256={str(p): digest(p) for p in [target, potential]},
                  cache_sha256=caches, elapsed_seconds=time.monotonic()-started)
    write_result(receipt, result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['parent_report', 'work', 'report']:
        p.add_argument(name, type=Path)
    p.add_argument('--minimum-log-q', type=float, default=-3.)
    p.add_argument('--jobs', type=int, default=8)
    p.add_argument('--limit', type=int)
    p.add_argument('--reservation', required=True)
    a = p.parse_args()
    if a.report.exists() or not 1 <= a.jobs <= 8 or not math.isfinite(a.minimum_log_q):
        raise ValueError('invalid request or existing report')
    parent = json.loads(a.parent_report.read_text())
    if (parent['outcome'] not in ['completed', 'extended_unselected_candidate']
            or parent['failures'] or parent['partial']):
        raise ValueError('a completed source family is required')
    manifest = Path(parent['manifest']).resolve()
    if digest(manifest) != parent['manifest_sha256']:
        raise ValueError('parent manifest changed')
    a.work.mkdir(parents=True, exist_ok=True)
    source_paths = [Path(__file__), a.parent_report, manifest,
                    Path(__file__).with_name('import_freeeos_potential.py'),
                    Path(__file__).with_name('eos_source_coverage.py'),
                    Path(__file__).with_name('metal_eos_composition.py'),
                    Path(__file__).with_name('stellar_composition.py')]
    inputs = {str(q.resolve()): digest(q) for q in source_paths}
    jobs = []
    for plane in parent['planes']:
        rel = Path(plane['potential']).resolve().relative_to(manifest.parent)
        raw, potential = Path(plane['raw']), Path(plane['potential'])
        for path in [raw, potential]:
            if digest(path) != plane['artifacts_sha256'][str(path)]:
                raise ValueError('parent source plane changed')
        probe = plane['job']['probe']
        inputs[probe] = plane['job']['probe_sha256']
        jobs.append(dict(raw=str(raw), potential=str(potential),
                         work=str(a.work/rel.parent/('plane-'+rel.stem.split('-')[-1])),
                         output_potential=str(a.work/rel), probe=probe,
                         importer=str(Path(__file__).with_name('import_freeeos_potential.py').resolve()),
                         minimum_logQ=a.minimum_log_q,
                         input_sha256={str(q): digest(q) for q in [raw, potential, Path(probe)]}))
    plan = dict(input_sha256=inputs, jobs=jobs)
    path = a.work/'plan.json'
    if path.exists() and json.loads(path.read_text()) != plan:
        raise ValueError('changed plan requires another directory')
    if not path.exists():
        write_result(path, plan)
    reservation = Path(__file__).resolve().parents[1]/'docs/research/fable/coordination/reservations'/a.reservation
    if reservation.exists():
        raise FileExistsError('use a fresh reservation')
    started_utc = datetime.now(timezone.utc)
    resources = dict(owner='primary', state='running', threads=a.jobs,
                     memory_bytes=8000000000, output_bytes=2500000000,
                     controller_pid=os.getpid(), started_utc=started_utc.isoformat(),
                     scratch_paths=[str(a.work)], command=[sys.executable, *sys.argv])
    write_result(reservation, resources)
    chosen = jobs if a.limit is None else jobs[:a.limit]
    results, failures = [], []
    with ProcessPoolExecutor(a.jobs) as pool:
        futures = {pool.submit(extend, job): job for job in chosen}
        for f in as_completed(futures):
            try:
                result = f.result()
                results.append(result)
                print('completed', len(results), '/', len(chosen), result['potential'], flush=True)
            except Exception as e:
                failures.append(dict(job=futures[f], error=repr(e)))
                print('failed', failures[-1], flush=True)
    for path, expected in inputs.items():
        if digest(path) != expected:
            failures.append(dict(changed_input=path))
    complete = len(chosen) == len(jobs) and not failures
    output_manifest = a.work/'variable_metal.dat'
    if complete:
        # The relative source-plane layout is retained in the new directory.
        for line in manifest.read_text().splitlines()[4:]:
            if not (a.work/shlex.split(line)[0]).is_file():
                raise ValueError('missing output plane')
        shutil.copyfile(manifest, output_manifest)
    ended = datetime.now(timezone.utc)
    report = dict(created_utc=ended.isoformat(), outcome='completed' if not failures else 'completed_with_failures',
                  partial=not complete, selected_for_evolution=False, input_sha256=inputs,
                  plan=str(a.work/'plan.json'), planes=results, failures=failures,
                  manifest=str(output_manifest) if complete else None,
                  manifest_sha256=digest(output_manifest) if complete else None,
                  added_source_rows=sum(v['added_source_rows'] for v in results),
                  reused_source_rows=sum(v['reused_source_rows'] for v in results),
                  preserved_valid_potential_nodes=sum(v['preserved_valid_potential_nodes'] for v in results),
                  failed_added_source_rows=sum(v['failed_added_source_rows'] for v in results),
                  elapsed_utc_seconds=(ended-started_utc).total_seconds(),
                  limitations=['Source completion does not select a stellar material table; independent source/runtime and finite stellar controls are required.'])
    write_result(a.report, report)
    resources.update(state='released', ended_utc=ended.isoformat(), outcome=report['outcome'])
    reservation.write_text(json.dumps(resources, indent=2)+'\n')


if __name__ == '__main__':
    main()
