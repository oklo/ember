#!/usr/bin/env python3
"""Extend a FreeEOS family using only source states not already calculated.

Existing source rows and every previously valid potential node must remain
exactly unchanged. Missing cells and failed source states retain the existing
importer's masks. A partial run can be resumed with the same immutable plan.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import time
import os

from eos_source_coverage import absent_source_rows
from write_scientific_result import write_result


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def compressed_json(path, value):
    payload = (json.dumps(value, separators=(',', ':'), allow_nan=False)+'\n').encode()
    temporary = Path(str(path)+'.part')
    temporary.write_bytes(gzip.compress(payload, mtime=0))
    temporary.replace(path)


def extend_plane(job):
    started = time.monotonic()
    old_path = Path(job['parent_raw'])
    work = Path(job['work'])
    work.mkdir(parents=True, exist_ok=True)
    receipt = work/'extended.json'
    if receipt.exists():
        saved = json.loads(receipt.read_text())
        if saved['job'] != job:
            raise ValueError('changed extension job')
        for path, expected in saved['artifacts_sha256'].items():
            if digest(path) != expected:
                raise ValueError('completed extension artifact changed')
        return saved
    if (digest(old_path) != job['parent_raw_sha256']
            or digest(job['parent_potential']) != job['parent_potential_sha256']):
        raise ValueError('parent source changed')
    old = json.loads(gzip.decompress(old_path.read_bytes()))
    absent_source_rows(old)
    if (old['probe_sha256'] != job['probe_sha256'] or old['options'] != [3, 223, -2]
            or old['source_archive_sha256'] != job['source_archive_sha256']
            or old.get('precision_fallback') or old.get('source_consistency_exclusions')):
        raise ValueError('unsupported or changed parent source physics')
    ts, qs = job['logT'], job['logQ']
    if old['logQ'] != qs or ts[:len(old['logT'])] != old['logT']:
        raise ValueError('extension must preserve old coordinates exactly')
    nq = len(qs)
    data = old['data'][:] + [None]*((len(ts)-len(old['logT']))*nq)
    coverage = job['source_coverage']
    added = 0
    cache_hashes = {}
    scale = old['source_mass_scale']
    for it, t in enumerate(ts):
        indices = [it*nq+iq for iq, q in enumerate(qs)
                   if data[it*nq+iq] is None and
                   (t >= coverage['minimum_added_logT'] or q <= coverage['original_logQ_max'])]
        if not indices:
            continue
        request = ' '.join(map(str, old['eps']))+'\n3 223 -2\n'+''.join(
            f'{math.log(scale)+math.log(10)*(qs[k%nq]+1.5*(t-6)):.17g} {math.log(10)*t:.17g}\n'
            for k in indices)
        fingerprint = hashlib.sha256((job['probe_sha256']+request).encode()).hexdigest()
        cache = work/f'added-temperature-{it:03d}.json.gz'
        if cache.exists():
            saved = json.loads(gzip.decompress(cache.read_bytes()))
            if (saved['input_sha256'] != fingerprint or saved['input'] != request
                    or saved['indices'] != indices):
                raise ValueError('cached extension request differs')
            rows = saved['data']
        else:
            run = subprocess.run([job['probe']], input=request, text=True,
                                 capture_output=True, timeout=600)
            try:
                rows = [list(map(float, line.split())) for line in run.stdout.splitlines()]
            except ValueError:
                rows = []
            if (run.returncode or len(rows) != len(indices)
                    or any(len(r) != 22 or not all(map(math.isfinite, r)) for r in rows)):
                compressed_json(work/f'failed-temperature-{it:03d}.json.gz',
                                dict(input=request, returncode=run.returncode,
                                     stdout=run.stdout, stderr=run.stderr))
                raise ValueError(f'failed source request at logT={t:g}')
            # Match the original source's conversion from atomic to baryonic grams.
            for row in rows:
                row[2] /= scale
                for k in [5, 6, 9, 10, 11, 12, 13, 16, 17]:
                    row[k] *= scale
            compressed_json(cache, dict(input_sha256=fingerprint, input=request,
                                       indices=indices, data=rows, stderr=run.stderr))
        if (len(rows) != len(indices)
                or any(len(r) != 22 or not all(map(math.isfinite, r)) for r in rows)):
            raise ValueError('invalid cached extension response')
        for k, row in zip(indices, rows, strict=True):
            data[k] = row
        added += len(rows)
        cache_hashes[str(cache)] = digest(cache)
    reused = 0
    for before, after in zip(old['data'], data, strict=False):
        if before is not None:
            if before != after:
                raise ValueError('extension changed an existing source response')
            reused += 1
    raw = {**old, 'logT': ts, 'logQ': qs, 'data': data,
           'source_coverage': coverage,
           'source_extension': {'parent': str(old_path),
                                'parent_sha256': job['parent_raw_sha256'],
                                'added_source_states': added,
                                'reused_source_states': reused,
                                'added_cache_sha256': cache_hashes}}
    absent_source_rows(raw)
    target = work/'source.json.gz'
    compressed_json(target, raw)
    potential = Path(job['potential'])
    with (work/'import.log').open('w') as log:
        subprocess.run([sys.executable, job['importer'], str(target), str(potential)],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    # Appending temperature rows leaves the old potential's flattened ordering
    # intact. Every valid old stencil contains only preserved source responses.
    old_lines = Path(job['parent_potential']).read_text().splitlines()
    new_lines = potential.read_text().splitlines()
    before = old_lines[old_lines.index('data')+1:]
    after = new_lines[new_lines.index('data')+1:]
    retained_nodes = 0
    for a, b in zip(before, after[:len(before)], strict=True):
        if a.split()[0] == '1':
            if a != b:
                raise ValueError('previously valid potential node changed')
            retained_nodes += 1
    if digest(old_path) != job['parent_raw_sha256'] or digest(job['probe']) != job['probe_sha256']:
        raise ValueError('source changed while extending plane')
    report = dict(job=job, raw=str(target), potential=str(potential),
                  added_source_states=added, reused_source_states=reused,
                  preserved_valid_potential_nodes=retained_nodes,
                  absent_source_states=sum(row is None for row in data),
                  failed_source_states=sum(row is not None and row[0] != 0 for row in data),
                  masked_nodes=sum(line.split()[0] == '0' for line in after),
                  artifacts_sha256={str(target): digest(target), str(potential): digest(potential)},
                  elapsed_seconds=time.monotonic()-started)
    write_result(receipt, report)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('parent_report', type=Path)
    p.add_argument('probe', type=Path)
    p.add_argument('work', type=Path)
    p.add_argument('report', type=Path)
    p.add_argument('--maximum-logT', type=float, default=7.35)
    p.add_argument('--cool-maximum-logQ', type=float, default=1.85)
    p.add_argument('--jobs', type=int, default=8)
    p.add_argument('--only', nargs='+', help='relative raw plane directories for a partial pilot')
    p.add_argument('--reservation', type=Path, required=True)
    a = p.parse_args()
    if a.report.exists() or a.reservation.exists() or not 1 <= a.jobs <= 8:
        raise ValueError('fresh report and reservation and one to eight workers required')
    parent = json.loads(a.parent_report.read_text())
    manifest = Path(parent['manifest'])
    if digest(manifest) != parent['manifest_sha256']:
        raise ValueError('parent manifest changed')
    base = manifest.parent.resolve()
    unique = {v['raw']: v for v in parent['planes']}
    sample = json.loads(gzip.decompress(Path(next(iter(unique))).read_bytes()))
    ts, qs = sample['logT'][:], sample['logQ'][:]
    step = ts[1]-ts[0]
    count = round((a.maximum_logT-ts[-1])/step)
    if count <= 0 or not math.isclose(ts[-1]+count*step, a.maximum_logT, abs_tol=1e-12):
        raise ValueError('temperature extension must follow the original grid')
    ts += [ts[-1]+(i+1)*step for i in range(count)]
    boundary = min(qs, key=lambda q: abs(q-a.cool_maximum_logQ))
    old_coverage = sample['source_coverage']
    if (abs(boundary-a.cool_maximum_logQ) > 1e-12 or boundary <= old_coverage['original_logQ_max']
            or old_coverage.get('maximum_added_hydrogen') is not None):
        raise ValueError('unsupported density extension')
    coverage = {**old_coverage, 'original_logQ_max': boundary}
    importer = Path(__file__).with_name('import_freeeos_potential.py').resolve()
    source_paths = [Path(__file__).resolve(), importer, a.parent_report.resolve(), a.probe.resolve(),
                    Path(__file__).with_name('eos_source_coverage.py').resolve(),
                    Path(__file__).with_name('metal_eos_composition.py').resolve(),
                    Path(__file__).with_name('stellar_composition.py').resolve(),
                    Path(__file__).with_name('generate_nongrey_grid.py').resolve(),
                    Path(__file__).resolve().parents[1]/'data/atmosphere/sources/synple-elements.json',
                    Path(__file__).resolve().parents[1]/'data/opacity/sources/tops_gs98_x070_z020.request.json']
    hashes = {str(path): digest(path) for path in source_paths}
    if digest(a.probe) != sample['probe_sha256']:
        raise ValueError('extension uses a different source executable')
    a.work.mkdir(parents=True, exist_ok=True)
    plan = dict(input_sha256=hashes, logT=ts, logQ=qs, source_coverage=coverage,
                parent_manifest=str(manifest), parent_manifest_sha256=digest(manifest))
    plan_path = a.work/'extension-plan.json'
    if plan_path.exists():
        if json.loads(plan_path.read_text()) != plan:
            raise ValueError('changed extension requires a fresh work directory')
    else:
        write_result(plan_path, plan)
        archive = a.work/'scripts'
        archive.mkdir()
        for path in source_paths:
            if path.suffix == '.py':
                shutil.copy2(path, archive/path.name)
    jobs = []
    for raw_path, entry in unique.items():
        relative = Path(raw_path).resolve().relative_to(base)
        if a.only and str(relative.parent) not in a.only:
            continue
        potential_relative = Path(entry['potential']).resolve().relative_to(base)
        jobs.append(dict(parent_raw=raw_path, parent_raw_sha256=entry['raw_sha256'],
                         parent_potential=entry['potential'],
                         parent_potential_sha256=entry['potential_sha256'],
                         work=str((a.work/relative.parent).resolve()),
                         potential=str((a.work/potential_relative).resolve()),
                         probe=str(a.probe.resolve()), probe_sha256=sample['probe_sha256'],
                         source_archive_sha256=sample['source_archive_sha256'],
                         importer=str(importer), logT=ts, logQ=qs, source_coverage=coverage))
    if not jobs or (a.only and len(jobs) != len(set(a.only))):
        raise ValueError('unknown or empty plane selection')
    reservation = dict(owner='primary', state='running', controller_pid=os.getpid(),
                       threads=a.jobs, memory_bytes=6000000000, output_bytes=4000000000,
                       started_utc=datetime.now(timezone.utc).isoformat(),
                       scratch_paths=[str(a.work.resolve())], command=[sys.executable, *sys.argv],
                       task_ids=['E-METAL-EOS-COVERAGE-001'])
    write_result(a.reservation, reservation)
    started = time.monotonic()
    completed, failures = [], []
    try:
        with ProcessPoolExecutor(a.jobs) as pool:
            futures = {pool.submit(extend_plane, job): job for job in jobs}
            for future in as_completed(futures):
                try:
                    result = future.result()
                    completed.append(result)
                    print('extended', len(completed), '/', len(jobs), result['raw'], flush=True)
                except Exception as e:
                    failures.append(dict(work=futures[future]['work'], error=repr(e)))
                    print('failed', failures[-1], flush=True)
        for path, expected in hashes.items():
            if digest(path) != expected:
                failures.append(dict(changed_source=path))
        output_bytes = sum(path.stat().st_size for path in a.work.rglob('*')
                           if path.is_file() and not path.is_symlink())
        if output_bytes > reservation['output_bytes']:
            raise ValueError('extension output budget exceeded')
        manifest_result = {}
        if not failures and not a.only:
            output = a.work/'variable_metal.dat'
            # The new family uses identical relative paths and composition axes.
            text = manifest.read_text()
            for line in text.splitlines()[4:]:
                if not (a.work/json.loads(line)).is_file():
                    raise ValueError('missing extended potential')
            output.write_text(text)
            manifest_result = dict(manifest=str(output.resolve()), manifest_sha256=digest(output))
    except Exception as e:
        failures.append(dict(error=repr(e)))
        manifest_result = {}
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  outcome='extended_unselected_candidate' if not failures else 'completed_with_failures',
                  selected_for_evolution=False, partial=bool(a.only),
                  input_sha256=hashes, plan=plan, **manifest_result,
                  planes=completed, failures=failures,
                  added_source_states=sum(r['added_source_states'] for r in completed),
                  reused_source_states=sum(r['reused_source_states'] for r in completed),
                  preserved_valid_potential_nodes=sum(r['preserved_valid_potential_nodes'] for r in completed),
                  failed_source_states=sum(r['failed_source_states'] for r in completed),
                  elapsed_seconds=time.monotonic()-started,
                  limitations=['Source extension only; native profile coverage and source accuracy still require checks.',
                               'Existing fixed GS98 pattern and gas EOS approximations are unchanged.'])
    write_result(a.report, report)
    reservation.update(state='released', ended_utc=datetime.now(timezone.utc).isoformat(),
                       outcome=report['outcome'])
    a.reservation.write_text(json.dumps(reservation, indent=2)+'\n')
    print(report['outcome'], len(completed), 'planes', len(failures), 'failures', flush=True)


if __name__ == '__main__':
    main()
