#!/usr/bin/env python3
"""Recompute EOS isotherms using checked numerical electron quadrature.

All other source rows are retained exactly. The same physical source is used
with the previously checked tighter electron quadrature. This prepares a
candidate source family; potential import and independent checks remain required.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import shutil

from metal_eos_composition import mixture


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def defects(row):
    rho, T, P = row[2:5]
    chit, Er, Et, Sr, St = row[8:13]
    return [rho*Er/P+chit-1, T*St/Et-1, rho*T*Sr/P+chit]


def repair(job):
    previous, output, record, precision, allow_exclusion, reuse, all_isotherms, reference = job
    i = record['plane']; folder = output/f'plane-{i:03d}'; folder.mkdir()
    source = previous/f'plane-{i:03d}'/'source.json.gz'
    if sha(source) != record['source_sha256']:
        raise ValueError('audited source changed')
    if reuse:
        parent, saved = reuse
        completed = parent/f'plane-{i:03d}'
        if sha(completed/'source.json.gz') != saved['source_sha256']:
            raise ValueError('completed precision repair changed')
        shutil.copytree(completed, folder, dirs_exist_ok=True)
        print(json.dumps(saved), flush=True)
        return saved
    raw = json.loads(gzip.decompress(source.read_bytes()))
    original = raw['data']
    rows = original.copy(); nq = len(raw['logQ']); m = mixture(raw['hydrogen'], raw['helium3'])
    existing = raw.get('precision_fallback_isotherms', [])
    if raw.get('precision_fallback') not in [None, precision]:
        raise ValueError('existing source precision differs')
    already_precise = {r['temperature_index'] for r in existing}
    if (len(already_precise) != len(existing)
            or any(type(r['temperature_index']) is not int or not 0 <= r['temperature_index'] < len(raw['logT'])
                   or r['logT'] != raw['logT'][r['temperature_index']] for r in existing)):
        raise ValueError('invalid existing source precision record')
    failed = {r['temperature_index'] for r in record['failures']}
    changed = sorted((set(range(len(raw['logT'])))-already_precise) | failed
                     if all_isotherms else failed)
    overrides = [r for r in existing if r['temperature_index'] not in changed]
    exclusions = []; maximum_change = 0.; maximum_defect = 0.; reference_checks = []
    for override in overrides:
        it = override['temperature_index']
        cache = json.loads(gzip.decompress((previous/f'plane-{i:03d}'/f'temperature-{it:03d}.json.gz').read_bytes()))
        fingerprint = hashlib.sha256((precision['probe_sha256']+cache['input']).encode()).hexdigest()
        if (cache.get('precision_fallback') != precision
                or cache['data'] != original[it*nq:(it+1)*nq]
                or cache.get('actual_probe_input_sha256') != fingerprint
                or override['actual_probe_input_sha256'] != fingerprint):
            raise ValueError('retained tighter isotherm lacks matching source provenance')
    for it in changed:
        cache_path = previous/f'plane-{i:03d}'/f'temperature-{it:03d}.json.gz'
        cached = json.loads(gzip.decompress(cache_path.read_bytes()))
        t = raw['logT'][it]; scale = m['source_mass_scale']
        request = ' '.join(map(str, m['eps']))+'\n'+' '.join(map(str, raw['options']))+'\n'+''.join(
            f'{math.log(scale)+math.log(10)*(q+1.5*(t-6)):.17g} {math.log(10)*t:.17g}\n' for q in raw['logQ'])
        if (cached['input'] != request or cached['data'] != original[it*nq:(it+1)*nq]
                or cached['input_sha256'] != hashlib.sha256((raw['probe_sha256']+request).encode()).hexdigest()):
            raise ValueError('cached isotherm differs from the source request or rows')
        result = subprocess.run([precision['probe']], input=request, text=True,
                                capture_output=True, timeout=600)
        failure_record = {'original_cache': str(cache_path), 'original_cache_sha256': sha(cache_path),
                          'source_validation_failures': [r for r in record['failures'] if r['temperature_index'] == it]}
        failure = folder/f'temperature-{it:03d}.nominal-identity-failure.json'
        failure.write_text(json.dumps(failure_record, indent=2)+'\n')
        saved = folder/f'temperature-{it:03d}.precision-source.json.gz'
        saved.write_bytes(gzip.compress(json.dumps(dict(input=request, stdout=result.stdout,
                         stderr=result.stderr, returncode=result.returncode)).encode(), mtime=0))
        replacement = [list(map(float, line.split())) for line in result.stdout.splitlines()]
        if (result.returncode or len(replacement) != nq
                or any(len(r) != 22 or r[0] != 0 or not all(map(math.isfinite, r)) for r in replacement)):
            raise ValueError('tighter source failed; evidence retained')
        unscaled = [r.copy() for r in replacement]
        local_change = 0.
        for j, r in enumerate(replacement):
            expected_lr = math.log(scale)+math.log(10)*(raw['logQ'][j]+1.5*(t-6))
            if abs(math.log(r[2])-expected_lr) > 1e-9 or abs(math.log(r[3])-math.log(10)*t) > 1e-10:
                raise ValueError('tighter source does not match requested coordinates')
            worst = max(map(abs, defects(r))); maximum_defect = max(maximum_defect, worst)
            if worst > 1e-7 and not allow_exclusion:
                raise ValueError('tighter source still fails the unchanged first-law criterion')
            r[2] /= scale
            for k in [5, 6, 9, 10, 11, 12, 13, 16, 17]:
                r[k] *= scale
            if max(map(abs, defects(r))) > 1e-7:
                exclusions.append(dict(index=it*nq+j, reason='first_law_defect', criterion=1e-7,
                                       maximum_absolute_defect=max(map(abs, defects(r)))))
            old = original[it*nq+j]
            local_change = max(local_change, max(abs(a-b)/max(1., abs(a), abs(b))
                                 for a, b in zip(r[2:], old[2:], strict=True)))
        maximum_change = max(maximum_change, local_change)
        if local_change > 1e-5:
            if not reference:
                raise ValueError('precision change requires an independent accuracy reference')
            check = subprocess.run([reference['probe']], input=request, text=True,
                                   capture_output=True, timeout=600)
            checked = folder/f'temperature-{it:03d}.reference-source.json.gz'
            checked.write_bytes(gzip.compress(json.dumps(dict(input=request, stdout=check.stdout,
                                stderr=check.stderr, returncode=check.returncode)).encode(), mtime=0))
            response = [list(map(float, line.split())) for line in check.stdout.splitlines()]
            if (check.returncode or len(response) != nq
                    or any(len(r) != 22 or r[0] != 0 or not all(map(math.isfinite, r)) for r in response)):
                raise ValueError('independent accuracy reference failed')
            error = max(abs(a-b)/max(1., abs(a), abs(b)) for first, second in
                        zip(unscaled, response, strict=True) for a, b in zip(first[2:], second[2:], strict=True))
            if error > 1e-8 or max(max(map(abs, defects(r))) for r in response) > 1e-7:
                raise ValueError('independent electron quadratures do not agree')
            reference_checks.append(dict(temperature_index=it, maximum_scaled_difference=error,
                                         reference_source_sha256=sha(checked)))
        rows[it*nq:(it+1)*nq] = replacement
        overrides.append(dict(temperature_index=it, logT=t,
                              actual_probe_input_sha256=hashlib.sha256((precision['probe_sha256']+request).encode()).hexdigest(),
                              nominal_failure_sha256=sha(failure),
                              reason='full electron precision calculation' if all_isotherms else 'nominal source validation failure',
                              precision_source_sha256=sha(saved)))
    if any(rows[k] != r for k, r in enumerate(original) if k//nq not in changed):
        raise ValueError('an unaffected source row changed')
    raw.update(data=rows, precision_fallback=precision, precision_fallback_isotherms=overrides,
               precision_repair_original_source_sha256=record['source_sha256'])
    if exclusions:
        raw['source_consistency_exclusions'] = exclusions
    target = folder/'source.json.gz'
    target.write_bytes(gzip.compress((json.dumps(raw, separators=(',', ':'), allow_nan=False)+'\n').encode(), mtime=0))
    (folder/'mixture.json').write_text(json.dumps(m, indent=2)+'\n')
    result = dict(plane=i, hydrogen=raw['hydrogen'], helium3=raw['helium3'],
                  recomputed_isotherms=changed, recomputed_states=len(changed)*nq,
                  original_states_retained_exactly=len(rows)-len(changed)*nq,
                  maximum_recomputed_first_law_defect=maximum_defect,
                  maximum_scaled_physical_change=maximum_change, source_sha256=sha(target))
    if reference_checks:
        result['independent_accuracy_checks'] = reference_checks
    if exclusions:
        result['inconsistent_source_exclusions'] = exclusions
    print(json.dumps(result), flush=True)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['previous', 'audit', 'precision_receipt', 'output']:
        p.add_argument(name, type=Path)
    p.add_argument('--jobs', type=int, choices=range(1, 5), default=4)
    p.add_argument('--reuse-repair', type=Path,
                   help='completed per-plane receipts from a preserved partial repair')
    p.add_argument('--exclude-inconsistent-source', action='store_true',
                   help='explicitly mask source states still failing the unchanged first-law criterion')
    p.add_argument('--all-isotherms', action='store_true',
                   help='recompute every isotherm not already using the same tighter quadrature')
    p.add_argument('--accuracy-reference', type=Path,
                   help='pinned still tighter source for independently checking changes exceeding 1e-5')
    a = p.parse_args()
    if a.output.exists():
        raise ValueError('use a fresh output directory')
    spec = json.loads((a.previous/'specification.json').read_text())
    audit = json.loads(a.audit.read_text()); precision = json.loads(a.precision_receipt.read_text())
    precision['receipt_sha256'] = sha(a.precision_receipt)
    controls = json.loads(Path(precision['source_control_report']).read_text())
    if (audit['source_specification_sha256'] != sha(a.previous/'specification.json')
            or audit['criterion'] != 1e-7 or spec['source_options'] != [3, 223, -2]
            or spec.get('precision_fallback') not in [None, precision] or spec.get('source_coverage')
            or spec.get('source_consistency_exclusion_limit')
            or precision['nominal_probe_sha256'] != spec['probe_sha256']
            or precision['source_archive_sha256'] != spec['source_archive_sha256']
            or precision['fallback_relative_integral_error'] != 1e-11
            or precision['changed_Fortran_files'] != ['src/fermi_dirac_direct.f90']
            or sha(Path(precision['probe'])) != precision['probe_sha256']
            or sha(Path(precision['library'])) != precision['library_sha256']
            or sha(Path(precision['source_control_report'])) != precision['source_control_sha256']
            or controls['queries'] < 100 or controls['maximum_scaled_physical_difference'] >= 1e-8):
        raise ValueError('incompatible source, audit or precision controls')
    mixtures = [(x, y) for x in spec['hydrogen'] for y in spec['helium3']]
    if len(audit['records']) != len(mixtures) or any(
            (r['plane'], r['hydrogen'], r['helium3']) != (i, *xy)
            for i, (r, xy) in enumerate(zip(audit['records'], mixtures, strict=True))):
        raise ValueError('audit does not cover the complete source family')
    inputs = {str(path.resolve()): sha(path) for path in
              [a.audit, a.precision_receipt, a.previous/'specification.json',
               Path(precision['probe']), Path(precision['library']), Path(__file__)]}
    reference = None
    if a.accuracy_reference:
        reference = json.loads(a.accuracy_reference.read_text())
        controls_path = Path(reference['convergence_report'])
        controls = json.loads(controls_path.read_text())
        if (reference['source_archive_sha256'] != spec['source_archive_sha256']
                or reference['previous_probe_sha256'] != precision['probe_sha256']
                or reference['relative_integral_error'] != 1e-13
                or reference['changed_Fortran_files'] != ['src/fermi_dirac_direct.f90']
                or sha(Path(reference['probe'])) != reference['probe_sha256']
                or sha(Path(reference['library'])) != reference['library_sha256']
                or sha(controls_path) != reference['convergence_report_sha256']
                or controls['probe_sha256'] != reference['probe_sha256']
                or sum(r['rows'] for r in controls['records']) < 605
                or any(not r['valid'] or r['maximum_scaled_physical_difference'] > 1e-8
                       for r in controls['records'])):
            raise ValueError('independent accuracy reference lacks matching convergence checks')
        for path in [a.accuracy_reference, controls_path, Path(reference['probe']), Path(reference['library'])]:
            inputs[str(path.resolve())] = sha(path)
    if a.all_isotherms and a.reuse_repair:
        raise ValueError('partial repair reuse is supported only for failed-isotherm repairs')
    reused = {}
    if a.reuse_repair:
        saved = json.loads(a.reuse_repair.read_text())
        parent = Path(saved['work'])
        parent_spec = json.loads((parent/'specification.json').read_text())
        if (parent_spec['precision_fallback'] != precision
                or any(parent_spec[k] != spec[k] for k in
                       ['hydrogen', 'helium3', 'logT', 'logQ', 'probe_sha256', 'source_archive_sha256'])):
            raise ValueError('partial repair has different source or precision')
        if saved['original_audit_sha256'] != sha(a.audit):
            old_audit_path = Path(saved['original_audit'])
            if sha(old_audit_path) != saved['original_audit_sha256']:
                raise ValueError('previous source audit changed')
            old_audit = json.loads(old_audit_path.read_text())
            if (old_audit['source_specification_sha256'] != audit['source_specification_sha256']
                    or [r['source_sha256'] for r in old_audit['records']] !=
                       [r['source_sha256'] for r in audit['records']]):
                raise ValueError('source rows changed between precision-repair audits')
            inputs[str(old_audit_path.resolve())] = sha(old_audit_path)
        inputs[str(a.reuse_repair.resolve())] = sha(a.reuse_repair)
        inputs[str((parent/'specification.json').resolve())] = sha(parent/'specification.json')
        for r in saved['records']:
            if r['plane'] in reused:
                raise ValueError('duplicate partial repair plane')
            # A stricter scan can identify additional failed isotherms. Reuse
            # a plane only if its completed repair covers all of those rows.
            required = {f['temperature_index'] for f in audit['records'][r['plane']]['failures']}
            if required.issubset(r['recomputed_isotherms']):
                if r.get('inconsistent_source_exclusions') and not a.exclude_inconsistent_source:
                    raise ValueError('reused plane has undeclared source consistency exclusions')
                reused[r['plane']] = (parent, r)
    a.output.mkdir()
    spec.update(precision_fallback=precision, source_identity_repair=dict(
        audit_sha256=sha(a.audit), original_source=str(a.previous.resolve()),
        first_law_criterion=1e-7, maximum_scaled_physical_change=1e-5))
    if a.all_isotherms:
        spec['source_identity_repair']['all_isotherms'] = True
    if reference:
        spec['source_identity_repair']['large_changes_require_independent_accuracy'] = dict(
            reference_receipt_sha256=sha(a.accuracy_reference), maximum_scaled_difference=1e-8)
    if a.exclude_inconsistent_source:
        spec['source_consistency_exclusion_limit'] = 1e-7
    if a.reuse_repair:
        spec['precision_repair_reuse_manifest_sha256'] = sha(a.reuse_repair)
    (a.output/'specification.json').write_text(json.dumps(spec, indent=2)+'\n')
    with ProcessPoolExecutor(a.jobs) as pool:
        records = list(pool.map(repair, [(a.previous, a.output, r, precision,
                         a.exclude_inconsistent_source, reused.get(r['plane']), a.all_isotherms,
                         reference) for r in audit['records']]))
    if any(sha(Path(name)) != value for name, value in inputs.items()):
        raise ValueError('a repair input changed')
    report = dict(scope=__doc__, input_sha256=inputs, records=records,
                  accepted_for_evolution=False,
                  total_recomputed_states=sum(r['recomputed_states'] for r in records),
                  total_original_states_retained_exactly=sum(r['original_states_retained_exactly'] for r in records))
    (a.output/'repair_manifest.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
