#!/usr/bin/env python3
"""Divide the conditional full-star interval without recomputing its saved step.

One native process serves all dependent steps. This tests time discretization
under the stated physics; it cannot validate the conditional cool transport.
"""
import argparse
import csv
import json
import math
import os
import selectors
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from audit_full_star_diffusion import digest


def compare(coarse, fine):
    a, b = np.asarray(coarse['model']), np.asarray(fine['model'])
    scale = max(abs(a[-1, 4]), abs(b[-1, 4]))
    metrics = {
        'max_dlnr': float(np.max(abs(np.log(a[:, 1] / b[:, 1])))),
        'max_dlnrho': float(np.max(abs(np.log(a[:, 2] / b[:, 2])))),
        'max_dlnT': float(np.max(abs(np.log(a[:, 3] / b[:, 3])))),
        'max_dL_over_surface_L': float(np.max(abs(a[:, 4] - b[:, 4])) / scale),
        'max_dX': float(np.max(abs(a[:, 5:7] - b[:, 5:7]))),
        'dlnTeff': abs(math.log(coarse['Teff'] / fine['Teff'])),
        'dln_surface_L': abs(math.log(a[-1, 4] / b[-1, 4])),
        'dln_surface_R': abs(math.log(a[-1, 1] / b[-1, 1])),
    }
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('probe', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--scratch', type=Path, required=True)
    parser.add_argument('--library', type=Path, required=True)
    parser.add_argument('--reuse', type=Path, help='Retained first timestep batch; reuse only verified converged steps.')
    parser.add_argument('--reuse-report', type=Path)
    parser.add_argument('--reuse-report-sha256')
    args = parser.parse_args()
    assert not args.output.exists() and not args.scratch.exists()
    args.scratch.mkdir()
    root = Path(__file__).resolve().parents[1]
    previous = Path('/tmp/ember-full-star-diffusion-v1')
    manifest = json.loads((previous / 'inputs.json').read_text())
    report_path = root / 'docs/results/full_star_diffusion_v1.json'
    assert digest(report_path) == 'ea225e43cd018c367b1d29a6f104454521e1520c2f557c60137d202211f8e29c'
    base_report = json.loads(report_path.read_text())
    assert base_report['outcome'] == 'passed_conditional_whole_star_controls'
    retained_sources = {}
    for name in ('/tmp/ember-full-star-diffusion-source-v1/manifest.json',
                 '/tmp/ember-full-star-timestep-source-v1/manifest.json',
                 '/tmp/ember-species-merit-source-v1/manifest.json',
                 '/tmp/ember-species-merit-source-v2/manifest.json'):
        for original, saved in json.loads(Path(name).read_text()).items():
            retained_sources.setdefault(original, []).append(saved)
    identities, relocations = {}, {}
    def retain_identity(name, expected):
        if digest(name) == expected:
            identities[name] = expected
        else:
            saved, = [r for r in retained_sources[name] if r['sha256'] == expected]
            assert saved['sha256'] == expected and digest(saved['retained']) == expected, name
            identities[saved['retained']] = expected
            relocations[name] = saved['retained']
    for name, expected in manifest['input_sha256'].items():
        retain_identity(name, expected)
    for name, expected in base_report['artifacts_sha256'].items():
        retain_identity(name, expected)
    identities[str(report_path)] = digest(report_path)
    for path in (Path(__file__), root / 'scripts/full_star_diffusion_step_probe.cpp',
                 root / 'src/species_transport.cpp', args.probe, args.library):
        identities[str(path.resolve())] = digest(path)
    for name, expected in identities.items():
        assert digest(name) == expected, name
    original_replies = [json.loads(line) for line in (previous / 'responses.jsonl').read_text().splitlines()]
    baseline, = [r for r in original_replies if r.get('kind') == 'step']
    assert baseline['converged'] and baseline['dt_seconds'] == 1e6 * 365.25 * 86400
    assert baseline['screening_ions'] == 1 and baseline['factor'] == 1
    cached = {}
    if args.reuse:
        assert args.reuse_report and args.reuse_report_sha256
        reuse_report_path = args.reuse_report.resolve()
        assert digest(reuse_report_path) == args.reuse_report_sha256
        reuse_report = json.loads(reuse_report_path.read_text())
        for name, expected in reuse_report['artifacts_sha256'].items():
            retain_identity(name, expected)
        reuse_inputs = json.loads((args.reuse / 'inputs.json').read_text())
        for name, expected in reuse_inputs['input_sha256'].items():
            retain_identity(name, expected)
        identities[str(reuse_report_path)] = digest(reuse_report_path)
        old_queries = (args.reuse / 'queries.txt').read_text().splitlines(keepends=True)
        old_replies = [json.loads(line) for line in (args.reuse / 'responses.jsonl').read_text().splitlines()]
        old_steps = [r for r in old_replies if r.get('kind') == 'step']
        assert len(old_steps) == len(old_queries) == len(reuse_report['steps'])
        for record, query, checked in zip(old_steps, old_queries, reuse_report['steps']):
            if record['converged']:
                cached[checked['subdivisions'], checked['index']] = (record, query)

    profile = root / 'docs/reports/2026-09-11/evolution_latest_profile.csv'
    history = root / 'docs/reports/2026-09-11/evolution_latest.csv'
    with profile.open() as stream:
        rows = list(csv.DictReader(stream))
    names = ('mass_g', 'radius_cm', 'density_g_cm3', 'temperature_K', 'luminosity_erg_s', 'X', 'Y3')
    initial = np.array([[float(row[k]) for k in names] for row in rows])
    with history.open() as stream:
        last = list(csv.DictReader(stream))[-1]
    # The physical equations have no explicit age dependence. Their only age
    # accesses validate and increment it; the saved single step used origin 0.
    age_key = next(k for k in ('age_yr', 'age_years') if k in last)
    initial_age = float(last[age_key]) * 365.25 * 86400
    mass = initial[:, 0]
    weights = np.zeros(512)
    weights[0] = mass[0]
    weights[:-1] += .5 * np.diff(mass)
    weights[1:] += .5 * np.diff(mass)
    assert abs(weights.sum() / mass[-1] - 1) < 1e-14

    command = [str(args.probe.resolve())] + manifest['command'][1:]
    thresholds = dict(max_dlnr=1e-5, max_dlnrho=1e-5, max_dlnT=1e-5,
                      max_dL_over_surface_L=1e-5, max_dX=1e-6)
    inputs = dict(command=command, input_sha256=identities, input_relocations=relocations,
                  reused_step=baseline['dt_seconds'], reused_substeps=list(map(list, cached)),
                  subdivisions=[2, 4, 8], duration_seconds=baseline['dt_seconds'],
                  initial_age_seconds=initial_age, numerical_comparison_targets=thresholds,
                  targets_scope='Diagnostic local time accuracy, not physical-model acceptance or accumulated track error.')
    (args.scratch / 'inputs.json').write_text(json.dumps(inputs, indent=2) + '\n')
    failures, steps, endpoints = [], [], {1: baseline}
    executed_steps = 0
    started = time.monotonic()
    query_file = (args.scratch / 'queries.txt').open('x')
    raw_file = (args.scratch / 'responses.jsonl').open('x')
    err_file = (args.scratch / 'stderr.txt').open('x')
    proc = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            stderr=err_file, text=True, bufsize=1)
    (args.scratch / 'native_process.json').write_text(json.dumps(dict(
        pid=proc.pid, start_utc=datetime.now(timezone.utc).isoformat(), command=command,
        controller_pid=os.getpid(), maximum_elapsed_seconds=540), indent=2) + '\n')
    print(json.dumps(dict(native_pid=proc.pid, subdivisions=[2, 4, 8])), flush=True)
    selector = selectors.DefaultSelector()
    selector.register(proc.stdout, selectors.EVENT_READ)
    try:
        for count in (2, 4, 8):
            state, old_physical, age, logarithmic = initial, initial, initial_age, 0
            dt = baseline['dt_seconds'] / count
            for index in range(count):
                line = f'1 1 {dt:.17g} 1e-14 1e-10 {age:.17g} {logarithmic} 512 ' + ' '.join(format(x, '.17g') for x in state.flat) + '\n'
                query_file.write(line)
                query_file.flush()
                record = None
                reused = (count, index) in cached
                if reused:
                    record, original_query = cached[count, index]
                    assert line == original_query, f'Changed input for reused step {count}:{index}'
                    raw_file.write(json.dumps(record, separators=(',', ':')) + '\n')
                    raw_file.flush()
                else:
                    proc.stdin.write(line)
                    proc.stdin.flush()
                    executed_steps += 1
                while record is None:
                    if not selector.select(120):
                        raise TimeoutError('No native record within 120 seconds; preserve partial output.')
                    reply = proc.stdout.readline()
                    if not reply:
                        raise RuntimeError(f'Native process ended before step {count}:{index}')
                    raw_file.write(reply)
                    raw_file.flush()
                    value = json.loads(reply)
                    if 'error' in value:
                        raise RuntimeError(value['error'])
                    if value.get('kind') == 'preflight':
                        if value['errors'] or value['supported_nodes'] != 512:
                            raise RuntimeError(f'Preflight failed: {value}')
                    elif value.get('kind') == 'step':
                        record = value
                    else:
                        raise RuntimeError(f'Unexpected native record: {value.get("kind")}')
                label = dict(subdivisions=count, index=index)
                compact = {k: v for k, v in record.items() if k not in
                           ('model', 'model_log', 'total_species_rates', 'nuclear_dXdt', 'mixing_regions')}
                compact.update(label)
                compact['reused'] = reused
                if not record['converged']:
                    failures.append(dict(**label, message=record['message']))
                    steps.append(compact)
                    raise RuntimeError('Native step did not converge; no dependent step attempted.')
                new = np.asarray(record['model'])
                sources = np.asarray(record['nuclear_dXdt'])
                rates = np.zeros((513, 2))
                rates[1:-1] = record['total_species_rates']
                residual = (weights[:, None] * (new[:, 5:7] - old_physical[:, 5:7] - dt * sources)
                            + dt * np.diff(rates, axis=0)) / mass[-1]
                compact['independent_local_species_balance'] = float(np.max(abs(residual)))
                compact['independent_global_species_balance'] = float(np.max(abs(np.sum(
                    weights[:, None] * (new[:, 5:7] - old_physical[:, 5:7] - dt * sources), axis=0) / mass[-1])))
                compact['mixed_regions'] = [r for r in record['mixing_regions'] if r[1] > r[0] + 1]
                expected_age = age + dt
                age_error = abs(record['age_seconds'] - expected_age)
                compact['age_rounding_error_seconds'] = (record['age_seconds'] - initial_age) - (index + 1) * dt
                bad = (not record['input_preserved'] or age_error != 0 or
                       record['initial_age_seconds'] != age or record['dt_seconds'] != dt or
                       abs(record['luminosity_balance']) > 2e-8 or
                       abs(record['nuclear_mass_balance']) > 2e-6 or
                       record['abundance_residual'] > 1e-14 or record['material_heat_residual'] > 1e-10 or
                       compact['independent_local_species_balance'] > 1e-13 or
                       compact['independent_global_species_balance'] > 1e-13 or
                       new.shape != (512, 7) or not np.isfinite(new).all() or
                       not np.array_equal(new[:, 0], mass) or np.min(new[:, 1:4]) <= 0 or
                       np.min(new[:, 5:7]) < 0 or np.max(new[:, 5:7].sum(axis=1)) > .98)
                for lo, hi in compact['mixed_regions']:
                    if np.max(abs(new[lo:hi, 5:7] - new[lo, 5:7])) != 0:
                        bad = True
                if bad:
                    failures.append(dict(**label, message='Age, state, conservation or coupling check failed.'))
                steps.append(compact)
                print(json.dumps(compact), flush=True)
                if bad:
                    raise RuntimeError('Numerical control failed; no dependent step attempted.')
                state = np.asarray(record['model_log'])
                if state.shape != (512, 7) or not np.isfinite(state).all():
                    raise RuntimeError('Invalid native logarithmic checkpoint.')
                old_physical, age, logarithmic = new, record['age_seconds'], 1
            endpoints[count] = record
        proc.stdin.close()
        proc.wait(timeout=10)
        if proc.returncode:
            failures.append(dict(native_exit=proc.returncode))
    except Exception as exc:
        failures.append(dict(exception=type(exc).__name__, message=str(exc)))
    finally:
        if proc.poll() is None:
            proc.terminate()
            try:
                proc.wait(timeout=3)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait(timeout=3)
        selector.close()
        query_file.close()
        raw_file.close()
        err_file.close()

    comparisons = []
    for a, b in ((1, 2), (2, 4), (4, 8)):
        if a in endpoints and b in endpoints:
            metrics = compare(endpoints[a], endpoints[b])
            comparisons.append(dict(coarse_subdivisions=a, fine_subdivisions=b, metrics=metrics,
                                    meets_targets=all(metrics[k] <= v for k, v in thresholds.items())))
    orders = []
    for a, b in zip(comparisons, comparisons[1:]):
        orders.append(dict(subdivisions=[a['coarse_subdivisions'], a['fine_subdivisions'], b['fine_subdivisions']],
                           measured_order={k: math.log2(a['metrics'][k] / b['metrics'][k])
                                           for k in thresholds if a['metrics'][k] > 0 and b['metrics'][k] > 0}))
    for name, expected in identities.items():
        if digest(name) != expected:
            failures.append(dict(changed_input=name))
    if len(steps) != 14:
        failures.append(dict(incomplete_steps=len(steps), expected=14))
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  outcome='passed_conditional_step_equations' if not failures else 'failed_conditional_timestep_controls',
                  accepted_for_stellar_evolution=False, selected_star_unchanged=True,
                  new_steps=executed_steps, checked_substeps=len(steps), reused_steps=1 + len(cached),
                  elapsed_seconds=time.monotonic() - started,
                  native_exit=proc.returncode, failures=failures, steps=steps,
                  targets=thresholds, comparisons=comparisons, observed_orders=orders,
                  input_manifest=dict(path=str(args.scratch / 'inputs.json'), sha256=digest(args.scratch / 'inputs.json')),
                  artifacts_sha256={str(p): digest(p) for p in args.scratch.iterdir()},
                  limitations=['Conditional cooler opacity and heat treatment is unchanged and not accepted.',
                               'Local time comparison spans 1 Myr from the saved star, not accumulated full-track error.',
                               'The smooth EOS and transport differ from the saved state, so initial thermal adjustment is included.',
                               'Reported convergence orders require differences above nonlinear-solve and roundoff error.',
                               'The retained single step has age origin zero; unchanged equations have no explicit age dependence.'])
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: report[k] for k in ('outcome', 'elapsed_seconds', 'failures', 'comparisons', 'observed_orders')}, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
