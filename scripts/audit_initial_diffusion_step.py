#!/usr/bin/env python3
"""Refine only the initial conditional stellar step, reusing its coarse state."""
import argparse
import csv
import json
import os
import selectors
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from audit_full_star_diffusion import digest


YEAR = 365.25 * 86400


def time_error(coarse, fine):
    a, b = np.asarray(coarse['model_log']), np.asarray(fine['model_log'])
    delta = abs(a[:, 1:4] - b[:, 1:4])
    dx = a[:, 5:7] - b[:, 5:7]
    terms = dict(log_structure=float(delta.max() / 1e-5),
                 absolute_abundance=float(max(abs(dx).max(), abs(dx.sum(axis=1)).max()) / 1e-8),
                 relative_surface_luminosity=float(abs(a[-1, 4] / b[-1, 4] - 1) / 1e-4))
    largest = np.unravel_index(np.argmax(delta), delta.shape)
    return dict(error_terms=terms, error_norm=max(terms.values()),
                meets_driver_time_accuracy=max(terms.values()) <= 1,
                largest_structure_difference=dict(node=int(largest[0]),
                    variable=('lnr', 'lnrho', 'lnT')[largest[1]],
                    difference=float(delta[largest])),
                structure_difference_by_node=delta.tolist(),
                surface_log_difference=(a[-1, 1:4] - b[-1, 1:4]).tolist())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--scratch', type=Path, required=True)
    parser.add_argument('--levels', type=int, default=3)
    parser.add_argument('--reuse-refinement', type=Path,
                        help='Continue halving from the last successful first half in this report.')
    args = parser.parse_args()
    assert 1 <= args.levels <= 6 and not args.output.exists() and not args.scratch.exists()
    args.scratch.mkdir()
    root = Path(__file__).resolve().parents[1]
    prior = Path('/tmp/ember-full-star-timestep-v4')
    report_path = root / 'docs/results/full_star_timestep_v4.json'
    report = json.loads(report_path.read_text())
    assert report['outcome'] == 'passed_conditional_step_equations' and not report['failures']
    for name, expected in report['artifacts_sha256'].items():
        assert digest(name) == expected, name
    manifest = json.loads((prior / 'inputs.json').read_text())
    command = manifest['command']
    source_manifest = Path('/tmp/ember-full-star-timestep-source-v2/manifest.json')
    identities = {}
    for name, record in json.loads(source_manifest.read_text()).items():
        assert digest(name) == record['sha256'] == digest(record['retained']), name
        identities[name] = record['sha256']
    for name in command:
        key = str(Path(name).resolve())
        expected = manifest['input_sha256'][key if key in manifest['input_sha256'] else name]
        actual = digest(name)
        assert actual == expected, name
        identities[name] = actual
    raw_steps = [json.loads(line) for line in (prior / 'responses.jsonl').read_text().splitlines()
                 if json.loads(line).get('kind') == 'step']
    queries = (prior / 'queries.txt').read_text().splitlines(keepends=True)
    assert len(raw_steps) == len(report['steps']) == len(queries) == 14
    index, = [i for i, c in enumerate(report['steps']) if c['subdivisions'] == 8 and c['index'] == 0]
    coarse = raw_steps[index]
    assert coarse['converged'] and coarse['dt_seconds'] == 125000 * YEAR
    profile = root / 'docs/reports/2026-09-11/evolution_latest_profile.csv'
    names = ('mass_g', 'radius_cm', 'density_g_cm3', 'temperature_K', 'luminosity_erg_s', 'X', 'Y3')
    with profile.open() as stream:
        initial = np.array([[float(row[k]) for k in names] for row in csv.DictReader(stream)])
    age0 = coarse['initial_age_seconds']
    def query(state, age, logarithmic, dt):
        return f'1 1 {dt:.17g} 1e-14 1e-10 {age:.17g} {logarithmic} 512 ' + ' '.join(format(x, '.17g') for x in state.flat) + '\n'
    assert query(initial, age0, 0, coarse['dt_seconds']) == queries[index]
    if args.reuse_refinement:
        reused = json.loads(args.reuse_refinement.read_text())
        assert reused['outcome'] == 'completed_conditional_initial_refinement' and not reused['failures']
        for name, expected in reused['artifacts_sha256'].items():
            assert digest(name) == expected, name
            identities[name] = expected
        prior_inputs, = [Path(p) for p in reused['artifacts_sha256'] if Path(p).name == 'inputs.json']
        for name, expected in json.loads(prior_inputs.read_text())['input_sha256'].items():
            if digest(name) == expected:
                identities[name] = expected
            else:
                # The driver was preserved before adding this resume option.
                saved = Path('/tmp/ember-initial-refinement-source-v1.py')
                assert Path(name).resolve() == Path(__file__).resolve() and digest(saved) == expected, name
                identities[str(saved)] = expected
        raw_path, = [Path(p) for p in reused['artifacts_sha256'] if Path(p).name == 'responses.jsonl']
        query_path, = [Path(p) for p in reused['artifacts_sha256'] if Path(p).name == 'queries.txt']
        completed = [json.loads(s) for s in raw_path.read_text().splitlines()
                     if json.loads(s).get('kind') == 'step']
        recorded_queries = query_path.read_text().splitlines(keepends=True)
        assert len(completed) == len(recorded_queries) == len(reused['steps'])
        take = max(i for i, s in enumerate(reused['steps']) if s['half'] == 0)
        coarse = completed[take]
        assert coarse['converged'] and coarse['initial_age_seconds'] == age0
        assert query(initial, age0, 0, coarse['dt_seconds']) == recorded_queries[take]
        identities[str(args.reuse_refinement.resolve())] = digest(args.reuse_refinement)
    for p in (Path(__file__).resolve(), report_path, prior / 'inputs.json',
              prior / 'responses.jsonl', prior / 'queries.txt', source_manifest, profile):
        identities[str(p)] = digest(p)
    stamps = {p: (Path(p).stat().st_size, Path(p).stat().st_mtime_ns) for p in identities}
    inputs = dict(command=command, input_sha256=identities,
                  initial_age_seconds=age0,
                  retained_coarse_case=str(args.reuse_refinement) if args.reuse_refinement else [8, 0],
                  levels=args.levels, initial_interval_years=coarse['dt_seconds'] / YEAR,
                  driver_tolerances=dict(log_structure=1e-5, absolute_abundance=1e-8,
                                         relative_surface_luminosity=1e-4))
    (args.scratch / 'inputs.json').write_text(json.dumps(inputs, indent=2) + '\n')
    mass = initial[:, 0]
    weights = np.zeros(512); weights[0] = mass[0]
    weights[:-1] += .5 * np.diff(mass); weights[1:] += .5 * np.diff(mass)
    assert abs(weights.sum() / mass[-1] - 1) < 1e-14
    failures, steps, comparisons = [], [], []
    streams = [(args.scratch / name).open('x') for name in ('queries.txt', 'responses.jsonl', 'stderr.txt')]
    queries_out, raw, err = streams
    start = time.monotonic()
    proc = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            stderr=err, text=True, bufsize=1)
    process = dict(pid=proc.pid, command=command, controller_pid=os.getpid(),
                   start_utc=datetime.now(timezone.utc).isoformat(), maximum_elapsed_seconds=480)
    (args.scratch / 'native_process.json').write_text(json.dumps(process, indent=2) + '\n')
    print(json.dumps(process), flush=True)
    selector = selectors.DefaultSelector(); selector.register(proc.stdout, selectors.EVENT_READ)
    try:
        for level in range(args.levels):
            interval = coarse['dt_seconds']
            dt = interval / 2
            state, previous, age, logarithmic = initial, initial, age0, 0
            first = None
            for half in range(2):
                line = query(state, age, logarithmic, dt)
                queries_out.write(line); queries_out.flush()
                proc.stdin.write(line); proc.stdin.flush()
                record = None
                while record is None:
                    if not selector.select(120):
                        raise TimeoutError('No native record within 120 seconds; preserve partial result.')
                    reply = proc.stdout.readline()
                    if not reply: raise RuntimeError('Native process ended before its stellar result.')
                    raw.write(reply); raw.flush(); value = json.loads(reply)
                    if 'error' in value: raise RuntimeError(value['error'])
                    if value.get('kind') == 'preflight':
                        assert value['supported_nodes'] == 512 and not value['errors'], value
                    elif value.get('kind') == 'step': record = value
                    else: raise RuntimeError('Unexpected native result kind.')
                compact = {k: v for k, v in record.items() if k not in
                           ('model', 'model_log', 'total_species_rates', 'nuclear_dXdt', 'mixing_regions')}
                compact.update(level=level, half=half)
                steps.append(compact)
                assert record['converged'], record['message']
                new = np.asarray(record['model']); state = np.asarray(record['model_log'])
                assert new.shape == state.shape == (512, 7)
                assert np.isfinite(new).all() and np.isfinite(state).all()
                assert np.array_equal(new[:, 0], mass) and np.min(new[:, 1:4]) > 0
                assert np.min(new[:, 5:7]) >= 0 and np.max(new[:, 5:7].sum(axis=1)) <= .98
                rates = np.zeros((513, 2)); rates[1:-1] = record['total_species_rates']
                source = np.asarray(record['nuclear_dXdt'])
                storage = weights[:, None] * (new[:, 5:7] - previous[:, 5:7] - dt * source) / mass[-1]
                local = storage + dt / mass[-1] * np.diff(rates, axis=0)
                compact.update(independent_local_species_balance=float(abs(local).max()),
                               independent_global_species_balance=float(abs(storage.sum(axis=0)).max()),
                               mixed_regions=[r for r in record['mixing_regions'] if r[1] > r[0] + 1])
                assert record['input_preserved'] and record['initial_age_seconds'] == age
                assert record['age_seconds'] == age + dt and record['dt_seconds'] == dt
                assert abs(record['luminosity_balance']) <= 2e-8
                assert abs(record['nuclear_mass_balance']) <= 2e-6
                assert record['abundance_residual'] <= 1e-14 and record['material_heat_residual'] <= 1e-10
                assert max(compact['independent_local_species_balance'], compact['independent_global_species_balance']) <= 1e-13
                for lo, hi in compact['mixed_regions']:
                    assert np.max(abs(new[lo:hi, 5:7] - new[lo, 5:7])) == 0
                print(json.dumps(compact), flush=True)
                previous, age, logarithmic = new, record['age_seconds'], 1
                if half == 0: first = record
            comparison = dict(interval_years=interval / YEAR, **time_error(coarse, record))
            comparisons.append(comparison)
            print(json.dumps({k: v for k, v in comparison.items() if k != 'structure_difference_by_node'}), flush=True)
            if comparison['meets_driver_time_accuracy']: break
            coarse = first
        proc.stdin.close(); proc.wait(timeout=10)
        assert proc.returncode == 0
    except Exception as exc:
        failures.append(dict(exception=type(exc).__name__, message=str(exc)))
    finally:
        if proc.poll() is None:
            proc.terminate()
            try: proc.wait(timeout=3)
            except subprocess.TimeoutExpired:
                proc.kill(); proc.wait(timeout=3)
        selector.close()
        for stream in streams: stream.close()
    for p, stamp in stamps.items():
        if (Path(p).stat().st_size, Path(p).stat().st_mtime_ns) != stamp and digest(p) != identities[p]:
            failures.append(dict(changed_input=p))
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  outcome='completed_conditional_initial_refinement' if not failures else 'failed_conditional_initial_refinement',
                  accepted_for_stellar_evolution=False, selected_star_unchanged=True,
                  initial_time_accuracy_met=bool(comparisons and comparisons[-1]['meets_driver_time_accuracy']),
                  new_steps=len(steps), retained_coarse_steps=1,
                  elapsed_seconds=time.monotonic() - start, native_exit=proc.returncode,
                  failures=failures, steps=steps, comparisons=comparisons,
                  artifacts_sha256={str(p): digest(p) for p in args.scratch.iterdir()},
                  limitations=['Conditional cooler opacity and heat remain unaccepted physical assumptions.',
                               'Initial EOS/transport adjustment is included; this is not an accepted track extension.',
                               'Only the initial interval is tested; no full-track error claim follows.'])
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: result[k] for k in ('outcome', 'initial_time_accuracy_met', 'new_steps', 'failures')}), flush=True)
    if failures: raise SystemExit(1)


if __name__ == '__main__':
    main()
