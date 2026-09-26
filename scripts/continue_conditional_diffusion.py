#!/usr/bin/env python3
"""Measure adaptive full-star diffusion after a verified conditional first step.

This does not select the conditional opacity or cool-envelope heat assumptions.
One native process retains the physics tables; rejected trials are preserved.
"""
import argparse
import json
import math
import os
import selectors
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from audit_initial_diffusion_step import YEAR, digest, time_error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--scratch', type=Path, required=True)
    parser.add_argument('--duration-years', type=float, default=1e6)
    parser.add_argument('--maximum-attempts', type=int, default=30)
    parser.add_argument('--resume-checkpoint', type=Path)
    parser.add_argument('--resume-report', type=Path)
    parser.add_argument('--resume-audit', type=Path, help='Independent audit required to resume a completed prefix of an interrupted run.')
    parser.add_argument('--maximum-output-bytes', type=int, default=35000000)
    parser.add_argument('--abundance-relative-tolerance', type=float, default=0.)
    parser.add_argument('--abundance-absolute-tolerance', type=float, default=1e-8)
    parser.add_argument('--structure-tolerance', type=float, default=1e-5)
    parser.add_argument('--surface-luminosity-tolerance', type=float, default=1e-4)
    parser.add_argument('--coupling-abundance-tolerance', type=float, default=1e-14)
    parser.add_argument('--coupling-heat-tolerance', type=float, default=1e-10)
    parser.add_argument('--probe', type=Path, help='Explicit alternative native executable for a numerical comparison.')
    parser.add_argument('--fuel-relative-tolerance', type=float, default=1e-6)
    args = parser.parse_args()
    assert not args.output.exists() and not args.scratch.exists()
    assert 15625 < args.duration_years <= 1e9 and 1 <= args.maximum_attempts <= 100
    assert 1000000 <= args.maximum_output_bytes <= 100000000
    assert bool(args.resume_checkpoint) == bool(args.resume_report)
    assert 0 <= args.abundance_relative_tolerance <= .01
    assert all(math.isfinite(x) and 0 < x <= .01 for x in
               (args.abundance_absolute_tolerance,args.structure_tolerance,args.surface_luminosity_tolerance))
    assert 1e-15 <= args.coupling_abundance_tolerance <= 1e-10
    assert 1e-12 <= args.coupling_heat_tolerance <= 1e-5
    assert 0 < args.fuel_relative_tolerance <= 1e-3
    args.scratch.mkdir()
    root = Path(__file__).resolve().parents[1]
    initial_report = args.resume_report.resolve() if args.resume_report else root / 'docs/results/full_star_initial_refinement_v2.json'
    r = json.loads(initial_report.read_text())
    if r['failures']:
        assert args.resume_audit and args.resume_checkpoint
        audit=json.loads(args.resume_audit.read_text())
        assert audit['outcome']=='passed' and audit['scope']=='completed_prefix'
        assert Path(audit['source_run_report']).resolve()==initial_report
        assert Path(audit['resumable_checkpoint']).resolve()==args.resume_checkpoint.resolve()
        pinned={Path(p).resolve():h for p,h in audit['input_sha256'].items()}
        assert pinned[initial_report]==digest(initial_report)
        assert pinned[args.resume_checkpoint.resolve()]==digest(args.resume_checkpoint)
        assert audit['total_conditional_elapsed_years']==r['elapsed_years']
    if not args.resume_checkpoint:assert r['initial_time_accuracy_met']
    identities = {str(initial_report): digest(initial_report)}
    if args.resume_audit:identities[str(args.resume_audit.resolve())]=digest(args.resume_audit)
    for name, expected in r['artifacts_sha256'].items():
        assert digest(name) == expected, name
        identities[name] = expected
    raw_path, = [Path(p) for p in identities if Path(p).name == 'responses.jsonl']
    input_path, = [Path(p) for p in identities if Path(p).name == 'inputs.json']
    prior = json.loads(input_path.read_text())
    for name, expected in prior['input_sha256'].items():
        actual = Path(name)
        if actual == Path(__file__).resolve() and digest(actual) != expected:
            actual, = [p for p in (Path('/tmp/ember-adaptive-continuation-source-v1.py'),
                                   Path('/tmp/ember-adaptive-continuation-source-v3.py'),
                                   Path('/tmp/ember-adaptive-continuation-source-v4.py'),
                                   Path('/tmp/ember-adaptive-continuation-source-v5.py'),
                                   Path('/tmp/ember-adaptive-continuation-source-v6.py'),
                                   Path('/tmp/ember-adaptive-continuation-source-v7.py'),
                                   Path('/tmp/ember-adaptive-continuation-source-v8.py')) if digest(p)==expected]
        if digest(actual) != expected:
            snapshot = json.loads(Path('/tmp/ember-parcel-radiation-source-v1/manifest.json').read_text())
            saved = snapshot[str(actual)]
            assert saved['sha256'] == expected
            actual = Path(saved['retained'])
        assert digest(actual) == expected, name
        identities[str(actual)] = expected
    identities[str(Path(__file__).resolve())] = digest(Path(__file__))
    if args.resume_checkpoint:
        cp = args.resume_checkpoint.resolve()
        assert cp in {Path(p).resolve() for p in r['artifacts_sha256']}
        checkpoint = json.loads(cp.read_text())
        current = checkpoint if r.get('experiment')=='cn' else checkpoint['model_record']
        assert current['converged']
        if r.get('experiment')=='cn':
            assert r['outcome']=='passed' and current['cn_carbon']==1
            assert cp.name=='CN-endpoint.json' and prior['command'][-1]=='select-cn'
            parent=Path(prior['initial_checkpoint'])
            assert str(parent) in prior['input_sha256']
            parent_state=json.loads(parent.read_text())
            initial_elapsed=parent_state['elapsed_seconds']+r['additional_interval_years']*YEAR
            assert initial_elapsed/YEAR==r['initial_conditional_elapsed_years']+r['additional_interval_years']
            initial_interval=r['additional_interval_years']*YEAR
            initial_error=max(r['comparisons']['CN']['time_error_terms'].values())
            age0=parent_state['model_record']['age_seconds']-parent_state['elapsed_seconds']
            command=prior['command'][:-1]+['cn']
        elif 'elapsed_seconds' in checkpoint:
            initial_elapsed = checkpoint['elapsed_seconds']
            initial_interval = checkpoint['numerical_history_row']['interval_years'] * YEAR
            initial_error = checkpoint['numerical_history_row']['error_norm']
            command = prior['command']
            age0 = prior['initial_age_seconds']
        else:
            assert r['outcome'] == 'completed_conditional_parcel_radiation_comparison'
            parent = Path(checkpoint['initial_checkpoint'])
            assert str(parent) in prior['input_sha256']
            initial_elapsed = json.loads(parent.read_text())['elapsed_seconds'] + checkpoint['additional_elapsed_seconds']
            initial_interval = r['additional_interval_years'] * YEAR
            initial_error = r['time_accuracy'][str(current['parcel_radiation'])]['error_norm']
            assert r['time_accuracy'][str(current['parcel_radiation'])]['meets_driver_time_accuracy']
            assert prior['command'][-1] == 'select-radiation'
            command = prior['command'][:-1] + (['radiation'] if current['parcel_radiation'] else [])
            parent_inputs = parent.parent / 'inputs.json'
            assert str(parent_inputs) in prior['input_sha256']
            age0 = json.loads(parent_inputs.read_text())['initial_age_seconds']
    else:
        raw_initial = [json.loads(s) for s in raw_path.read_text().splitlines()
                       if json.loads(s).get('kind') == 'step']
        assert len(raw_initial) == 2 and all(s['converged'] for s in raw_initial)
        current = raw_initial[-1]
        initial_elapsed = initial_interval = r['comparisons'][-1]['interval_years'] * YEAR
        initial_error = r['comparisons'][-1]['error_norm']
        age0 = prior['initial_age_seconds']
        command = prior['command']
    if args.probe:
        command = [str(args.probe.resolve())] + command[1:]
        identities[command[0]] = digest(command[0])
    assert initial_elapsed < args.duration_years * YEAR and 0 < initial_error <= 1
    elapsed = initial_elapsed
    dt = initial_interval * float(np.clip(.9 / math.sqrt(initial_error), .5, 2))
    mass = np.asarray(current['model'])[:, 0]
    weights = np.zeros(512); weights[0] = mass[0]
    weights[:-1] += .5 * np.diff(mass); weights[1:] += .5 * np.diff(mass)
    assert abs(weights.sum() / mass[-1] - 1) < 1e-14
    inputs = dict(command=command, input_sha256=identities,
                  initial_age_seconds=age0, reused_initial_interval_seconds=initial_elapsed,
                  target_elapsed_seconds=args.duration_years * YEAR,
                  resumed_checkpoint=str(args.resume_checkpoint) if args.resume_checkpoint else None,
                  maximum_attempts=args.maximum_attempts, maximum_output_bytes=args.maximum_output_bytes,
                  coupling_tolerances=dict(abundance=args.coupling_abundance_tolerance,heat=args.coupling_heat_tolerance),
                  numerical_tolerances=dict(log_structure=args.structure_tolerance,
                     absolute_abundance=args.abundance_absolute_tolerance,relative_abundance=args.abundance_relative_tolerance,
                     relative_surface_luminosity=args.surface_luminosity_tolerance,
                     relative_integrated_absolute_hydrogen=args.fuel_relative_tolerance,
                     convective_boundary_cells=1))
    (args.scratch / 'inputs.json').write_text(json.dumps(inputs, indent=2) + '\n')
    stamps = {p: (Path(p).stat().st_size, Path(p).stat().st_mtime_ns) for p in identities}
    streams = [(args.scratch / p).open('x') for p in ('queries.txt', 'responses.jsonl', 'stderr.txt', 'steps.jsonl', 'history.jsonl')]
    queries_out, raw, err, journal, history = streams
    start = time.monotonic()
    proc = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            stderr=err, text=True, bufsize=1)
    process = dict(pid=proc.pid, command=command, controller_pid=os.getpid(),
                   start_utc=datetime.now(timezone.utc).isoformat(), maximum_elapsed_seconds=840)
    (args.scratch / 'native_process.json').write_text(json.dumps(process, indent=2) + '\n')
    print(json.dumps(process), flush=True)
    selector = selectors.DefaultSelector(); selector.register(proc.stdout, selectors.EVENT_READ)
    cache, attempted, accepted, failures = {}, [], [], []
    counts = dict(native_steps=0, reused_trials=0, numerical_rejections=0, equation_rejections=0,
                  verification_retries=0,verification_rejections=0)

    def assess_time_error(coarse,fine):
        result=time_error(coarse,fine)
        a,b=np.asarray(coarse['model_log']),np.asarray(fine['model_log'])
        xa=np.column_stack((a[:,5:7],.98-a[:,5:7].sum(axis=1)))
        xb=np.column_stack((b[:,5:7],.98-b[:,5:7].sum(axis=1)))
        independent=a[:,5:7]-b[:,5:7]
        delta=np.column_stack((independent,-independent.sum(axis=1)))
        scale=args.abundance_absolute_tolerance+args.abundance_relative_tolerance*np.maximum(abs(xa),abs(xb))
        error=abs(delta)/scale;location=np.unravel_index(np.argmax(error),error.shape)
        terms=dict(log_structure=float(abs(a[:,1:4]-b[:,1:4]).max()/args.structure_tolerance),
                   scaled_abundance=float(error.max()),
                   relative_surface_luminosity=float(abs(a[-1,4]/b[-1,4]-1)/args.surface_luminosity_tolerance),
                   integrated_hydrogen=float(np.dot(weights,abs(independent[:,0]))/
                       max(np.dot(weights,b[:,5]),1e-300)/args.fuel_relative_tolerance))
        regions = [np.asarray([r for r in s['mixing_regions'] if r[1]>r[0]+1],dtype=int).reshape(-1,2)
                   for s in (coarse,fine)]
        terms['convective_boundary_cells'] = (float(abs(regions[0]-regions[1]).max()) if regions[0].size else 0.) \
            if regions[0].shape == regions[1].shape else 2.
        result.update(error_terms=terms,error_norm=max(terms.values()),meets_driver_time_accuracy=max(terms.values())<=1,
            largest_abundance_error=dict(node=int(location[0]),species=('H','He3','He4')[location[1]],
                absolute_error=float(abs(delta)[location]),normalizing_scale=float(scale[location])))
        return result

    def evolve(old, span, attempt, part, abundance_tolerance=None, verification_retry=0):
        if abundance_tolerance is None:abundance_tolerance=args.coupling_abundance_tolerance
        state = np.asarray(old['model_log'])
        age = old['age_seconds']
        line = f'1 1 {span:.17g} {abundance_tolerance:.17g} {args.coupling_heat_tolerance:.17g} {age:.17g} 1 512 ' + ' '.join(format(x, '.17g') for x in state.flat) + '\n'
        reuse = line in cache
        if reuse:
            record = cache[line]; counts['reused_trials'] += 1
        else:
            queries_out.write(line); queries_out.flush()
            proc.stdin.write(line); proc.stdin.flush()
            record = None
            while record is None:
                if not selector.select(120): raise TimeoutError('No native result in 120 seconds.')
                reply = proc.stdout.readline()
                if not reply: raise RuntimeError('Native process ended before the stellar result.')
                raw.write(reply); raw.flush(); value = json.loads(reply)
                if 'error' in value: raise RuntimeError(value['error'])
                if value.get('kind') == 'preflight':
                    if value['supported_nodes'] != 512 or value['errors']:
                        raise RuntimeError('Unsupported initial state: ' + json.dumps(value))
                elif value.get('kind') == 'step': record = value
                else: raise RuntimeError('Unexpected native result.')
            counts['native_steps'] += 1
            if record['converged']: cache[line] = record
        compact = {k: v for k, v in record.items() if k not in
                   ('model', 'model_log', 'total_species_rates', 'nuclear_dXdt', 'mixing_regions')}
        compact.update(attempt=attempt, part=part, reused=reuse,verification_retry=verification_retry)
        verification_failures=[]
        if record['converged']:
            new = np.asarray(record['model']); log = np.asarray(record['model_log'])
            assert new.shape == log.shape == (512, 7) and np.isfinite(new).all() and np.isfinite(log).all()
            assert np.array_equal(new[:, 0], mass) and np.min(new[:, 1:4]) > 0
            assert np.min(new[:, 5:7]) >= 0 and np.max(new[:, 5:7].sum(axis=1)) <= .98
            assert record['input_preserved'] and record['initial_age_seconds'] == age
            if 'parcel_radiation' in current:assert record['parcel_radiation'] == current['parcel_radiation']
            if 'cn_carbon' in current:assert record['cn_carbon']==current['cn_carbon']
            assert record['age_seconds'] == age + span and record['dt_seconds'] == span
            rates = np.zeros((513, 2)); rates[1:-1] = record['total_species_rates']
            source = np.asarray(record['nuclear_dXdt'])
            storage = weights[:, None] * (new[:, 5:7] - np.asarray(old['model'])[:, 5:7] - span * source) / mass[-1]
            local = storage + span / mass[-1] * np.diff(rates, axis=0)
            compact.update(independent_local_species_balance=float(abs(local).max()),
                           independent_global_species_balance=float(abs(storage.sum(axis=0)).max()),
                           mixed_regions=[r for r in record['mixing_regions'] if r[1] > r[0] + 1])
            criteria=[('local species',compact['independent_local_species_balance'],1e-13),
                      ('global species',compact['independent_global_species_balance'],1e-13),
                      ('energy',abs(record['luminosity_balance']),2e-8),
                      ('nuclear mass',abs(record['nuclear_mass_balance']),2e-6),
                      ('abundance iteration',record['abundance_residual'],abundance_tolerance),
                      ('heat iteration',record['material_heat_residual'],args.coupling_heat_tolerance)]
            verification_failures=[dict(check=name,value=value,limit=limit) for name,value,limit in criteria if not math.isfinite(value) or value>limit]
            for lo, hi in compact['mixed_regions']:
                assert np.max(abs(new[lo:hi, 5:7] - new[lo, 5:7])) == 0
        compact.update(verification_passed=record['converged'] and not verification_failures,
                       verification_failures=verification_failures)
        journal.write(json.dumps(compact, separators=(',', ':')) + '\n'); journal.flush()
        if sum(s.tell() for s in (queries_out, raw, journal, history)) > args.maximum_output_bytes:
            raise RuntimeError('Scientific output limit reached; preserve the current checkpoint.')
        if verification_failures:
            cache.pop(line,None)
            counts['verification_rejections']+=1
            print(json.dumps(dict(attempt=attempt,part=part,verification_failures=verification_failures,
                                  abundance_tolerance=abundance_tolerance)),flush=True)
            if verification_retry<2 and abundance_tolerance>1e-12:
                counts['verification_retries']+=1
                return evolve(old,span,attempt,part,max(1e-12,.1*abundance_tolerance),verification_retry+1)
            # Keep the native reply/journal intact. Only the step controller
            # sees a rejected result and reduces its trial interval.
            return dict(record,converged=False,message='Independent checks failed after tighter abundance solves.')
        return record

    def save_state(record, interval, error):
        row = dict(elapsed_years=elapsed / YEAR, age_seconds=record['age_seconds'],
                   age_rounding_error_seconds=record['age_seconds'] - age0 - elapsed,
                   interval_years=interval / YEAR, error_norm=error,
                   Teff=record['Teff'], central_hydrogen=record['model'][0][5],
                   surface_hydrogen=record['model'][-1][5],
                   convective_mass_fraction=record['convective_mass_fraction'])
        history.write(json.dumps(row) + '\n'); history.flush()
        checkpoint = dict(accepted_for_stellar_evolution=False,
                          elapsed_seconds=elapsed, model_record=record, numerical_history_row=row)
        tmp = args.scratch / 'checkpoint.json.tmp'; tmp.write_text(json.dumps(checkpoint) + '\n')
        os.replace(tmp, args.scratch / 'checkpoint.json')
        return row

    save_state(current, initial_interval, initial_error)
    try:
        for index in range(args.maximum_attempts):
            if time.monotonic() - start > 840:
                raise TimeoutError('Bounded continuation time limit reached.')
            if elapsed >= args.duration_years * YEAR: break
            dt = min(dt, args.duration_years * YEAR - elapsed)
            if dt < YEAR: raise RuntimeError('Minimum timestep reached.')
            entry = dict(attempt=index, initial_elapsed_years=elapsed / YEAR, interval_years=dt / YEAR)
            attempted.append(entry)
            full = evolve(current, dt, index, 'full')
            if full['converged']:
                first = evolve(current, .5 * dt, index, 'first_half')
                if first['converged']:
                    second = evolve(first, .5 * dt, index, 'second_half')
                else: second = first
            else: first = second = full
            if not all(s['converged'] for s in (full, first, second)):
                entry.update(accepted_numerically=False, equation_failure=next(s['message'] for s in (full, first, second) if not s['converged']))
                counts['equation_rejections'] += 1; dt *= .5
            else:
                comparison = assess_time_error(full, second)
                entry.update({k: v for k, v in comparison.items() if k != 'structure_difference_by_node'})
                entry['accepted_numerically'] = comparison['meets_driver_time_accuracy']
                if not entry['accepted_numerically']:
                    counts['numerical_rejections'] += 1; dt *= .5
                else:
                    elapsed += dt; current = second
                    accepted.append(save_state(current, dt, comparison['error_norm']))
                    dt *= float(np.clip(.9 / math.sqrt(max(comparison['error_norm'], 1e-6)), .5, 2))
            print(json.dumps(entry), flush=True)
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
    reached = elapsed >= args.duration_years * YEAR
    outcome = ('failed_conditional_adaptive_control' if failures else
               'completed_conditional_adaptive_interval' if reached else 'conditional_attempt_limit_reached')
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(), outcome=outcome,
                  accepted_for_stellar_evolution=False, selected_star_unchanged=True,
                  target_reached=reached, elapsed_years=elapsed / YEAR, target_years=args.duration_years,
                  elapsed_run_seconds=time.monotonic() - start, native_exit=proc.returncode,
                  **counts, new_accepted_intervals=len(accepted), reused_initial_intervals=0 if args.resume_checkpoint else 1,
                  resumed_checkpoint=str(args.resume_checkpoint) if args.resume_checkpoint else None,
                  initial_elapsed_years=initial_elapsed / YEAR,
                  failures=failures, attempts=attempted, accepted_history=accepted,
                  artifacts_sha256={str(p): digest(p) for p in args.scratch.iterdir()},
                  limitations=['Conditional opacity and envelope heat assumptions remain unaccepted.',
                               'This measures local adaptive time accuracy, not accumulated lifetime error.',
                               'The initial adjustment includes an EOS and transport change.'])
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: result[k] for k in ('outcome', 'elapsed_years', 'native_steps', 'reused_trials', 'new_accepted_intervals', 'failures')}), flush=True)
    if failures: raise SystemExit(1)


if __name__ == '__main__':
    main()
