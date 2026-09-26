#!/usr/bin/env python3
"""Adaptive stellar evolution with physical metal transport and opacity feedback.

Each accepted interval uses two half steps, checked against one full step.
Save the complete isotope state and independent conservation measurements.
Material selections are explicit in a JSON configuration and copied binaries
and source identities accompany each segment. Domain failures stop progress.
"""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import time
import numpy as np
from audit_cn_thermal import physical as fixed_metal_physical
from audit_metal_evolution import YEAR, physical, weights, audit, compare, abundance_step_change
from write_scientific_result import write_result


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def restore(data, convert):
    record = data.get('model_record', data)
    model = np.array(record['model'])
    age = data.get('age_seconds', record.get('age_seconds'))
    if age is None or not np.isfinite(age) or age < 0:
        raise ValueError('restart requires its absolute age')
    if model.shape == (512, 10):
        if not convert:
            raise ValueError('fixed-metal material restart requires explicit physical conversion')
        z = .02+(model[:, 7:10]*[12, 13, 14]).sum(axis=1)-.02*(.171836+.050335)
        converted = np.column_stack((model, z))
        if abs(physical(converted)[:, :6]-fixed_metal_physical(model)).max() > 3e-16:
            raise ValueError('restart conversion changed physical isotope inventories')
        model = converted
    elif model.shape != (512, 11) or convert:
        raise ValueError('invalid physical restart or redundant conversion request')
    if not np.isfinite(model).all() or physical(model).min() < 0:
        raise ValueError('invalid restart composition')
    return model, age, record


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--configuration', type=Path, required=True)
    ap.add_argument('--checkpoint', type=Path, required=True)
    ap.add_argument('--convert-fixed-metal', action='store_true')
    ap.add_argument('--scratch', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--years', type=float, required=True)
    ap.add_argument('--initial-step-years', type=float, default=1e8)
    ap.add_argument('--max-intervals', type=int, default=5000)
    a = ap.parse_args()
    if not np.isfinite([a.years, a.initial_step_years]).all() or min(a.years, a.initial_step_years) <= 0:
        raise ValueError('finite positive time intervals required')
    if a.max_intervals < 1 or a.output.exists():
        raise ValueError('invalid interval limit or existing report')
    config = json.loads(a.configuration.read_text())
    command = list(config['native_command'])
    if (len(command) not in (13, 15, 17, 18, 19, 20, 21, 22, 23, 24)
            or command[11] not in ('variable-opacity', 'variable-opacity-hydrogen-share')):
        raise ValueError('explicit variable interior and atmospheric metal response required')
    hydrogen_share_opacity = command[11] == 'variable-opacity-hydrogen-share'
    if len(command) >= 15 and command[14] not in ('oplib-ratio', 'linear-kappa'):
        raise ValueError('unknown upper-metal opacity extension')
    if len(command) >= 17 and command[16] not in ('linear-kappa', 'log-kappa'):
        raise ValueError('unknown lower-metal opacity continuation')
    if len(command) >= 21:
        maximum_warm_hydrogen = float(command[20])
        allowed_hydrogen = .99999 if hydrogen_share_opacity else .97
        if not np.isfinite(maximum_warm_hydrogen) or not .90 <= maximum_warm_hydrogen <= allowed_hydrogen:
            raise ValueError('warm hydrogen continuation must be limited to [.90,.99999]'
                             if hydrogen_share_opacity else
                             'warm hydrogen continuation must be limited to [.90,.97]')
    bridge_hydrogen_method = command[21] if len(command) >= 22 else 'oplib-ratio'
    if bridge_hydrogen_method not in ('oplib-ratio', 'source-log'):
        raise ValueError('unknown bridge hydrogen opacity method')
    atmosphere_reference_mode = command[22] if len(command) >= 23 else 'nested-reference'
    if atmosphere_reference_mode not in ('nested-reference', 'low-z-reference', 'continuous-reference', 'trace-metal-interval'):
        raise ValueError('unknown atmosphere reference mode')
    if (atmosphere_reference_mode in ('continuous-reference', 'trace-metal-interval')) != (len(command) == 24):
        raise ValueError('continuous atmosphere mode requires exactly one interval manifest')
    maximum_metal = float(command[17]) if len(command) >= 18 else .05
    minimum_metal = float(command[18]) if len(command) >= 19 else .004
    if (not np.isfinite([minimum_metal, maximum_metal]).all()
            or not 0 < minimum_metal < .01 < .03 < maximum_metal < 1):
        raise ValueError('invalid opacity metallicity interval')
    relative = config['composition_relative_tolerance']
    if not np.isfinite(relative) or not 0 < relative <= .01:
        raise ValueError('invalid composition time-error tolerance')
    requested_abundance_tolerance = config.get('preferred_abundance_tolerance')
    if requested_abundance_tolerance is not None and (
            not np.isfinite(requested_abundance_tolerance)
            or not 1e-12 <= requested_abundance_tolerance <= 1e-10):
        raise ValueError('preferred abundance solver tolerance must be in [1e-12,1e-10]')
    structure_relative = config.get('structure_relative_tolerance', 1e-4)
    helium3_relative = config.get('helium3_relative_tolerance', relative)
    if (not np.isfinite([structure_relative, helium3_relative]).all()
            or not 0 < structure_relative <= .01 or not 0 < helium3_relative <= .01):
        raise ValueError('invalid structure or helium-3 time-error tolerance')
    # Record effective settings even for configurations using old defaults.
    config.update(structure_relative_tolerance=structure_relative,
                  helium3_relative_tolerance=helium3_relative)
    if config['maximum_abundance_change'] != .005:
        raise ValueError('abundance proposal must agree with the selected native guard')
    threads = int(command[10])
    if not 1 <= threads <= 12:
        raise ValueError('invalid native thread allocation')
    if not config.get('limitations'):
        raise ValueError('explicit physical approximations required')
    core_trial = config.get('core_disappearance_trial_years', 0.)*YEAR
    core_mass_limit = config.get('core_disappearance_mass_fraction', 0.)
    if (not np.isfinite([core_trial, core_mass_limit]).all() or core_trial < 0
            or (core_trial > 0 and not 0 < core_mass_limit <= 1e-4)):
        raise ValueError('invalid core-disappearance proposal')
    a.scratch.mkdir(exist_ok=False)
    sources = a.scratch/'sources'; sources.mkdir()
    paths = [a.configuration, a.checkpoint, Path(__file__),
             Path(__file__).with_name('audit_metal_evolution.py'),
             Path(__file__).with_name('audit_cn_thermal.py'),
             *map(Path, command[:10]), Path(command[12]),
             *map(Path, config['source_paths'])]
    if len(command) >= 15:
        paths.append(Path(command[13]))
    if len(command) >= 17:
        paths.append(Path(command[15]))
    if len(command) >= 20:
        paths.append(Path(command[19]))
    if len(command) >= 24:
        paths.append(Path(command[23]))
    # Retain every referenced EOS and opacity body, including the warm,
    # molecular and hot mixture families, not just the hydrogen correction.
    pending = list(map(Path, command[1:10])); seen = set()
    pending.append(Path(command[12]))
    if len(command) >= 17:
        pending.append(Path(command[15]))
    if len(command) >= 20:
        pending.append(Path(command[19]))
    if len(command) >= 24:
        pending.append(Path(command[23]))
    if len(command) >= 15:
        pending.append(Path(command[13]))
    while pending:
        manifest = pending.pop().resolve(strict=True)
        if manifest in seen: continue
        seen.add(manifest)
        with manifest.open() as stream: header = stream.readline()
        if header.startswith(('EMBER_VARIABLE_METAL_HELMHOLTZ', 'EMBER_OPACITY_MIXTURE',
                              'EMBER_METAL_ATMOSPHERE_CHAIN',
                              'EMBER_HYDROGEN_ATMOSPHERE_INTERVAL')):
            for line in manifest.read_text().splitlines()[1:]:
                if '"' in line:
                    child = (manifest.parent/json.loads(line[line.index('"'):])).resolve(strict=True)
                    paths.append(child); pending.append(child)
    hashes = {str(p.resolve()): digest(p) for p in dict.fromkeys(paths)}
    for p in [Path(__file__), Path(__file__).with_name('audit_metal_evolution.py'),
              Path(__file__).with_name('audit_cn_thermal.py'), *map(Path, config['source_paths'])]:
        shutil.copy2(p, sources/p.name)
    native = a.scratch/'native'; shutil.copy2(command[0], native); command[0] = str(native.resolve())
    data = json.loads(a.checkpoint.read_text()); old, base_age, previous = restore(data, a.convert_fixed_metal)
    initial = old.copy(); age = base_age; w = weights(old[:, 0]); mass = old[-1, 0]
    initial_species = w@physical(initial)/mass
    write_result(a.scratch/'plan.json', dict(configuration=config, input_sha256=hashes,
        native_command=command, initial_age_seconds=age, requested_years=a.years,
        converted_fixed_metal=a.convert_fixed_metal, initial_global_species=initial_species.tolist()))
    elapsed = 0.; target = a.years*YEAR; h = a.initial_step_years*YEAR
    preferred = data.get('next_abundance_tolerance', 1e-12)
    fully_mixed = previous.get('mixing_regions') == [[0, len(old)]]
    if fully_mixed: preferred = 1e-12
    elif requested_abundance_tolerance is not None:
        preferred = requested_abundance_tolerance
    accumulated = np.zeros((len(old), 7)); history = []; attempts = []
    final_record = None; count = 0; solver_seconds = 0.; solver_cpu = 0.
    started_utc = datetime.now(timezone.utc); start = time.monotonic()
    child_start = resource.getrusage(resource.RUSAGE_CHILDREN)
    stderr = (a.scratch/'native.stderr').open('x')
    environment = os.environ.copy()
    for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS']:
        environment[key] = '1'
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=stderr, text=True, env=environment)
    (a.scratch/'native_pid.txt').write_text(str(process.pid)+'\n')
    def step(model, absolute_age, interval, tolerance):
        nonlocal count, solver_seconds, solver_cpu
        values = [absolute_age, interval, tolerance, 1e-6, len(model), *model.flat]
        query = ' '.join(format(v, '.17g') for v in values)+'\n'
        process.stdin.write(query); process.stdin.flush(); raw = process.stdout.readline(); count += 1
        if not raw: raise RuntimeError('native process ended without a result')
        result = json.loads(raw)
        solver_seconds += result.get('seconds', 0); solver_cpu += result.get('cpu_seconds', 0)
        try:
            if not result.get('variable_interior_opacity') or not result.get('variable_atmospheric_response'):
                raise ValueError(result.get('error', 'requested opacity feedback is absent'))
            if bool(result.get('hydrogen_share_opacity', False)) != hydrogen_share_opacity:
                raise ValueError('requested hydrogen-share opacity selection differs')
            if bool(result.get('upper_metal_opacity_extension')) != (len(command) >= 15):
                raise ValueError('requested upper-metal opacity selection differs')
            if (bool(result.get('lower_metal_opacity_continuation')) != (len(command) >= 17)
                    or bool(result.get('lower_metal_atmospheric_response')) != (len(command) >= 17)):
                raise ValueError('requested lower-metal material selection differs')
            if bool(result.get('depleted_metal_atmospheric_response')) != (len(command) >= 20):
                raise ValueError('requested depleted-metal atmosphere selection differs')
            if (result.get('maximum_opacity_metallicity') != maximum_metal
                    or result.get('minimum_opacity_metallicity') != minimum_metal):
                raise ValueError('requested opacity metallicity interval differs')
            if len(command) >= 21 and result.get('maximum_warm_hydrogen') != maximum_warm_hydrogen:
                raise ValueError('requested warm hydrogen interval differs')
            if len(command) >= 22 and result.get('bridge_hydrogen_method') != bridge_hydrogen_method:
                raise ValueError('requested bridge hydrogen opacity method differs')
            if len(command) >= 23 and result.get('atmosphere_reference_mode') != atmosphere_reference_mode:
                raise ValueError('requested atmosphere reference mode differs')
            new, conservation, balance = audit(model, result, absolute_age, interval)
        except ValueError:
            with gzip.open(a.scratch/f'failed-{count:05d}.json.gz', 'wt') as stream:
                json.dump(dict(query=query, response=result), stream)
            raise
        return new, result, conservation, balance
    def checkpoint():
        payload = dict(accepted_for_stellar_evolution=False, composition_layout='physical_metal_11',
            age_seconds=age, segment_elapsed_seconds=elapsed, model_record=final_record,
            next_step_seconds=h, next_abundance_tolerance=preferred,
            initial_source=str(a.checkpoint.resolve()), input_sha256=hashes,
            numerical_history_row=history[-1], accumulated_species_balance=accumulated.tolist(),
            limitations=config['limitations'])
        temporary = a.scratch/'checkpoint.tmp'; temporary.write_text(json.dumps(payload)+'\n')
        temporary.replace(a.scratch/'checkpoint.json')
    outcome = 'running'; failure = None; core_trial_age = None
    profiles = gzip.open(a.scratch/'accepted_models.jsonl.gz', 'wt', compresslevel=1)
    try:
        while elapsed < target and len(history) < a.max_intervals:
            h = min(h, target-elapsed)
            item = dict(start_age_seconds=age); accepted = False
            # A node-count boundary rule can drive the last, negligible core
            # toward zero timestep. Try crossing its disappearance at finite
            # dt once per accepted state. All full/half, local composition,
            # fuel and energy criteria remain in force, including the
            # pointwise fallback when the number of convective regions changes.
            regions = previous.get('mixing_regions', [])
            if core_trial > h and core_trial_age != age and regions:
                begin, end = regions[0]
                if begin == 0 and end > 1 and w[:end].sum()/mass <= core_mass_limit:
                    item['core_disappearance_proposal'] = dict(
                        previous_interval_years=h/YEAR,
                        core_mass_fraction=float(w[:end].sum()/mass))
                    h = min(core_trial, target-elapsed); core_trial_age = age
            item['interval_years'] = h/YEAR
            tolerances = sorted({preferred, 1e-12}, reverse=True)
            for tolerance in tolerances:
                try:
                    _, full, c1, _ = step(old, age, h, tolerance)
                    mid, half1, c2, b1 = step(old, age, h/2, tolerance)
                    fine, half2, c3, b2 = step(mid, half1['age_seconds'], h/2, tolerance)
                    diagnostics = {}
                    error, metrics = compare(full, half2, 'moving-boundary', diagnostics,
                        previous, half1, relative, structure_relative, helium3_relative)
                    item.update(error_norm=error, criteria=metrics, tolerance=tolerance,
                        conservation=[c1, c2, c3], comparison_diagnostics=diagnostics)
                    accepted = error <= 1
                    break
                except ValueError as e:
                    item.setdefault('solve_failures', []).append(dict(tolerance=tolerance, message=str(e)))
                    # A smaller step addresses failed coupling convergence.
                    # Tightening its convergence target can repeat a failed
                    # large step. Retain the stricter retry for conservation or
                    # other failures of a converged solve.
                    if (config.get('skip_stricter_retry_after_coupling_limit', False)
                            and str(e) == 'evolve_step: coupling iteration limit'):
                        item['stricter_retry_skipped'] = 'reduce timestep after coupling iteration limit'
                        break
            item['accepted'] = accepted; attempts.append(item)
            with (a.scratch/'attempts.jsonl').open('a') as stream: stream.write(json.dumps(item)+'\n')
            if not accepted:
                if len(attempts) >= 4 and all(t.get('solve_failures') and not t['accepted'] for t in attempts[-4:]):
                    messages = [t['solve_failures'][-1]['message'] for t in attempts[-4:]]
                    if len(set(messages)) == 1:
                        raise RuntimeError('unchanged rejection after four reduced steps: '+messages[0])
                h *= max(.1, min(.5, .8/np.sqrt(item.get('error_norm', 4.))))
                if h < YEAR or (len(attempts) >= 30 and not any(t['accepted'] for t in attempts[-30:])):
                    raise RuntimeError('repeated rejected steps: '+str(item))
                continue
            change = abundance_step_change(old, fine); old = fine; elapsed += h
            age = half2['age_seconds']; accumulated += b1+b2; final_record = half2; previous = half2
            now_mixed = half2['mixing_regions'] == [[0, len(old)]]
            preferred = 1e-10 if fully_mixed and not now_mixed else tolerance
            # A stricter retry is local to its interval. An explicitly selected
            # normal target is retried on the next stratified stellar model;
            # independent time and conservation checks still decide acceptance.
            if not now_mixed and requested_abundance_tolerance is not None:
                preferred = requested_abundance_tolerance
            fully_mixed = now_mixed
            global_species = w@physical(fine)/mass
            row = dict(age_seconds=age, age_years=age/YEAR, elapsed_years=elapsed/YEAR,
                interval_years=h/YEAR, error_norm=error, Teff=half2['Teff'],
                central_T=float(fine[0, 3]), central_H=float(fine[0, 5]), surface_H=float(fine[-1, 5]),
                central_Z=float(fine[0, 10]), surface_Z=float(fine[-1, 10]),
                luminosity=float(fine[-1, 4]), radius=float(fine[-1, 1]),
                nuclear_luminosity=half2['nuclear_luminosity'],
                mixed_mass_fraction=half2['convective_mass_fraction'], fully_mixed=now_mixed,
                convective_regions=[r for r in half2['mixing_regions'] if r[1] > r[0]+1],
                maximum_abundance_step_change=change, global_species=global_species.tolist(),
                accumulated_global_species_balance=accumulated.sum(axis=0).tolist(),
                total_inert_mass_fraction_change=float(global_species[-1]-initial_species[-1]),
                represented_age_difference_seconds=age-(base_age+elapsed))
            factor = min(2., max(.5, .85/np.sqrt(max(error, 1e-6))))
            factor = min(factor, .85*.005/max(change, 1e-30))
            history.append(row); h *= factor; checkpoint()
            profiles.write(json.dumps(dict(age_seconds=age, model_record=half2))+'\n'); profiles.flush()
            with (a.scratch/'history.jsonl').open('a') as stream: stream.write(json.dumps(row)+'\n')
            print(f'accepted {len(history)}: age {age/YEAR/1e12:.4g} Tyr, Teff {half2["Teff"]:.4g} K, Xs {fine[-1,5]:.4g}, Zs {fine[-1,10]:.4g}, error {error:.4g}', flush=True)
            threshold = config.get('stop_surface_hydrogen')
            if threshold is not None and fine[-1, 5] >= threshold:
                outcome = 'surface_hydrogen_review_trigger'; break
            if (a.scratch/'stop.request').exists(): outcome = 'requested_checkpoint_stop'; break
        else: outcome = 'requested_interval_reached' if elapsed >= target else 'interval_count_limit'
    except Exception as e:
        outcome = 'stopped'; failure = str(e); print('STOP', failure, flush=True)
    finally:
        profiles.close(); process.stdin.close(); code = process.wait(timeout=60); stderr.close()
    changed_inputs = [p for p, expected in hashes.items() if digest(p) != expected]
    if changed_inputs: outcome = 'changed_inputs'; failure = str(changed_inputs)
    child_end = resource.getrusage(resource.RUSAGE_CHILDREN)
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(), outcome=outcome, failure=failure,
        accepted_for_stellar_evolution=False, initial_age_seconds=base_age, final_age_seconds=age,
        elapsed_years=elapsed/YEAR, requested_years=a.years, accepted_intervals=len(history),
        rejected_intervals=sum(not t['accepted'] for t in attempts), native_solves=count,
        solver_seconds=solver_seconds, solver_cpu_seconds=solver_cpu,
        wall_seconds=time.monotonic()-start,
        elapsed_utc_seconds=(datetime.now(timezone.utc)-started_utc).total_seconds(),
        started_utc=started_utc.isoformat(),
        timing_scope='Python monotonic time excludes system sleep on macOS; UTC elapsed includes it. Native steady-clock and OS CPU timing are reported separately.',
        native_total_cpu_seconds=child_end.ru_utime+child_end.ru_stime-child_start.ru_utime-child_start.ru_stime,
        native_exit_code=code, zone_threads=threads, configuration=config,
        checkpoint=str(a.scratch/'checkpoint.json') if history else None,
        history=history, accumulated_global_species_balance=accumulated.sum(axis=0).tolist(),
        input_sha256=hashes, artifacts_sha256={str(p):digest(p) for p in a.scratch.iterdir() if p.is_file()},
        limitations=config['limitations'])
    write_result(a.output, report)
    return int(outcome in ['stopped', 'changed_inputs'] or code != 0)


if __name__ == '__main__':
    raise SystemExit(main())
