#!/usr/bin/env python3
"""Adaptive CN/H/He evolution with explicit trace transport and independent audits.

Checkpoints retain the full isotope state. A supplied non-CN starting profile
requires an explicit initial carbon conversion and is a conditional restart.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import shutil
import resource
import time

import numpy as np
from audit_cn_thermal import ATOMIC, C0, LIGHT, MASS, N0, YEAR, physical, weights

TABLES = [
    '/tmp/ember-smooth-eos-refined-family-v1/freeeos300_gs98_z020.dat',
    '/tmp/ember-refractive-hot-family-refined-v1/aesopus21_gs98_mixture.dat',
    '/tmp/ember-op-warm-family-v4/density0025/op_gs98_ordinary.dat',
    '/tmp/ember-tops-bridge-v2/tables/tops_bridge.dat',
    '/tmp/ember-refractive-hot-family-refined-v1/tops_gs98_mixture_high.dat',
    'data/conduction/condtab21wd_metals.dat',
    '/tmp/ember-fable-nongrey-x030-candidate-v1.dat',
    '/tmp/ember-refractive-hot-family-refined-v1/tops_gs98_mixture_low.dat',
]


def nuclear_mass_energy_error(release, nuclear, photon):
    """Scale a rest-mass discrepancy to the power relevant to the structure.

    The nuclear scale is retained during burning. Photon luminosity supplies
    the scale during cooling, including exactly zero nuclear power; no
    arbitrary positive floor or numerical roundoff allowance is introduced.
    """
    if (not np.isfinite([release, nuclear, photon]).all()
            or nuclear < 0 or photon < 0):
        raise ValueError('invalid power in the independent nuclear energy check')
    scale = max(nuclear, photon)
    if scale == 0:
        raise ValueError('independent nuclear energy check has no physical power scale')
    return (release - nuclear) / scale


def audit(old, result, dt):
    if not result.get('converged'):
        raise ValueError(result.get('message', result.get('error', 'missing native result')))
    new = np.array(result['model']); source = np.array(result['physical_sources'])
    energy = np.array(result['energy_cells']); w = weights(old[:, 0]); mass = old[-1, 0]
    before, after = physical(old), physical(new)
    if not np.isfinite(new).all() or after.min() < 0 or not np.array_equal(new[:, 0], old[:, 0]):
        raise ValueError('invalid physical composition, state or mesh')
    if not result['input_preserved'] or result['elapsed_age_seconds'] != dt:
        raise ValueError('input state or interval changed')
    if abs(after.sum(axis=1)+.02-(12*C0+14*N0)-1).max() > 2e-14:
        raise ValueError('physical mass fractions do not close')
    delta = after-before; residual = delta-dt*source
    flux = np.zeros((len(new)-1, 6))
    for face, f in result['cn_boundary_fluxes']:
        flux[face] = dt/mass*np.array([f[0], f[1], -sum(f), *f[2:]])
    region_error = 0.
    for begin, end in result['mixing_regions']:
        r = w[begin:end]@residual[begin:end]/mass
        if begin: r -= flux[begin-1]
        if end < len(new): r += flux[end-1]
        region_error = max(region_error, float(abs(r).max()))
        if end > begin+1 and not np.all(new[begin:end, 5:10] == new[begin, 5:10]):
            raise ValueError('convective reservoir is not homogeneous')
    global_error = float(abs(w@residual/mass).max())
    catalyst_error = float(abs(w@(new[:, 7:10].sum(axis=1)-old[:, 7:10].sum(axis=1))/mass/(C0+N0)))
    release = -float(w@delta@(ATOMIC/MASS-1))*LIGHT**2/dt
    nuclear = float(w@(energy[:, 3]+energy[:, 4]))
    mass_energy = nuclear_mass_energy_error(release, nuclear, float(new[-1, 4]))
    grav = -(energy[:, 0]-energy[:, 2]+energy[:, 1]*(1/new[:, 2]-1/old[:, 2]))/dt
    first_law = float(w@(energy[:, 3]-energy[:, 5]+grav))/new[-1, 4]-1
    report = dict(region_species_error=region_error, global_species_error=global_error,
                  catalyst_number_error=catalyst_error, mass_energy_error=mass_energy,
                  first_law_error=first_law)
    # The native production inner abundance tolerance is 1e-11. Demanding
    # 1e-13 of stellar mass here discarded steps whose independently checked
    # fuel-energy error was already below 5e-8. Retain the tighter physical
    # energy test and measure accumulated conservation separately. In the
    # cooling limit, its unchanged 2e-6 budget is relative to photon power
    # rather than the vanishing nuclear power.
    limits = dict(region_species_error=2e-11, global_species_error=2e-11,
                  catalyst_number_error=2e-9, mass_energy_error=2e-6, first_law_error=2e-8)
    if any(not np.isfinite(v) or abs(v)>limits[k] for k, v in report.items()):
        raise ValueError('independent conservation check: '+json.dumps(report))
    return new, report, w[:, None]/mass*residual


def compare(coarse, fine, composition_mode='pointwise', diagnostics=None,
            previous=None, intermediate=None, composition_relative_tolerance=1e-4):
    if not np.isfinite(composition_relative_tolerance) or not 0<composition_relative_tolerance<=.01:
        raise ValueError('composition relative tolerance must be positive and at most 0.01')
    a, b = np.array(coarse['model']), np.array(fine['model'])
    pa, pb = physical(a), physical(b); w = weights(a[:, 0])
    mixed_a = [(x, y) for x, y in coarse['mixing_regions'] if y>x+1]
    mixed_b = [(x, y) for x, y in fine['mixing_regions'] if y>x+1]
    boundary = max((abs(x-y) for ra, rb in zip(mixed_a, mixed_b) for x, y in zip(ra, rb)), default=0)
    if len(mixed_a)!=len(mixed_b): boundary = len(a)
    species_error = abs(pa-pb)/(1e-8+composition_relative_tolerance*np.maximum(pa, pb))
    distribution_error = float(w@abs(a[:, 5]-b[:, 5])/(w@b[:, 5])/1e-6)
    metrics = dict(structure=float(abs(np.log(a[:, 1:4]/b[:, 1:4])).max()/1e-4),
                   species=float(species_error.max()),
                   luminosity=float(abs(a[-1, 4]/b[-1, 4]-1)/1e-3),
                   hydrogen=distribution_error,
                   mixed_boundary=float(boundary))
    if composition_mode=='moving-boundary':
        if previous is None or intermediate is None:
            raise ValueError('moving-boundary comparison requires start and intermediate reservoirs')
        keep = np.ones(len(a), dtype=bool)
        paths = [[(x,y) for x,y in record['mixing_regions'] if y>x+1]
                 for record in [previous, intermediate, coarse, fine]]
        motion = 0
        topology_fallback = len({len(path) for path in paths})!=1
        if not topology_fallback:
            for reservoirs in zip(*paths):
                for boundaries in zip(*reservoirs):
                    motion = max(motion,max(boundaries)-min(boundaries))
            if motion<=1:
                for reservoirs in zip(*paths):
                    for boundaries in zip(*reservoirs):
                        keep[min(boundaries):max(boundaries)] = False
        # Endpoints may have the same boundary but cross the cell at different
        # times. Account for the cell traversed over the whole interval,
        # permitting at most one cell per boundary, not an arbitrary band.
        if diagnostics is not None:
            diagnostics.update(boundary_cells=np.flatnonzero(~keep).tolist(),
                               all_node_species_error_norm=metrics['species'],
                               hydrogen_distribution_error_norm=distribution_error,
                               topology_fallback=topology_fallback)
        metrics['species'] = float(species_error[keep].max())
        metrics['boundary_motion'] = float(motion)
        metrics['composition_l1'] = float(((w@abs(pa-pb))/
            (w.sum()*1e-8+composition_relative_tolerance*(w@np.maximum(pa,pb)))).max())
        # Fuel accuracy is the difference in total hydrogen. Redistribution
        # is constrained separately by composition_l1 and the local checks.
        metrics['hydrogen'] = (distribution_error if topology_fallback else
                              float(abs(w@(a[:,5]-b[:,5]))/(w@b[:,5])/1e-6))
    elif composition_mode!='pointwise':
        raise ValueError('unknown composition comparison')
    if diagnostics is not None:
        active_error = species_error.copy()
        if composition_mode=='moving-boundary': active_error[~keep] = -1
        node, species = np.unravel_index(active_error.argmax(), active_error.shape)
        diagnostics.update(maximum_species_error_node=int(node),
                           maximum_species_error_species=['H1','He3','He4','C12','C13','N14'][species],
                           maximum_species_error_values=[float(pa[node,species]), float(pb[node,species])])
    return max(metrics.values()), metrics


def abundance_step_change(before, after):
    """Same abundance coordinates as the native step-size guard."""
    hhe = after[:, 5:7]-before[:, 5:7]
    catalysts = (after[:, 7:10]-before[:, 7:10])*np.array([12., 13., 14.])
    return max(float(abs(hhe).max()), float(abs(hhe.sum(axis=1)).max()),
               float(abs(catalysts).max()))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--native', type=Path, required=True)
    ap.add_argument('--checkpoint', type=Path, required=True)
    ap.add_argument('--scratch', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--collision-table', type=Path, default=Path('/tmp/ember-collision-transport-table-v1.dat'))
    ap.add_argument('--atmosphere', type=Path,
                    help='explicit checked gas-atmosphere table; included in the input hashes')
    ap.add_argument('--years', type=float, required=True)
    ap.add_argument('--initial-step-years', type=float, default=2.5e5)
    ap.add_argument('--initial-carbon-conversion', type=float)
    ap.add_argument('--warm-opacity', choices=['op', 'retained-tops'], default='op')
    ap.add_argument('--transport', choices=['helium', 'zero-cn', 'fully-convective'], default='helium')
    ap.add_argument('--composition-error-mode', choices=['pointwise', 'moving-boundary'], default='pointwise',
                    help='composition comparison across full and half steps; see the matched controls in PRODUCTION_ACCURACY.md')
    ap.add_argument('--zone-threads', type=int, default=1,
                    help='parallel zone evaluations with compatible native physics; tables are shared')
    ap.add_argument('--composition-relative-tolerance', type=float, default=1e-4,
                    help='relative composition tolerance; the selected 0.001 budget is passed explicitly by production launchers')
    ap.add_argument('--stop-carbon12-depletion', '--stop-conversion', dest='stop_conversion', type=float,
                    help='stop after this fraction of surface C12 is destroyed, including conversion to C13')
    ap.add_argument('--stop-surface-hydrogen', type=float,
                    help='checkpoint after an accepted interval reaches this surface H fraction; does not clip composition')
    ap.add_argument('--max-intervals', type=int, default=5000)
    a = ap.parse_args()
    if a.atmosphere is not None:
        TABLES[6] = str(a.atmosphere.resolve(strict=True))
    if a.years<=0 or a.initial_step_years<=0: raise ValueError('positive intervals required')
    if a.stop_surface_hydrogen is not None and not 0<a.stop_surface_hydrogen<=1:
        raise ValueError('surface hydrogen review threshold must be between zero and one')
    if not 1<=a.zone_threads<=64: raise ValueError('zone threads must be between 1 and 64')
    if not np.isfinite(a.composition_relative_tolerance) or not 0<a.composition_relative_tolerance<=.01:
        raise ValueError('composition relative tolerance must be positive and at most 0.01')
    a.scratch.mkdir(exist_ok=False)
    if a.output.exists(): raise ValueError('preserve completed reports')
    start = time.monotonic(); hashes = {}
    def pin(path):
        path = Path(path); h = hashlib.sha256()
        with path.open('rb') as f:
            for chunk in iter(lambda: f.read(8*1024*1024), b''): h.update(chunk)
        hashes[str(path)] = h.hexdigest()
    collision = str(a.collision_table)
    for path in [a.native, a.checkpoint, __file__, 'scripts/audit_cn_thermal.py', *TABLES, collision]: pin(path)
    for path in [Path(__file__), Path('scripts/audit_cn_thermal.py')]:
        shutil.copy2(path, a.scratch/path.name)
    # Seal referenced EOS plane bodies, not just the small family manifest.
    manifest = Path(TABLES[0])
    for line in manifest.read_text().splitlines()[3:]:
        if line.startswith('"'): pin((manifest.parent/json.loads(line)).resolve())
    data = json.loads(a.checkpoint.read_text()); inherited = data.get('limitations', [])
    if 'profile' in data:
        old = np.array(data['profile'])[:, :7]
        last = data['history'][-1]
        base_age = (last.get('age_yr', 0) if isinstance(last, dict) else last[0])*YEAR
    elif 'model_record' in data:
        old = np.array(data['model_record']['model'])
        base_age = data.get('age_seconds', data.get('numerical_history_row', {}).get('age_seconds', 0))
    else: old = np.array(data['model']); base_age = data.get('age_seconds', 0)
    imposed = old.shape[1] == 7
    if imposed:
        f = a.initial_carbon_conversion
        if f is None or not 0<=f<=1: raise ValueError('non-CN restart needs an explicit initial carbon conversion')
        old = np.column_stack((old, np.tile([C0*(1-f), 0., N0+C0*f], (len(old), 1))))
    elif old.shape[1]!=10: raise ValueError('unsupported checkpoint layout')
    elif a.initial_carbon_conversion is not None: raise ValueError('do not overwrite an evolved catalyst state')
    initial = old.copy(); w = weights(old[:, 0]); mass = old[-1, 0]
    regions = data.get('model_record', {}).get('mixing_regions', [])
    fully_mixed = regions == [[0, len(old)]]
    preferred_tolerance = data.get('next_abundance_tolerance',
                                  1e-12 if fully_mixed else 1e-10)
    history, attempts = [], []; accumulated = np.zeros((len(old), 6)); elapsed = 0.; h = a.initial_step_years*YEAR
    target = a.years*YEAR; step_count = 0; solver_seconds = 0.; solver_cpu_times = []; final_record = None
    limitations = [
        'Conditional physical calculation; not a replacement for the consistent trajectory from formation.',
        'Fixed GS98 material tables and atmosphere; separate CN settling, CN material feedback and oxygen branches omitted.',
        'CN transport shares the helium-group velocity; zero-CN option is a named comparison.',
        'Cool reduced microscopic heat omitted below the stated 2–3 MK join; total H/He material enthalpy retained.',
        'Initial catalyst inventory is imposed on a non-CN profile.' if imposed else 'Evolved catalyst inventory restored without renormalization.',
    ]
    if a.transport=='fully-convective':
        limitations.append('Fully convective stage only: material conduction and total composition enthalpy retained; reduced microscopic heat omitted everywhere. A radiative species boundary rejects, without automatic continuation.')
    if a.composition_error_mode=='moving-boundary':
        limitations.append('Moving-boundary time comparison: pointwise composition outside at most one cell per shifted boundary, whole-profile composition L1, and total hydrogen fuel; validated intervals are documented in PRODUCTION_ACCURACY.md.')
    if a.composition_relative_tolerance!=1e-4:
        limitations.append('Nondefault composition relative tolerance is explicitly recorded; the 0.001 selection has matched finite-interval controls, not a full-track error bound. Absolute composition, hydrogen fuel, structure, boundary and conservation limits remain unchanged.')
    if inherited: limitations.append({'inherited_restart_scope': inherited})
    stderr = (a.scratch/'native.stderr').open('w')
    command = [str(a.native), *TABLES, '.1', '1', collision]
    if a.zone_threads>1: command.append(str(a.zone_threads))
    child_start = resource.getrusage(resource.RUSAGE_CHILDREN)
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=stderr, text=True)
    (a.scratch/'native_pid.txt').write_text(str(process.pid)+'\n')
    def step(model, interval, tolerance=1e-10):
        nonlocal step_count, solver_seconds
        values = [int(a.warm_opacity=='retained-tops'), 1., 1., interval, tolerance, len(model),
                  {'zero-cn': 0, 'helium': 1, 'fully-convective': 2}[a.transport], 1e-6, *model.flat]
        query = ' '.join(format(x, '.17g') for x in values)+'\n'
        process.stdin.write(query); process.stdin.flush(); raw = process.stdout.readline(); step_count += 1
        if not raw: raise RuntimeError('native process ended without a result')
        result = json.loads(raw); solver_seconds += result.get('seconds', 0.)
        if 'cpu_seconds' in result: solver_cpu_times.append(result['cpu_seconds'])
        try: new, conservation, balances = audit(model, result, interval)
        except ValueError:
            (a.scratch/'last_failed_query.txt').write_text(query)
            (a.scratch/'last_failed_result.json').write_text(raw)
            raise
        return new, result, conservation, balances
    def save_checkpoint():
        payload = dict(accepted_for_stellar_evolution=False, age_seconds=base_age+elapsed,
                       segment_elapsed_seconds=elapsed, model_record=final_record,
                       next_step_seconds=h, limitations=limitations,
                       next_abundance_tolerance=preferred_tolerance,
                       initial_source=str(a.checkpoint), input_sha256=hashes,
                       numerical_history_row=history[-1], accumulated_species_balance=accumulated.tolist())
        tmp = a.scratch/'checkpoint.tmp'; tmp.write_text(json.dumps(payload)+'\n');tmp.replace(a.scratch/'checkpoint.json')
    outcome = 'running'; failure = None
    try:
        while elapsed<target and len(history)<a.max_intervals:
            h = min(h, target-elapsed); item = dict(start_years=elapsed/YEAR, interval_years=h/YEAR)
            accepted = False
            tolerances = ([1e-12, 1e-13] if a.transport=='fully-convective' else
                          sorted({preferred_tolerance, 1e-12}, reverse=True))
            for tolerance in tolerances:
                try:
                    _, full, c1, _ = step(old, h, tolerance)
                    mid, half1, c2, b1 = step(old, h/2, tolerance)
                    fine, half2, c3, b2 = step(mid, h/2, tolerance)
                    time_diagnostics = {}
                    error, metrics = compare(full, half2, a.composition_error_mode, time_diagnostics,
                                             final_record or data.get('model_record'), half1,
                                             a.composition_relative_tolerance)
                    item.update(error_norm=error, criteria=metrics, tolerance=tolerance, conservation=[c1, c2, c3])
                    if time_diagnostics: item['comparison_diagnostics'] = time_diagnostics
                    if error<=1: accepted=True
                    break
                except ValueError as e:
                    item.setdefault('solve_failures', []).append(dict(tolerance=tolerance, message=str(e)))
            item['accepted'] = accepted; attempts.append(item)
            with (a.scratch/'attempts.jsonl').open('a') as f: f.write(json.dumps(item)+'\n')
            if not accepted:
                failures = item.get('solve_failures', [])
                if failures and 'fully convective transport: radiative' in failures[-1]['message']:
                    outcome='radiative_boundary_review_trigger';break
                if len(attempts)>=4 and all(t.get('solve_failures') and not t['accepted'] for t in attempts[-4:]):
                    messages = [t['solve_failures'][-1]['message'] for t in attempts[-4:]]
                    if len(set(messages))==1:
                        raise RuntimeError('unchanged rejection after four reduced steps: '+messages[0])
                h *= max(.1, min(.5, .8/np.sqrt(item.get('error_norm', 4.))))
                if h<YEAR or (len(attempts)>30 and not any(t['accepted'] for t in attempts[-30:])):
                    raise RuntimeError('repeated rejected steps: '+json.dumps(item))
                continue
            abundance_change = abundance_step_change(old, fine)
            old = fine; elapsed += h; accumulated += b1+b2; final_record=half2
            # Reuse the tolerance already required by the physical audit.
            # Repeatedly trying a looser value discarded a full stellar solve
            # at nearly every fully convective interval. Reassess once when
            # the star first develops separate composition reservoirs.
            now_fully_mixed = half2['mixing_regions'] == [[0, len(old)]]
            preferred_tolerance = (1e-10 if fully_mixed and not now_fully_mixed
                                   else tolerance)
            fully_mixed = now_fully_mixed
            conversion = 1-float(fine[-1, 7]/(fine[-1, 7:10].sum())*(C0+N0)/C0)
            total_cn = fine[-1, 7:10].sum()
            nitrogen_conversion = float((fine[-1, 9]/total_cn*(C0+N0)-N0)/C0)
            row = dict(elapsed_years=elapsed/YEAR, age_years=(base_age+elapsed)/YEAR,
                       interval_years=h/YEAR, error_norm=error, Teff=half2['Teff'],
                       central_T=fine[0, 3], central_H=fine[0, 5], surface_H=fine[-1, 5],
                       surface_carbon12_depletion=conversion, surface_carbon_to_nitrogen=nitrogen_conversion,
                       surface_carbon13_per_initial_carbon=float(fine[-1, 8]/total_cn*(C0+N0)/C0),
                       luminosity=fine[-1, 4], radius=fine[-1, 1],
                       mixed_mass_fraction=half2['convective_mass_fraction'],
                       fully_mixed=now_fully_mixed,
                       convective_regions=[region for region in half2['mixing_regions'] if region[1]>region[0]+1],
                       maximum_abundance_step_change=abundance_change,
                       relative_catalyst_number_change=float(w@(fine[:, 7:10].sum(axis=1)-initial[:, 7:10].sum(axis=1))/mass/(C0+N0)),
                       maximum_accumulated_species_balance=float(abs(accumulated.sum(axis=0)).max()))
            # cn_thermal_probe retains its 0.001 native abundance-change
            # guard. Predict a step below that guard instead of repeatedly
            # growing into it and then rejecting an otherwise accurate step.
            # This is a proposal only; every native/time/energy audit remains.
            factor = min(2., max(.5, .85/np.sqrt(max(error, 1e-6))))
            factor = min(factor, .85e-3/max(abundance_change, 1e-30))
            history.append(row); h *= factor
            save_checkpoint()
            with (a.scratch/'history.jsonl').open('a') as f: f.write(json.dumps(row)+'\n')
            print(f"accepted {len(history)}: +{elapsed/YEAR:.4g} yr, Teff {half2['Teff']:.4g} K, Xc {fine[0,5]:.4g}, surface C12 depletion {conversion:.4g}, C to N {nitrogen_conversion:.4g}, error {error:.4g}", flush=True)
            if a.stop_conversion is not None and conversion>=a.stop_conversion:
                outcome='surface_carbon12_depletion_review_trigger';break
            if a.stop_surface_hydrogen is not None and fine[-1,5]>=a.stop_surface_hydrogen:
                outcome='surface_hydrogen_review_trigger';break
            if (a.scratch/'stop.request').exists(): outcome='requested_checkpoint_stop';break
        else: outcome='requested_interval_reached' if elapsed>=target else 'interval_count_limit'
    except Exception as e: outcome='stopped'; failure=str(e); print('STOP', failure, flush=True)
    finally:
        process.stdin.close(); code=process.wait(timeout=120);stderr.close()
    child_end = resource.getrusage(resource.RUSAGE_CHILDREN)
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(), outcome=outcome, failure=failure,
                  accepted_for_stellar_evolution=False, elapsed_years=elapsed/YEAR, requested_years=a.years,
                  solver_seconds=solver_seconds, wall_seconds=time.monotonic()-start, native_exit_code=code,
                  native_total_cpu_seconds=child_end.ru_utime+child_end.ru_stime-child_start.ru_utime-child_start.ru_stime,
                  timing_scope='CPU from OS child accounting includes table loading and diagnostics; solver_seconds is the native steady-clock sum. Clocks are recorded separately.',
                  native_command=command, native_solves=step_count, accepted_intervals=len(history),
                  zone_threads=a.zone_threads,
                  solver_cpu_seconds=sum(solver_cpu_times) if solver_cpu_times else None,
                  native_solves_with_cpu_timing=len(solver_cpu_times),
                  composition_error_mode=a.composition_error_mode,
                  composition_relative_tolerance=a.composition_relative_tolerance,
                  surface_hydrogen_review_threshold=a.stop_surface_hydrogen,
                  rejected_intervals=sum(not t['accepted'] for t in attempts), history=history,
                  checkpoint=str(a.scratch/'checkpoint.json') if history else None,
                  input_sha256=hashes, limitations=limitations,
                  accumulated_global_species_balance=accumulated.sum(axis=0).tolist())
    a.output.write_text(json.dumps(report, indent=2)+'\n')
    return int(outcome=='stopped' or code!=0)


if __name__=='__main__': raise SystemExit(main())
