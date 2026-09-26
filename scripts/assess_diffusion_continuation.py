#!/usr/bin/env python3
"""Audit a continuation or the accepted prefix of an interrupted run."""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('report', type=Path)
    ap.add_argument('output', type=Path)
    ap.add_argument('--completed-prefix', action='store_true',
                    help='Verify saved accepted intervals without accepting the failed or unfinished tail.')
    args = ap.parse_args()
    assert not args.output.exists()
    identities = {}

    def read(p):
        p = Path(p).resolve(); data = p.read_bytes()
        identities[str(p)] = hashlib.sha256(data).hexdigest()
        return data

    report = json.loads(read(args.report))
    if not args.completed_prefix:
        assert report['target_reached'] and not report['failures']
    assert report['new_accepted_intervals'] > 0
    for p, h in report['artifacts_sha256'].items():
        assert hashlib.sha256(read(p)).hexdigest() == h, p
    paths = {Path(p).name: Path(p) for p in report['artifacts_sha256']}
    inputs = json.loads(read(paths['inputs.json']))
    records = [json.loads(s) for s in read(paths['responses.jsonl']).splitlines()
               if json.loads(s).get('kind') == 'step']
    journal = [json.loads(s) for s in read(paths['steps.jsonl']).splitlines()]
    iterator = iter(records); cache = {}; paired = []

    def key(s):
        return tuple(s[k] for k in ('initial_age_seconds','dt_seconds','seconds','converged','Teff'))

    for j in journal:
        if j['reused']:
            s = cache[key(j)]
        else:
            s = next(iterator); assert key(s) == key(j); cache[key(s)] = s
        paired.append((s,j))
    unmatched = list(iterator)
    if not args.completed_prefix:
        assert not unmatched
    accepted = {a['attempt'] for a in report['attempts'] if a.get('accepted_numerically',False)}
    verified = [(s,j) for s,j in paired if j['attempt'] in accepted and j.get('verification_passed',True)]
    selected = [s for s,j in verified if j['part'] != 'full']
    assert len(selected) == 2*report['new_accepted_intervals']
    initial_checkpoint = json.loads(read(report['resumed_checkpoint']))
    initial = initial_checkpoint.get('model_record',initial_checkpoint)
    a = np.asarray(initial['model']); m = a[:,0]
    w = np.zeros(512); w[0] = m[0]
    w[:-1] += .5*np.diff(m); w[1:] += .5*np.diff(m)
    assert a.shape == (512,7) and np.all(np.diff(m)>0)
    assert abs(w.sum()/m[-1]-1)<1e-14
    maximum_checks = {}

    def check_step(old, s, j):
        b = np.asarray(s['model']); log = np.asarray(s['model_log'])
        assert s['converged'] and s['input_preserved']
        assert b.shape == log.shape == (512,7) and np.isfinite(b).all() and np.isfinite(log).all()
        assert np.array_equal(b[:,0],m) and np.min(b[:,1:4])>0
        assert np.min(b[:,5:7])>=0 and np.max(b[:,5:7].sum(axis=1))<=.98
        assert s['initial_age_seconds']==old['age_seconds']
        assert s['dt_seconds']>0 and s['age_seconds']==old['age_seconds']+s['dt_seconds']
        assert s.get('cn_carbon',-1)==initial.get('cn_carbon',-1)
        assert s.get('parcel_radiation',False)==initial.get('parcel_radiation',False)
        rates = np.zeros((513,2)); rates[1:-1] = s['total_species_rates']
        storage = w[:,None]*(b[:,5:7]-np.asarray(old['model'])[:,5:7]
                            -s['dt_seconds']*np.asarray(s['nuclear_dXdt']))/m[-1]
        local = storage+s['dt_seconds']/m[-1]*np.diff(rates,axis=0)
        abundance_tolerance=inputs['coupling_tolerances']['abundance']
        for _ in range(j.get('verification_retry',0)):
            abundance_tolerance=max(1e-12,.1*abundance_tolerance)
        checks = dict(local_species=(float(abs(local).max()),1e-13),
                      global_species=(float(abs(storage.sum(axis=0)).max()),1e-13),
                      energy=(abs(s['luminosity_balance']),2e-8),
                      nuclear_mass=(abs(s['nuclear_mass_balance']),2e-6),
                      abundance_iteration=(s['abundance_residual'],abundance_tolerance),
                      heat_iteration=(s['material_heat_residual'],inputs['coupling_tolerances']['heat']))
        for name,(value,limit) in checks.items():
            assert np.isfinite(value) and value<=limit,(j['attempt'],j['part'],name,value,limit)
            maximum_checks[name]=max(maximum_checks.get(name,0.),value)
        for lo,hi in s['mixing_regions']:
            if hi>lo+1:assert np.max(abs(b[lo:hi,5:7]-b[lo,5:7]))==0

    def check_time(coarse,fine):
        tol=inputs['numerical_tolerances']
        a,b=np.asarray(coarse['model_log']),np.asarray(fine['model_log'])
        xa=np.column_stack((a[:,5:7],.98-a[:,5:7].sum(axis=1)))
        xb=np.column_stack((b[:,5:7],.98-b[:,5:7].sum(axis=1)))
        delta=xa-xb
        scale=tol['absolute_abundance']+tol['relative_abundance']*np.maximum(abs(xa),abs(xb))
        terms=dict(log_structure=float(abs(a[:,1:4]-b[:,1:4]).max()/tol['log_structure']),
                   scaled_abundance=float((abs(delta)/scale).max()),
                   relative_surface_luminosity=float(abs(a[-1,4]/b[-1,4]-1)/tol['relative_surface_luminosity']),
                   integrated_hydrogen=float(w@abs(a[:,5]-b[:,5])/(w@b[:,5])/tol['relative_integrated_absolute_hydrogen']))
        regions=[np.asarray([r for r in s['mixing_regions'] if r[1]>r[0]+1],dtype=int).reshape(-1,2)
                 for s in (coarse,fine)]
        cells=(float(abs(regions[0]-regions[1]).max()) if regions[0].size else 0.) \
            if regions[0].shape==regions[1].shape else float('inf')
        terms['convective_boundary_cells']=cells/tol['convective_boundary_cells']
        assert all(np.isfinite(v) and v<=1 for v in terms.values()),terms
        return terms

    old=initial; time_checks=[]
    for attempt in report['attempts']:
        if attempt['attempt'] not in accepted:continue
        parts=[(s,j) for s,j in verified if j['attempt']==attempt['attempt']]
        assert len(parts)==3 and {j['part'] for s,j in parts}=={'full','first_half','second_half'}
        by_part={j['part']:(s,j) for s,j in parts}
        full,jfull=by_part['full']; first,jfirst=by_part['first_half']; second,jsecond=by_part['second_half']
        assert full['dt_seconds']==2*first['dt_seconds']==2*second['dt_seconds']
        assert abs(full['dt_seconds']/(365.25*86400)-attempt['interval_years'])<1e-6
        check_step(old,full,jfull); check_step(old,first,jfirst); check_step(first,second,jsecond)
        time_checks.append(check_time(full,second))
        old=second
    assert all(s['initial_age_seconds']>=old['age_seconds'] for s in unmatched)
    source = np.zeros((512,2)); flux = np.zeros((513,2)); elapsed = 0.
    old = initial
    for s in selected:
        assert s['converged'] and s['initial_age_seconds'] == old['age_seconds']
        assert s.get('cn_carbon',-1)==initial.get('cn_carbon',-1)
        dt = s['dt_seconds']; elapsed += dt
        source += dt*np.asarray(s['nuclear_dXdt'])
        flux[1:-1] += dt*np.asarray(s['total_species_rates'])
        old = s
    b = np.asarray(old['model'])
    storage = w[:,None]*(b[:,5:7]-a[:,5:7]-source)/m[-1]
    local = storage + np.diff(flux,axis=0)/m[-1]
    global_balance = storage.sum(axis=0)
    # Accumulation is checked against the sum of the unchanged per-substep
    # balance allowances, as well as reported without normalization by it.
    allowance = len(selected)*1e-13
    assert max(float(abs(local).max()),float(abs(global_balance).max())) <= allowance
    assert abs(elapsed/(365.25*86400) - report['elapsed_years'] + report['initial_elapsed_years']) < 1e-6
    cp = json.loads(read(paths['checkpoint.json']))
    assert cp['model_record'] == old
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed',
                  scope='completed_prefix' if args.completed_prefix else 'completed_run',
                  source_run_report=str(args.report.resolve()),resumable_checkpoint=str(paths['checkpoint.json'].resolve()),
                  source_run_outcome=report['outcome'],source_run_failures=report['failures'],
                  unpaired_tail_native_replies=len(unmatched),
                  excluded_attempts=[a['attempt'] for a in report['attempts'] if a['attempt'] not in accepted],
                  verification_rejections=report.get('verification_rejections',0),
                  maximum_per_step_checks=maximum_checks,
                  maximum_time_error_terms={k:max(t[k] for t in time_checks) for k in time_checks[0]},
                  total_conditional_elapsed_years=report['elapsed_years'],new_elapsed_years=elapsed/(365.25*86400),
                  accepted_intervals=report['new_accepted_intervals'],accepted_substeps=len(selected),
                  native_solves=len(records),reused_trials=report['reused_trials'],
                  equation_rejections=report['equation_rejections'],numerical_rejections=report['numerical_rejections'],
                  complete_trial_solver_seconds=sum(s['seconds'] for s in records),wall_seconds=report['elapsed_run_seconds'],
                  accepted_interval_range_years=[min(h['interval_years'] for h in report['accepted_history']),
                                                 max(h['interval_years'] for h in report['accepted_history'])],
                  accumulated_local_species_balance=float(abs(local).max()),accumulated_global_species_balance=global_balance.tolist(),
                  sum_per_substep_balance_allowance=allowance,
                  maximum_relative_energy_imbalance=max(abs(s['luminosity_balance']) for s in selected),
                  endpoint=dict(Teff=old['Teff'],central_hydrogen=float(b[0,5]),surface_hydrogen=float(b[-1,5]),
                                central_temperature=float(b[0,3]),central_density=float(b[0,2]),
                                total_hydrogen_mass_g=float(w@b[:,5]),radius_cm=float(b[-1,1]),luminosity_erg_s=float(b[-1,4]),
                                convective_mass_fraction=old['convective_mass_fraction']),
                  total_hydrogen_change_g=float(w@(b[:,5]-a[:,5])),
                  accepted_for_stellar_evolution=False,selected_star_unchanged=True,
                  input_sha256=identities,
                  cn_carbon=initial.get('cn_carbon',-1),
                  limitations=['Local time control and species/energy balances pass; no independent long-interval endpoint refinement has yet been run.',
                               'The envelope approximation has finite-interval sensitivity checks.',
                               'The closed CN option does not restore historical carbon-conversion fuel or CN burning.'])
    read(__file__)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'input_sha256'},indent=2))


if __name__ == '__main__':
    main()
