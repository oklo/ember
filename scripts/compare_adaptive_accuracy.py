#!/usr/bin/env python3
"""Compare completed trajectories and their accumulated species balances."""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from audit_full_star_timestep import compare


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('output', type=Path)
    ap.add_argument('--candidate', type=int, nargs='+', required=True)
    args = ap.parse_args()
    assert not args.output.exists()
    identities = {}

    def read(p):
        p = Path(p)
        data = p.read_bytes()
        identities[str(p)] = hashlib.sha256(data).hexdigest()
        return data

    initial_path = Path('/tmp/ember-parcel-star-comparison-v1/checkpoint-radiation-1.json')
    initial = json.loads(read(initial_path))['model_record']
    model = np.asarray(initial['model'])
    m = model[:, 0]
    w = np.zeros(512); w[0] = m[0]
    w[:-1] += .5*np.diff(m); w[1:] += .5*np.diff(m)

    def collect(versions):
        steps, all_steps, intervals = [], [], []
        wall = 0.
        for v in versions:
            report = json.loads(read(f'docs/results/full_star_adaptive_diffusion_v{v}.json'))
            assert not report['failures']
            raw = Path(f'/tmp/ember-adaptive-diffusion-v{v}')
            for p, expected in report['artifacts_sha256'].items():
                assert hashlib.sha256(read(p)).hexdigest() == expected, p
            replies = [json.loads(s) for s in read(raw/'responses.jsonl').splitlines()]
            records = [s for s in replies if s.get('kind') == 'step']
            journal = [json.loads(s) for s in read(raw/'steps.jsonl').splitlines()]
            assert len(records) == len(journal) and not any(s['reused'] for s in journal)
            history = [json.loads(s) for s in read(raw/'history.jsonl').splitlines()][1:]
            ages = {h['age_seconds'] for h in history}
            accepted = {s['attempt'] for s in journal if s['part'] == 'second_half' and s['age_seconds'] in ages}
            assert len(accepted) == len(history)
            steps.extend(s for s,j in zip(records,journal) if j['attempt'] in accepted and j['part'] != 'full')
            all_steps.extend(records); intervals.extend(history)
            receipt = json.loads(read(f'/tmp/ember-adaptive-diffusion-run-v{v}/receipt.json'))
            wall += receipt['jobs'][0]['elapsed_seconds']
        source = np.zeros((512,2)); flux = np.zeros((513,2))
        old = initial
        duration = 0.
        for s in steps:
            assert s['converged'] and s['initial_age_seconds'] == old['age_seconds']
            duration += s['dt_seconds']
            source += s['dt_seconds']*np.asarray(s['nuclear_dXdt'])
            flux[1:-1] += s['dt_seconds']*np.asarray(s['total_species_rates'])
            old = s
        a = np.asarray(old['model'])
        storage = w[:,None]*(a[:,5:7]-model[:,5:7]-source)/m[-1]
        local = storage + np.diff(flux,axis=0)/m[-1]
        balance = dict(local=float(abs(local).max()),global_species=storage.sum(axis=0).tolist())
        assert max(balance['local'],max(abs(x) for x in balance['global_species'])) <= 1e-13
        assert abs(duration/(365.25*86400)-8.9e6) < 1e-7
        summary = dict(versions=versions, accepted_intervals=len(intervals), accepted_substeps=len(steps),
                       complete_trial_solves=len(all_steps), equation_failures=sum(not s['converged'] for s in all_steps),
                       all_trial_solver_seconds=sum(s['seconds'] for s in all_steps),
                       successful_trial_solver_seconds=sum(s['seconds'] for s in all_steps if s['converged']),
                       wall_seconds_including_load_and_interrupted_work=wall,
                       interval_range_years=[min(s['interval_years'] for s in intervals),max(s['interval_years'] for s in intervals)],
                       accumulated_species_balance=balance,
                       maximum_relative_energy_imbalance=max(abs(s['luminosity_balance']) for s in steps))
        return old, summary

    strict, baseline = collect([3])
    candidate, result = collect(args.candidate)
    a, b = np.asarray(strict['model']), np.asarray(candidate['model'])
    quantities = dict(central_hydrogen=[a[0,5],b[0,5]],surface_hydrogen=[a[-1,5],b[-1,5]],
                      total_hydrogen_mass=[w@a[:,5],w@b[:,5]],total_helium3_mass=[w@a[:,6],w@b[:,6]],
                      luminosity=[a[-1,4],b[-1,4]],radius=[a[-1,1],b[-1,1]],
                      effective_temperature=[strict['Teff'],candidate['Teff']],
                      mixed_mass_fraction=[strict['convective_mass_fraction'],candidate['convective_mass_fraction']])
    relative = {k:float(v[1]/v[0]-1) for k,v in quantities.items()}
    passed = all(abs(relative[k]) < 5e-4 for k in ('luminosity','radius','effective_temperature'))
    passed &= abs(relative['total_hydrogen_mass']) < 1e-3 and abs(relative['central_hydrogen']) < 2e-3
    passed &= strict['mixing_regions'] == candidate['mixing_regions']
    read(__file__)
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed' if passed else 'failed',
                  elapsed_comparison_years=8.9e6,total_conditional_elapsed_years=1e7,
                  baseline=baseline,candidate=result,endpoint_values=quantities,relative_endpoint_changes=relative,
                  full_structure_comparison=compare(strict,candidate),
                  speedup_all_solver_time=baseline['all_trial_solver_seconds']/result['all_trial_solver_seconds'],
                  speedup_successful_solver_time=baseline['successful_trial_solver_seconds']/result['successful_trial_solver_seconds'],
                  same_mixed_regions=strict['mixing_regions']==candidate['mixing_regions'],
                  input_sha256=identities,
                  limitations=['This is a short endpoint convergence comparison, not a lifetime-error bound.',
                               'Conditional opacity and envelope heat physics are identical on both paths.',
                               'Wall time includes table loads and stopped work; solver time includes complete native replies only.'])
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('input_sha256','endpoint_values')},indent=2))
    assert passed


if __name__ == '__main__':
    main()
