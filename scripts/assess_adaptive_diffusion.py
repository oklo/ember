#!/usr/bin/env python3
"""Reconstruct accepted substeps and audit the accumulated species balance."""
import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from audit_full_star_diffusion import digest
from audit_full_star_timestep import compare


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    assert not args.output.exists()
    root = Path(__file__).resolve().parents[1]
    paths = [root / 'docs/results/full_star_initial_refinement_v2.json',
             root / 'docs/results/full_star_adaptive_diffusion_v1.json',
             root / 'docs/results/full_star_timestep_v4.json']
    reports = [json.loads(p.read_text()) for p in paths]
    identities = {str(p): digest(p) for p in paths}
    raw_sets = []
    for r in reports:
        assert not r['failures']
        for p, expected in r['artifacts_sha256'].items():
            assert digest(p) == expected, p
            identities[p] = expected
        raw_path, = [Path(p) for p in r['artifacts_sha256'] if Path(p).name == 'responses.jsonl']
        raw_sets.append([json.loads(s) for s in raw_path.read_text().splitlines()
                         if json.loads(s).get('kind') == 'step'])
    journal_path = Path('/tmp/ember-adaptive-diffusion-v1/steps.jsonl')
    journal = [json.loads(s) for s in journal_path.read_text().splitlines()]
    assert len(journal) == len(raw_sets[1]) == 30
    assert not any(s['reused'] for s in journal)
    accepted_ids = {s['attempt'] for s in reports[1]['attempts'] if s['accepted_numerically']}
    selected = raw_sets[0] + [s for s, j in zip(raw_sets[1], journal)
                             if j['attempt'] in accepted_ids and j['part'] in ('first_half', 'second_half')]
    assert len(selected) == 22
    profile = root / 'docs/reports/2026-09-11/evolution_latest_profile.csv'
    with profile.open() as stream:
        initial = np.array([[float(r[k]) for k in ('mass_g', 'radius_cm', 'density_g_cm3',
                                                   'temperature_K', 'luminosity_erg_s', 'X', 'Y3')]
                            for r in csv.DictReader(stream)])
    mass = initial[:, 0]
    w = np.zeros(512); w[0] = mass[0]
    w[:-1] += .5 * np.diff(mass); w[1:] += .5 * np.diff(mass)
    integrated_source = np.zeros((512, 2)); integrated_flux = np.zeros((513, 2))
    duration = 0.; previous = initial; age = selected[0]['initial_age_seconds']
    for step in selected:
        assert step['converged'] and step['initial_age_seconds'] == age
        dt = step['dt_seconds']; duration += dt
        integrated_source += dt * np.asarray(step['nuclear_dXdt'])
        integrated_flux[1:-1] += dt * np.asarray(step['total_species_rates'])
        previous = np.asarray(step['model']); age = step['age_seconds']
    residual = (w[:, None] * (previous[:, 5:7] - initial[:, 5:7] - integrated_source)
                + np.diff(integrated_flux, axis=0)) / mass[-1]
    local = float(abs(residual).max())
    global_balance = (w[:, None] * (previous[:, 5:7] - initial[:, 5:7] - integrated_source)).sum(axis=0) / mass[-1]
    assert local <= 1e-13 and abs(global_balance).max() <= 1e-13
    assert abs(duration / (365.25 * 86400) - 1e6) < 1e-8
    uniform = [s for s, c in zip(raw_sets[2], reports[2]['steps']) if c['index'] == c['subdivisions'] - 1]
    endpoint_comparisons = [dict(uniform_subdivisions=c, metrics=compare(s, selected[-1]))
                            for c, s in zip((2, 4, 8), uniform)]
    assert len(endpoint_comparisons) == 3
    identities.update({str(profile): digest(profile), str(Path(__file__).resolve()): digest(Path(__file__))})
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(), outcome='verified_conditional_adaptive_path',
                  accepted_for_stellar_evolution=False, selected_star_unchanged=True,
                  accepted_implicit_substeps=len(selected), elapsed_years=duration / (365.25 * 86400),
                  accumulated_local_species_balance=local, accumulated_global_species_balance=global_balance.tolist(),
                  total_hydrogen_change_over_stellar_mass=float(np.dot(w, previous[:, 5] - initial[:, 5]) / mass[-1]),
                  uniform_endpoint_comparisons=endpoint_comparisons,
                  all_trial_solver_seconds=sum(s['seconds'] for s in raw_sets[1]),
                  accepted_interval_range_years=[min(s['interval_years'] for s in reports[1]['accepted_history']),
                                                 max(s['interval_years'] for s in reports[1]['accepted_history'])],
                  input_sha256=identities,
                  limitations=['Conditional physical assumptions remain unchanged and unaccepted.',
                               'Local step accuracy does not bound accumulated full-track error.'])
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'input_sha256'}, indent=2))


if __name__ == '__main__':
    main()
