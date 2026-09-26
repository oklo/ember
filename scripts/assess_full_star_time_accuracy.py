#!/usr/bin/env python3
"""Apply the existing evolution driver's time-error tests to retained states."""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from audit_full_star_diffusion import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    assert not args.output.exists()
    root = Path(__file__).resolve().parents[1]
    report_path = root / 'docs/results/full_star_timestep_v4.json'
    raw_path = Path('/tmp/ember-full-star-timestep-v4/responses.jsonl')
    original_raw = Path('/tmp/ember-full-star-diffusion-v1/responses.jsonl')
    report = json.loads(report_path.read_text())
    assert report['outcome'] == 'passed_conditional_step_equations' and not report['failures']
    records = [json.loads(line) for line in raw_path.read_text().splitlines()]
    steps = [r for r in records if r.get('kind') == 'step']
    assert len(steps) == len(report['steps']) == 14
    states = {(c['subdivisions'], c['index']): r for c, r in zip(report['steps'], steps)}
    original = [json.loads(line) for line in original_raw.read_text().splitlines()]
    states[1, 0], = [r for r in original if r.get('kind') == 'step']
    comparisons = []
    for coarse in (1, 2, 4):
        a = np.asarray(states[coarse, 0]['model'])
        b = np.asarray(states[2 * coarse, 1]['model'])
        log_difference = abs(np.log(a[:, 1:4] / b[:, 1:4]))
        abundance_difference = a[:, 5:7] - b[:, 5:7]
        # With fixed metals, helium-4 compensates H1 and He3 exactly.
        absolute_abundance = max(float(np.max(abs(abundance_difference))),
                                 float(np.max(abs(abundance_difference.sum(axis=1)))))
        surface_luminosity = abs(a[-1, 4] / b[-1, 4] - 1)
        terms = dict(log_structure=float(log_difference.max() / 1e-5),
                     absolute_abundance=absolute_abundance / 1e-8,
                     relative_surface_luminosity=float(surface_luminosity / 1e-4))
        maximum = np.unravel_index(np.argmax(log_difference), log_difference.shape)
        comparisons.append(dict(interval_years=1e6 / coarse,
                                coarse_case=[coarse, 0], fine_case=[2 * coarse, 1],
                                error_terms=terms, error_norm=max(terms.values()),
                                meets_driver_time_accuracy=max(terms.values()) <= 1,
                                largest_structure_difference=dict(node=int(maximum[0]),
                                    variable=('lnr', 'lnrho', 'lnT')[maximum[1]],
                                    difference=float(log_difference[maximum]))))
    inputs = [Path(__file__), root / 'apps/evolve.cpp', report_path, raw_path,
              original_raw, root / 'docs/results/full_star_diffusion_v1.json']
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  outcome='completed_existing_driver_accuracy_assessment',
                  accepted_for_stellar_evolution=False, selected_star_unchanged=True,
                  native_steps_recomputed=0, comparisons=comparisons,
                  driver_tolerance_scale=1,
                  driver_tolerances=dict(log_structure=1e-5, absolute_abundance=1e-8,
                                         relative_surface_luminosity=1e-4),
                  next_initial_trial_years=125000,
                  next_required_comparison='Reuse the existing 125000-year coarse state and compute two consecutive 62500-year steps from the original profile.',
                  input_sha256={str(p): digest(p) for p in inputs},
                  limitations=['The earlier diagnostic targets differ from the production driver; they do not replace its tolerances.',
                               'A common initial interval is compared here, unlike full-million-year accumulation with different subdivisions.',
                               'The largest short-interval structural difference is at the atmospheric boundary; further initial refinement is needed.',
                               'The states use conditional cooler opacity and heat physics and are not selected trajectory points.',
                               'The retained one-million-year state stores physical variables; taking logarithms may add roundoff far below the measured differences.'])
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
