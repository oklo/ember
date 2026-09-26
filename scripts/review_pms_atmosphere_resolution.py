"""Measure atmosphere resolution error and cost at identical physical inputs.

This is a local comparison with a reviewed fine source. It does not select a
production resolution or establish an error bound throughout the atmosphere grid.
"""
import datetime
import json
from pathlib import Path
import sys

from generate_nongrey_grid import (
    completed_initial_structure, input_fingerprint, resample_initial_structure,
    sequence, temperatures,
)
from import_nongrey_grid import source_inputs, source_state
from prepare_nongrey_sources import digest


def review(work):
    work = Path(work).resolve()
    plan = json.loads((work / 'plan.json').read_text())
    assert digest(work / 'run.py') == plan['script_sha256']
    reference_path = Path(plan['reference_review'])
    assert digest(reference_path) == plan['reference_review_sha256']
    reference = json.loads(reference_path.read_text())
    assert reference['all_checks_pass']
    for path, expected in reference['input_sha256'].items():
        assert digest(path) == expected, path
    fine = max(reference['records'], key=lambda r: r['state']['optical_depth_range'][1])
    donor = Path(fine['work']).resolve()
    old_plan = json.loads((donor.parent / 'plan.json').read_text())
    old_case = next(c for c in old_plan['cases'] if c['name'] == donor.name)
    old_spec = old_case['specification']
    assert plan['prepared'] == old_plan['prepared']
    prepared = plan['prepared']
    assert digest(prepared['tlusty']) == prepared['executables']['tlusty']
    initial = completed_initial_structure(donor, old_spec, fine['Teff_K'], fine['log_g'], 300)
    hashes = {str(p): digest(p) for p in [work / 'plan.json', work / 'run.py',
                                         reference_path, Path(__file__)]}
    records, missing = [], []
    for case in plan['cases']:
        p = work / case['name']
        if not (p / 'validated.json').exists():
            missing.append(case['name'])
            continue
        spec = case['specification']
        T, g = case['teff_K'], case['log_g']
        assert (T, g) == (fine['Teff_K'], fine['log_g'])
        assert (case['XH'], case['X3']) == (.7, 0.)
        assert Path(case['donor']).resolve() == donor
        assert digest(donor / 'completed.json') == case['donor_completed_sha256']
        assert digest(donor / 'validated.json') == case['donor_validation_sha256']
        for key in set(spec) | set(old_spec):
            if key not in ('depths', 'atmosphere_frequencies', 'source'):
                assert spec[key] == old_spec[key], key
        assert case['opacity_sha256'] == old_plan['opacity_sha256']
        assert digest(case['opacity']) == case['opacity_sha256']
        reconstructed = resample_initial_structure(initial, spec['depths'], allow_coarsen=True)
        assert (p / 'fort.8').read_text() == reconstructed
        assert digest(p / 'fort.8') == case['initial_sha256']
        completed = json.loads((p / 'completed.json').read_text())
        assert input_fingerprint(prepared['executables']['tlusty'], p) == completed['input_sha256']
        for filename, expected in completed['outputs'].items():
            assert digest(p / filename) == expected
        validation = json.loads((p / 'validated.json').read_text())
        assert validation['plan_sha256'] == digest(work / 'plan.json')
        log = (p / 'run.log').read_text()
        source_inputs({k: (p / f).read_text() for k, f in [
            ('atmosphere_input', 'fort.5'), ('parameters', 'tas'),
            ('element_masses', 'ember-masses.dat'), ('initial_structure', 'fort.8')]},
            spec, .7, 0., T, g, log)
        state = source_state(log, (p / 'fort.9').read_text(), T, g,
                             dict(temperature_K=temperatures(spec),
                                  density_g_cm3=sequence(spec['log_density'])), 100.)
        assert state == validation['diagnostics']
        # Depth coarsening must preserve the physical column endpoints.
        for before, after in zip(fine['state']['column_mass_range'], state['column_mass_range']):
            assert abs(after / before - 1) < 2e-6
        records.append(dict(case=case['name'], depths=spec['depths'],
                            atmosphere_frequencies=spec['atmosphere_frequencies'], state=state,
                            relative_changes={k: state[k] / fine['state'][k] - 1
                                              for k in ('T', 'Pgas', 'source_density')},
                            source_cpu_seconds=validation['source_cpu_seconds']))
        for name in ['validated.json', 'completed.json', 'run.log', 'fort.5', 'fort.7',
                     'fort.8', 'fort.9', 'tas', 'ember-masses.dat', 'opacity.sha256', 'physics.json']:
            hashes[str(p / name)] = digest(p / name)
        hashes[case['opacity']] = case['opacity_sha256']
    return dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                completed_source_checks_pass=bool(records), all_cases_complete=not missing,
                missing_cases=missing, records=records, reference=fine,
                selected_for_physical_PMS_track=False, input_sha256=hashes,
                limitations=[
                    'Only one Teff/gravity/composition point; not a global atmosphere error bound.',
                    'Depth count and transfer frequency count change together; their individual effects are not separated.',
                    'Candidate columns start from the converged fine profile, so measured CPU is warm-start cost.',
                    'Matching depth is tau=100 with the same column endpoints; candidate depth sensitivity is not independently tested.',
                    'No scientific acceptance threshold is chosen by this comparison.'])


if __name__ == '__main__':
    result = review(sys.argv[1])
    destination = Path(sys.argv[2])
    assert not destination.exists(), 'use a fresh review path'
    destination.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['all_cases_complete', 'missing_cases', 'records']}, indent=2))
