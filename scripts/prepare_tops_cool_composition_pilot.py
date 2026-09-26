#!/usr/bin/env python3
"""Prepare a bounded composition-only opacity comparison near the envelope.

Hold density and native temperature fixed while testing a finer hydrogen grid.
Use 993 photon boundaries from 1e-8 keV to 800 kT. All direct X=.1/.3/.7,
Z=.02 controls remain withheld from the candidate X/Z interpolation. This
pilot does not construct or accept a stellar opacity family.
"""
import json
from pathlib import Path
from datetime import datetime, timezone

from fetch_tops_composition import digest
from reduce_tops_group_factors import add_inputs, verify


def write(path, value):
    with Path(path).open('x') as stream:
        stream.write(json.dumps(value, indent=2, allow_nan=False) + '\n')


def main():
    root = Path(__file__).resolve().parents[1]
    reservations = root / 'docs/research/fable/coordination/reservations'
    for path in reservations.glob('*.json'):
        if json.loads(path.read_text()).get('status') not in ('released', 'complete', 'completed', 'cancelled', 'failed', 'superseded'):
            raise ValueError('inspect existing reservation before preparing sources: ' + str(path))
    inputs = {}
    for name in ('tops_cool_transport_regime_v1.json', 'tops_cool_full_spectra_v4.json'):
        path = root / 'docs/results' / name
        report = json.loads(path.read_text())
        add_inputs(inputs, report['input_sha256'])
        add_inputs(inputs, {str(path.resolve()): digest(path)})
    source_plan = root / 'data/opacity/sources/tops_native_uncut_family_batched_form_v4_specification.json'
    source = json.loads(source_plan.read_text())
    add_inputs(inputs, {str(source_plan.resolve()): digest(source_plan)})
    seed_by_composition = {}
    for job in source['requests']:
        key = job['X'], job['Z']
        if key not in seed_by_composition:
            seed = Path(job['work'])
            receipt = json.loads((seed / 'receipt.json').read_text())
            if (receipt['X'], receipt['Z']) != key:
                raise ValueError('source composition differs from its plan')
            for name, expected in (('source.txt', receipt['sha256']), ('request.json', receipt['request_sha256'])):
                if digest(seed / name) != expected:
                    raise ValueError('baseline source changed')
                add_inputs(inputs, {str((seed / name).resolve()): expected})
            add_inputs(inputs, {str((seed / 'receipt.json').resolve()): digest(seed / 'receipt.json')})
            seed_by_composition[key] = str(seed)
    temperatures = { .003: [.1797, .60509], .005: [.1797, .60509],
                     .008: [.60509, 2.0375], .0125: [2.0375, 6.8606] }
    candidates = [(x,z) for x in (.05,.15,.25,.4,.65,.75) for z in (.01,.03)]
    controls = [(x,.02) for x in (.1,.3,.7)]
    # Inspect existing requests at the chosen frequency spacing before creating
    # any new source. All planned states differ from those existing requests.
    existing = []
    for path in (root / 'data/opacity/sources').glob('tops*_specification.json'):
        plan = json.loads(path.read_text())
        if plan.get('photon_boundaries') != 993 or 'baseline' not in plan:
            continue
        receipt = json.loads((Path(plan['baseline']) / 'receipt.json').read_text())
        key = receipt['X'], receipt['Z']
        for job in plan['requests']:
            for t in job['temperatures_keV']:
                if key not in candidates + controls or t not in temperatures:
                    continue
                if (plan['photon_min_keV'], plan['photon_max_keV']) != (1e-8, 800*t):
                    continue
                for rho in job['densities_g_cm3']:
                    existing.append({'X':key[0], 'Z':key[1], 'temperature_keV':t,
                                     'density_atomic_g_cm3':rho, 'work':job['work'], 'plan':str(path)})
                    if any(abs(rho/r-1) < 5e-5 for r in temperatures[t]):
                        raise ValueError('potential reusable source: ' + str(job))
    verify(inputs)
    work = Path('/tmp/ember-tops-cool-composition-pilot-v1')
    monitor = Path('/tmp/ember-tops-cool-composition-pilot-run-v1')
    monitor.mkdir()
    tasks, plans, reductions = [], [], []
    for role, compositions in (('candidate', candidates), ('control', controls)):
        for x, z in compositions:
            seed = seed_by_composition[x,z]
            for t, rr in temperatures.items():
                label = f'x{round(x*1000):03d}-z{round(z*1000):03d}-t{round(t*1e6):06d}'
                plan_path = root / 'data/opacity/sources' / f'tops_cool_composition_pilot_{label}_v1_specification.json'
                report_path = root / 'docs/results' / f'tops_cool_composition_pilot_{label}_reduction_v1.json'
                directory = work / label
                job = {'temperatures_keV':[t], 'densities_g_cm3':rr,
                       'density_spacing':'log', 'work':str(directory)}
                plan = {'scope':__doc__, 'baseline':seed, 'role':role,
                        'photon_min_keV':1e-8, 'photon_max_keV':800*t,
                        'photon_boundaries':993, 'maximum_response_bytes':20000000,
                        'maximum_group_states':100000, 'requests':[job]}
                write(plan_path, plan)
                add_inputs(inputs, {str(plan_path.resolve()):digest(plan_path)})
                plans.append(str(plan_path.relative_to(root)))
                tasks.append({'kind':'groups', 'plan':str(plan_path.relative_to(root)),
                              'request_index':0, 'work':str(directory)})
                reductions.append({'role':role, 'X':x, 'Z':z, 'temperature_keV':t,
                                   'plan':str(plan_path.relative_to(root)),
                                   'report':str(report_path.relative_to(root))})
    for path in (Path(__file__), *(root/'scripts'/name for name in
                ('prepare_tops_form_lease_v3.py','fetch_tops_form_plan_v3.py',
                 'tops_groups.py','reduce_tops_group_controls.py'))):
        add_inputs(inputs, {str(path.resolve()):digest(path)})
    comparison_path = root / 'data/opacity/sources/tops_cool_composition_pilot_comparison_v1_specification.json'
    write(comparison_path, {'scope':__doc__, 'candidate_X':[.05,.15,.25,.4,.65,.75],
                           'candidate_Z':[.01,.03], 'control_X':[.1,.3,.7],
                           'control_Z':[.02], 'interpolation':'bilinear in log ratio over X and Z',
                           'relative_criterion':.005, 'reductions':reductions})
    add_inputs(inputs, {str(comparison_path.resolve()):digest(comparison_path)})
    tf, ef = monitor/'tasks.json', monitor/'existing_sources.json'
    write(tf,tasks); write(ef,existing)
    reservation = reservations/'primary_tops_cool_composition_pilot_v1.json'
    write(reservation, {'owner':'primary','task_ids':['E-OPACITY-COOL-COMPOSITION-PILOT'],
        'utc':datetime.now(timezone.utc).isoformat(),'status':'reserved','threads':4,
        'expected_memory_bytes':1200000000,'data_cap_bytes':100000000,'maximum_elapsed_seconds':900,
        'maximum_attempts_per_task':3,'stop_after_consecutive_failed_tasks':3,
        'scratch_paths':[str(work),str(monitor)],'work_directory':str(monitor),
        'task_file':str(tf),'task_file_sha256':digest(tf),
        'existing_sources_file':str(ef),'existing_sources_sha256':digest(ef),
        'input_sha256':inputs,'plans':plans})
    controller = Path('/tmp/ember-tops-cool-frequency-scaled-run-v1.py').read_text().replace(
        'primary_tops_cool_frequency_scaled_v1.json', reservation.name)
    with Path('/tmp/ember-tops-cool-composition-pilot-run-v1.py').open('x') as stream:
        stream.write(controller)
    print(json.dumps({'requests':len(tasks),'source_states':2*len(tasks),
                      'candidate_compositions':len(candidates),'control_compositions':len(controls),
                      'existing_matching_frequency_states_inspected':len(existing)}))


if __name__ == '__main__':
    main()
