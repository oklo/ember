#!/usr/bin/env python3
"""Check the refined pair table on retained and newly withheld stellar states."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import time

import numpy as np

from audit_electron_pair_interpolation import response, relative_operator_error
from electron_pair_table import ElectronPairTable


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    assert not args.output.exists()
    root = Path(__file__).resolve().parents[1]
    inputs = {}
    def read(path):
        path = Path(path)
        raw = path.read_bytes()
        inputs[str(path)] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)
    report_path = root/'docs/results/electron_pair_table_sources_v2.json'
    data = read(report_path)
    retained = read(root/'docs/results/electron_electron_energy_refinement_v1.json')
    material = read(root/'docs/results/diffusion_material_regime_3890gyr_v1.json')
    heat = read(root/'docs/results/stellar_electron_heat_3890gyr_v1.json')
    prior = read(root/'docs/results/electron_pair_interpolation_v1.json')
    source_pairs = {(s['eta'], s['b_thermal']) for s in
        (json.loads(Path(e['path']).read_text()) for e in data['completed_sources'])}
    cases = []
    for zone in (0, 300, 396):
        cases.append(dict(zone=zone, screening='total_effective_static_screening',
            independent_of_refinement=False, source=retained['calculations'][f'zone_{zone}']))
    for entry in data['independent_controls']:
        source = read(entry['path'])
        assert inputs[entry['path']] == entry['sha256']
        assert (source['eta'], source['b_thermal']) not in source_pairs
        cases.append(dict(zone=entry['zone'], screening=entry['screening'],
                          independent_of_refinement=True, source=source))
    metal_path = root/'include/ember/gs98_mixture.hpp'
    inputs[str(metal_path)] = hashlib.sha256(metal_path.read_bytes()).hexdigest()
    metals = [tuple(map(float, s.split(','))) for s in
              re.findall(r'^\s*\{([^{}]+)\}, //', metal_path.read_text(), re.M)]
    masses = np.array([1., 3., 4.]+[m[1] for m in metals])
    charges = np.array([1., 2., 2.]+[m[0] for m in metals])
    for name in ['electron_pair_table.py', 'audit_electron_pair_interpolation.py',
                 'audit_electron_pair_interpolation_v2.py', 'electron_mixture_collision.py']:
        path = root/'scripts'/name
        inputs[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    table = ElectronPairTable(report_path)
    grid_error = max(relative_operator_error(table.matrices(s['eta'], s['b_thermal'])[1:, 1:],
                     np.asarray(s['collision_matrix'])[1:, 1:]) for s in table.sources)
    controls = []
    for case in cases:
        zone, source = case['zone'], case['source']
        mr = material['records'][zone]
        hr = heat['records'][zone]
        sample = next(c for c in hr['cases'] if c['screening'] == case['screening'])
        assert mr['zone'] == hr['zone'] == zone
        assert source['eta'] == hr['eta_nonrelativistic']
        assert source['b_thermal'] == sample['b_thermal']
        ztotal = 1-mr['X']-mr['Y3']-mr['Y4']
        fractions = np.array([mr['X'], mr['Y3'], mr['Y4']]+[ztotal*m[3] for m in metals])
        populations = fractions/masses
        actual_ratio = float(populations@charges/(populations@(charges*charges)))
        ei, ee = np.asarray(source['electron_ion_matrix']), np.asarray(source['collision_matrix'])
        trial = table.matrices(source['eta'], source['b_thermal'])
        for ratio in (.25, .5, 1., actual_ratio):
            reference = response(ei, ee, ratio)
            interpolated = response(ei, trial, ratio)
            eight = response(ei[:8, :8], ee[:8, :8], ratio)
            controls.append(dict(zone=zone, screening=case['screening'],
                independent_of_refinement=case['independent_of_refinement'],
                eta=source['eta'], b_thermal=source['b_thermal'], prefactor_ratio=ratio,
                response_relative_operator_error=relative_operator_error(interpolated, reference),
                eight_to_ten_response_relative_error=relative_operator_error(eight, reference),
                heat_response_relative_difference=float(
                    (interpolated[1, 1]-interpolated[0, 1]**2/interpolated[0, 0])/
                    (reference[1, 1]-reference[0, 1]**2/reference[0, 0])-1)))
    etas = [r['eta_nonrelativistic'] for r in heat['records'] for c in r['cases']]
    bs = [c['b_thermal'] for r in heat['records'] for c in r['cases']]
    start = time.perf_counter()
    matrices = table.matrices(etas, bs)
    elapsed = time.perf_counter()-start
    np.linalg.cholesky(matrices[:, 1:, 1:])
    assert np.all(matrices[:, 0] == 0) and np.all(matrices[:, :, 0] == 0)
    worst = max(c['response_relative_operator_error'] for c in controls)
    new_worst = max(c['response_relative_operator_error'] for c in controls if c['independent_of_refinement'])
    mode_worst = max(c['eight_to_ten_response_relative_error'] for c in controls)
    passed = worst <= .005 and mode_worst <= .005 and grid_error < 2e-10
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(),
        outcome='passed_independent_controls' if passed else 'needs_refinement',
        accepted_for_stellar_evolution=False, interpolation_degree=3,
        grid_sources=len(data['completed_sources']), reused_sources=data['reused_source_count'],
        independent_stellar_states=len(data['independent_controls']), controls=controls,
        source_reproduction_error=grid_error, maximum_response_relative_error=worst,
        maximum_new_control_response_relative_error=new_worst,
        maximum_eight_to_ten_response_relative_error=mode_worst,
        response_tolerance=.005, stellar_queries=len(etas),
        stellar_query_elapsed_seconds=elapsed, every_stellar_heat_block_positive=True,
        exact_momentum_null_mode=True, input_sha256=inputs, new_pair_integrals=0,
        previous_domain_rejections_reused=prior['results'][1]['outside_domain_rejections'],
        limitations=['Eight new withheld source states sample both screening prescriptions; they are not exhaustive error bounds.',
                     'Interpolation accuracy is relative to the prescribed nonrelativistic Born/Pauli collision model.',
                     'No conservative stellar abundance, energy or structure evolution is performed.'])
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({k: report[k] for k in ['outcome', 'maximum_response_relative_error',
        'maximum_new_control_response_relative_error', 'maximum_eight_to_ten_response_relative_error',
        'stellar_query_elapsed_seconds']}, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
