#!/usr/bin/env python3
"""Check completed stellar runs and export the working paper's numerical data."""
import argparse
import csv
import hashlib
import json
import math
import subprocess
from pathlib import Path

from summarize_remnant_endpoints import summarize

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / 'docs/reports/2026-09-10'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path, columns, rows):
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(columns)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mixing-probe', type=Path, required=True)
    args = parser.parse_args()
    records = []
    for points, suffix in [(512, 'v2'), (1024, 'continuation-v1')]:
        path = ROOT / f'out/evolution-cold-remnant-x015-transition-{points}-3560gyr-{suffix}.json'
        receipt_path = path.with_suffix('.receipt.json')
        data = json.loads(path.read_text())
        receipt = json.loads(receipt_path.read_text())
        assert receipt['returncode'] == 0 and data['converged']
        assert digest(path) == receipt['output_sha256']
        for key in ['data_changed_during_run', 'source_changed_during_run', 'restart_changed_during_run']:
            assert not receipt[key], key
        assert data['points'] == points and data['mass_Msun'] == .1
        assert data['history'][-1][0] == 3.56e12
        assert len(data['profile']) == points
        for row in data['history']:
            for key, value in zip(data['columns'], row, strict=True):
                if value is None:
                    assert key == 'last_halfstep_gravothermal_Lsun'
                else:
                    assert math.isfinite(value)
        endpoint = dict(zip(data['columns'], data['history'][-1], strict=True))
        assert endpoint['central_X'] == data['profile'][0][5]
        assert endpoint['central_T_K'] == data['profile'][0][3]
        assert endpoint['surface_X'] == data['profile'][-1][5]
        for row in data['profile']:
            assert all(math.isfinite(v) for v in row)
            assert abs(sum(row[5:8]) + .02 - 1) < 1e-12
        profile_text = ''.join(' '.join(format(v, '.17g') for v in row)+'\n' for row in data['profile'])
        mixing = subprocess.run([str(args.mixing_probe),
                                 'data/eos/exhaustion_refined_v2/freeeos300_gs98_z020.dat',
                                 'data/opacity/hydrogen_poor_refined_v4'],
                                cwd=ROOT, input=profile_text, text=True, capture_output=True, check=True)
        regions = []
        for region in map(json.loads, mixing.stdout.splitlines()):
            region['inner_radius_fraction'] = data['profile'][region['begin']][1]/data['profile'][-1][1]
            region['outer_radius_fraction'] = data['profile'][region['end_exclusive']-1][1]/data['profile'][-1][1]
            if regions and not regions[-1]['convective'] and not region['convective']:
                previous = regions[-1]
                previous['mass_fraction'] += region['mass_fraction']
                for key in ['end_exclusive', 'outer_enclosed_fraction', 'outer_radius_fraction']:
                    previous[key] = region[key]
            else:
                regions.append(region)
        assert abs(sum(r['mass_fraction'] for r in regions if r['convective'])-endpoint['convective_mass_fraction']) < 1e-12
        history_name = 'evolution_history.csv' if points == 512 else 'evolution_1024_continuation.csv'
        write_csv(PAPER / history_name, data['columns'], data['history'])
        profile_name = f'evolution_profile_{points}.csv'
        write_csv(PAPER / profile_name, data['profile_columns'], data['profile'])
        records.append({
            'points': points, 'input': str(path.relative_to(ROOT)),
            'input_sha256': digest(path), 'receipt_sha256': digest(receipt_path),
            'executable_sha256': receipt['executable_sha256'],
            'atmosphere_sha256': receipt['data_sha256'][str(Path('/tmp/ember-nongrey-exhaustion-warm-x015-v1.dat').resolve())],
            'history_csv': history_name, 'history_csv_sha256': digest(PAPER / history_name),
            'profile_csv': profile_name, 'profile_csv_sha256': digest(PAPER / profile_name),
            'restart': data.get('restart'), 'rejected_steps_lifetime': data['rejected_steps'],
            'segment_cpu_seconds': receipt['child_user_seconds'] + receipt['child_system_seconds'],
            'segment_awake_elapsed_seconds': receipt['awake_elapsed_seconds'],
            'endpoint': endpoint, 'milestones': summarize(data),
            'regions': regions,
            'step_error_tolerances': data['step_error_tolerances'],
            'largest_absolute_discrete_luminosity_balance': max(abs(r[8]) for r in data['history']),
            'helium3_maximum_in_segment': dict(zip(data['columns'], max(data['history'], key=lambda r:r[6]), strict=True)),
        })
    assert records[0]['atmosphere_sha256'] == records[1]['atmosphere_sha256']
    assert records[0]['step_error_tolerances'] == records[1]['step_error_tolerances']
    low, high = [r['endpoint'] for r in records]
    quantities = ['R_Rsun', 'L_Lsun', 'Teff_K', 'central_X', 'central_Y3',
                  'central_T_K', 'central_rho', 'convective_mass_fraction', 'hydrogen_mass_Msun']
    result = {
        'scope': 'Completed runs at the same age and with the same physical inputs and time tolerances. Different executable bytes reflect checkpoint-counter handling, not a physical change. Two meshes do not establish a limiting result; tighter time-step controls remain needed.',
        'runs': records,
        'mixing_probe_sha256': digest(args.mixing_probe),
        'mixing_probe_source_sha256': digest(ROOT / 'scripts/evolution_mixing_probe.cpp'),
        'fractional_difference_1024_over_512_minus_one': {k:high[k]/low[k]-1 for k in quantities},
        'absolute_difference_1024_minus_512': {k:high[k]-low[k] for k in quantities},
    }
    output = ROOT / 'docs/results/evolution_transition_3560gyr_v1.json'
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result['fractional_difference_1024_over_512_minus_one'], indent=2))


if __name__ == '__main__':
    main()
