"""Extract checked late models, preserving the retained time sequence."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import gzip
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parents[1]


def collect(last_run, output, last_checkpoint=None):
    pins = {}

    def read(path):
        path = ROOT / path
        data = path.read_bytes()
        pins[str(path)] = hashlib.sha256(data).hexdigest()
        return json.loads(data)

    first = read('docs/results/pulse_late_temperature_turn_v33_analysis.json')
    rows = first['retained_points'].copy()
    for row in rows:
        row['run'] = 33
    previous = read('out/pulse-front-evolution-review-v33/checkpoint.json')
    cpu_seconds = 0.
    years = 0.
    for version in range(34, last_run + 1):
        base = ROOT / f'out/pulse-front-evolution-v{version}'
        plan = read(base / 'plan.json')
        review = read(f'docs/results/pulse_front_evolution_primary_review_v{version}.json')
        run = read(f'docs/results/pulse_front_evolution_v{version}.json')
        assert review['maximum_error_norm'] <= 1
        assert review['accepted_intervals'] == len(run['accepted'])
        assert math.isclose(plan['prefix_years'], previous['total_pulse_years'], abs_tol=1e-7)
        old = previous['model_record']
        elapsed = 0.
        for entry in run['accepted']:
            assert entry['accepted'] and entry['error_norm'] <= 1
            queries = entry.get('retained_queries', entry['queries'][1:])
            for number in queries:
                path = base / 'queries' / f'{number:05d}.json.gz'
                data = path.read_bytes()
                digest = hashlib.sha256(data).hexdigest()
                assert review['input_sha256'][str(path)] == digest
                pins[str(path)] = digest
                raw = json.loads(gzip.decompress(data))
                reply = json.loads(raw['raw_reply'])
                tokens = raw['query'].split()
                count = int(tokens[5])
                assert count == len(old['model']) == 2115
                assert float(tokens[1]) == old['age_seconds']
                values = [float(x) for x in tokens[6:6 + count * 11]]
                assert values == [x for cell in old['model'] for x in cell]
                elapsed += float(tokens[2]) / 31557600.
                surface = reply['model'][-1]
                rows.append(dict(run=version, pulse_years=plan['prefix_years'] + elapsed,
                                 Teff_K=reply['Teff'], radius_cm=surface[1],
                                 photon_luminosity_erg_s=surface[4],
                                 nuclear_power_erg_s=reply['nuclear_luminosity'],
                                 logg=math.log10(6.6743e-8*surface[0]/surface[1]**2)))
                old = reply
        if version == last_run and last_checkpoint is not None:
            checkpoint = last_checkpoint
        else:
            checkpoint = ROOT / f'out/pulse-front-evolution-review-v{version}/checkpoint.json'
            if not checkpoint.is_file():
                # A corrected review keeps its failed predecessor intact. Its
                # checkpoint must still reproduce every retained model below.
                checkpoint = ROOT / f'out/pulse-front-evolution-review-v{version}-checked-ceiling/checkpoint.json'
                assert checkpoint.is_file(), f'missing reviewed checkpoint for run {version}'
        previous = read(checkpoint)
        assert old == previous['model_record']
        assert math.isclose(elapsed, review['elapsed_years'], abs_tol=1e-7)
        cpu_seconds += review['reported_completed_native_cpu_seconds']
        years += review['elapsed_years']
    assert all(b['pulse_years'] > a['pulse_years'] for a, b in zip(rows, rows[1:]))
    peak_index = max(range(len(rows)), key=lambda i: rows[i]['Teff_K'])
    assert all(b['Teff_K'] < a['Teff_K'] for a, b in zip(rows[peak_index:], rows[peak_index+1:]))
    end = rows[-1]
    maxima = read('docs/results/pulse_late_radius_turnaround_v1.json')
    facts = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                 first=rows[0], last=end, temperature_maximum=rows[peak_index],
                 retained_points=len(rows), consecutive_declining_half_steps=len(rows)-peak_index-1,
                 temperature_drop_K=rows[peak_index]['Teff_K']-end['Teff_K'],
                 decline_years=end['pulse_years']-rows[peak_index]['pulse_years'],
                 radius_decline=1-end['radius_cm']/maxima['radius_maximum_saved']['radius_cm'],
                 luminosity_decline=1-end['photon_luminosity_erg_s']/maxima['surface_luminosity_maximum_saved']['surface_luminosity'],
                 nuclear_fraction=end['nuclear_power_erg_s']/end['photon_luminosity_erg_s'],
                 appended_segment_years=years, appended_segment_cpu_minutes=cpu_seconds/60,
                 input_sha256=pins,
                 scope='Checked retained late models. No interpolated evolution, connection to the distinct onset history, or claim of a completed white-dwarf cooling track.')
    output.mkdir(parents=True, exist_ok=True)
    with (output/'late_cooling.csv').open('x') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    facts['csv_sha256'] = hashlib.sha256((output/'late_cooling.csv').read_bytes()).hexdigest()
    with (output/'late_cooling_inputs.json').open('x') as stream:
        json.dump(facts, stream, indent=2)
        stream.write('\n')
    print(json.dumps({k:v for k,v in facts.items() if k != 'input_sha256'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--last-run', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--last-checkpoint', type=Path)
    args = parser.parse_args()
    collect(args.last_run, args.output, args.last_checkpoint)
