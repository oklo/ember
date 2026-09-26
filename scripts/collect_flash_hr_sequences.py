"""Collect saved surface states for the two alternative flash histories.

Ordinary evolution contributes each accepted interval endpoint. Short comparison
experiments and changes of mass grid contribute their reviewed restart state.
This is a figure extraction, not a new validation of either complete history.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import gzip
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'out/flash-hr-sequences-sept20-v1'
YEAR = 31557600.
LSUN = 3.828e33
pins = {}


def read(path):
    path = ROOT / path
    raw = path.read_bytes()
    pins[str(path)] = hashlib.sha256(raw).hexdigest()
    return json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)


seed_path = Path('out/pulse-face-matched-bisected-v1/checkpoint_10yr.json')
seed = read(seed_path)['model_record']
epoch = read('out/pulse-matched-face-controls-v1/seed_bisected.json')
epoch = epoch.get('model_record', epoch)['age_seconds']


def point(record, source, kind, years=None):
    surface = record['model'][-1]
    result = dict(age_seconds=record['age_seconds'],
                  pulse_years=(record['age_seconds']-epoch)/YEAR if years is None else years,
                  Teff_K=record['Teff'], log_L_Lsun=math.log10(surface[4]/LSUN),
                  radius_cm=surface[1], nuclear_power_erg_s=record.get('nuclear_luminosity', 0.),
                  source=str(source), kind=kind)
    assert all(math.isfinite(result[k]) for k in ('Teff_K', 'log_L_Lsun', 'pulse_years'))
    return result


def collect_run(root_name, review_name):
    root = ROOT / 'out' / root_name
    report = ROOT / 'docs/results' / (root_name.replace('-', '_') + '.json')
    if not (root/'plan.json').exists() or not report.exists():
        return [], None
    run = read(report)
    if not run['accepted']:
        return [], None
    review = read('docs/results/' + review_name + '.json')
    assert review.get('maximum_error_norm', 0.) <= 1
    plan = read(root/'plan.json')
    parents = [Path(p) for p in plan.get('input_sha256', {}) if '/checkpoint' in p]
    assert parents, root_name
    parent_path = parents[0]
    parent = read(parent_path)
    assert pins[str(parent_path)] == plan['input_sha256'][str(parent_path)]
    first = read(root/'queries'/f"{run['accepted'][0]['queries'][0]:05d}.json.gz")
    tokens = first['query'].split()
    n = int(tokens[5])
    assert float(tokens[1]) == parent['model_record']['age_seconds']
    assert [float(x) for x in tokens[6:6+n*11]] == [x for cell in parent['model_record']['model'] for x in cell]
    # Use integrated local timesteps where available. The common absolute
    # clock locates older controls that used a different naming convention.
    prefix = run.get('prefix_years', (parent['model_record']['age_seconds']-epoch)/YEAR)
    rows = [point(parent['model_record'], parent_path, 'reviewed restart', prefix)]
    for entry in run['accepted']:
        assert entry['accepted'] and entry['error_norm'] <= 1
        number = entry['queries'][-1]
        path = root/'queries'/f'{number:05d}.json.gz'
        raw = read(path)
        assert pins[str(path)] == review['input_sha256'][str(path)], path
        reply = json.loads(raw['raw_reply'])
        assert reply['converged']
        years = prefix + entry['start_elapsed_years'] + entry['dt_years']
        rows.append(point(reply, path, 'accepted interval endpoint', years))
    print(root_name, len(rows), flush=True)
    return rows, dict(run=root_name, review=review_name, points=len(rows),
                     first=rows[0]['pulse_years'], last=rows[-1]['pulse_years'])


def assemble(specs, extra_checkpoints=()):
    rows = [point(seed, seed_path, 'common starting model', 10.)]
    stages = []
    for root, review in specs:
        new, stage = collect_run(root, review)
        rows.extend(new)
        if stage:
            stages.append(stage)
    for p in extra_checkpoints:
        d = read(p)
        if d.get('review'):
            read(d['review'])
        rows.append(point(d['model_record'], p, 'reviewed restart', d.get('total_pulse_years')))
    rows.sort(key=lambda r: r['age_seconds'])
    # Identical boundary states occur at the end of one run and start of the
    # next. Preserve distinct states at the same age (for example a remap).
    unique = []
    for r in rows:
        if unique and all(r[k] == unique[-1][k] for k in ('age_seconds','Teff_K','log_L_Lsun')):
            continue
        unique.append(r)
    return unique, stages


def main(last_onset=28, last_older=39, onset_extra_runs=(), onset_checkpoint=None):
    OUT.mkdir(exist_ok=True)
    onset = [(f'pulse-preonset-current-v{n}', f'pulse_preonset_current_primary_review_v{n}') for n in range(1,last_onset+1)]
    for name in onset_extra_runs:
        prefix = 'pulse-preonset-current-'
        assert name.startswith(prefix) and '/' not in name, name
        onset.append((name, 'pulse_preonset_current_primary_review_' + name[len(prefix):]))
    assert len(onset) == len(set(onset)), 'duplicate onset run'
    old = [(f'pulse-face-matched-bisected-evolution-v{n}', f'pulse_face_matched_bisected_review_v{r}')
           for n,r in [(1,1),(2,3),(3,5),(4,6)]]
    old += [('pulse-hybrid-boundary-bisected-evolution-v1','pulse_hybrid_boundary_bisected_evolution_review_v1')]
    old += [(f'pulse-microscopic-domain-bisected-evolution-v{n}',f'pulse_microscopic_domain_bisected_evolution_review_v{n}') for n in range(1,4)]
    old += [(f'pulse-helium-response-evolution-v{n}',f'pulse_helium_response_evolution_review_v{n}') for n in range(2,21)]
    old += [(f'pulse-front-evolution-v{n}',f'pulse_front_evolution_primary_review_v{n}') for n in range(1,last_older+1)]
    data = {}
    stages = {}
    data['onset'], stages['onset'] = assemble(onset)
    data['older'], stages['older'] = assemble(old, [
        'out/pulse-second-remap-evolution-review-v1/checkpoint.json',
        'out/pulse-front-evolution-review-hot6000_nominal_v2/checkpoint.json'])
    for name, cp in [('onset',onset_checkpoint or f'out/pulse-preonset-current-review-v{last_onset}/checkpoint.json'),
                     ('older',f'out/pulse-front-evolution-review-v{last_older}/checkpoint.json')]:
        record = read(cp)
        final = point(record['model_record'], cp, 'reviewed endpoint', record['total_pulse_years'])
        assert all(data[name][-1][k] == final[k] for k in ['age_seconds','Teff_K','log_L_Lsun'])
        data[name][-1]['pulse_years'] = record['total_pulse_years']
    snapshot = ROOT/'out/full-physics-hr-comparison-v15'
    main_meta = read(snapshot/'inputs.json')
    for source,dest in [('ember.csv','main_sequence.csv'),('pulse_prefix.csv','pre_flash.csv')]:
        raw = (snapshot/source).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == main_meta['snapshot_sha256'][source]
        pins[str(snapshot/source)] = hashlib.sha256(raw).hexdigest()
        (OUT/dest).write_bytes(raw)
    for name, rows in data.items():
        with (OUT/(name+'.csv')).open('w') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
            writer.writeheader(); writer.writerows(rows)
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  endpoints={k:v[-1] for k,v in data.items()},
                  number_of_states={k:len(v) for k,v in data.items()}, stages=stages,
                  common_seed=point(seed,seed_path,'common starting model',10.),
                  main_to_flash_connection='Retain the computed pre-flash prefix. Dotted connections across atmosphere and starting-model changes are not evolved interpolations.',
                  plotted_interval_policy='Each accepted ordinary or compound interval contributes its endpoint. Additional reviewed restart models represent intervening short controls and remaps; their internal time dependence is not drawn as resolved.',
                  scope='Alternative evolved histories, not consecutive portions of a single validated track. No smoothing or invented intermediate stellar states.',
                  input_sha256=pins,
                  csv_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.glob('*.csv')})
    (OUT/'inputs.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['endpoints'],indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--last-onset', type=int, default=28)
    parser.add_argument('--last-older', type=int, default=39)
    parser.add_argument('--onset-extra-run', action='append', default=[],
                        help='Reviewed named continuation belonging to the plotted onset history')
    parser.add_argument('--onset-checkpoint', type=Path,
                        help='Reviewed final checkpoint when the last continuation has a named version')
    parser.add_argument('--output', type=Path, default=OUT)
    args=parser.parse_args()
    OUT=args.output.resolve()
    main(args.last_onset,args.last_older,args.onset_extra_run,args.onset_checkpoint)
