"""Read retained birth-clock histories, preserving reviewed handovers and domain limits."""
import csv
import hashlib
import json
from pathlib import Path
from cold_restart_check import audited_steps, same_restart

def latest_model(path):
    # A worker can be writing its last record while this command reads.
    with path.open('rb') as stream:
        stream.seek(0, 2)
        stream.seek(max(0, stream.tell() - 65536))
        for line in reversed(stream.read().splitlines()):
            try:
                return json.loads(line)
            except json.JSONDecodeError:
                continue
    return None


def frontier(root, mass):
    case = root / 'warm-marcs' / f'extended-start-m{mass:02}'
    if mass == 7:
        case = root / 'm07/extended-start/low-gravity-pilot'
    elif mass == 30:
        case = case / 'low-gravity-pilot'
    patterns = ['evolution', 'contraction*', 'native-continuation/evolution',
                'native-continuation/nuclear-*', 'native-continuation/cold-edge-*',
                'native-continuation/late-overlap', 'native-continuation/gravity-edge', 'native-continuation/gravity-continuation', 'native-continuation/gravity-near-edge']
    for pilot in ('later-evolution-hydrogen-trial', 'native-continuation/gravity-pilot', 'native-continuation/direct-envelope', 'native-continuation/900K-pilot', 'native-continuation/900K-continuation', 'native-continuation/900K-archive', 'native-continuation/900K-near-edge', 'hydrogen-edge-column/pilot', 'hydrogen-edge-column/continuation', 'hydrogen-edge-column/wider-composition/pilot'):
        if (case/pilot/'review.json').is_file():
            patterns.append(pilot+'/evolution')
    candidates = []
    for base in (case, case/'deeper-base'):
        active_patterns = patterns.copy()
        for cold in ('native-continuation/below-1100K',
                     'native-continuation/below-1100K/budget-continuation',
                     'native-continuation/below-1100K/below-1000K'):
            overlap = base/cold/'atmosphere-overlap.json'
            if overlap.is_file() and json.loads(overlap.read_text()).get('passed'):
                boundary_file = base/cold/'first-domain-exit.json'
                if boundary_file.exists():
                    boundary = json.loads(boundary_file.read_bytes())
                    with (base/cold/'evolution/history.jsonl').open('rb') as stream:
                        prefix = stream.read(boundary['history_prefix_bytes'])
                    if hashlib.sha256(prefix).hexdigest() != boundary['history_prefix_sha256']:
                        raise ValueError(f'Atmosphere domain boundary history changed in {cold}')
                    inside, outside = map(json.loads, prefix.splitlines()[-2:])
                    if (inside['Teff_K'], outside['Teff_K'], outside['years']) != (
                            boundary['last_inside_K'], boundary['first_outside_K'], boundary['age_yr']):
                        raise ValueError(f'Atmosphere domain boundary differs from saved history in {cold}')
                    candidates.append(dict(inside, source=str(boundary_file.relative_to(root))))
                else:
                    active_patterns.append(cold+'/evolution')
        for pattern in active_patterns:
            for path in base.glob(pattern+'/history.jsonl'):
                row = latest_model(path)
                if row and row['accepted'] > 0:
                    candidates.append(dict(row, source=str(path.relative_to(root))))
    if mass == 6:
        direct = root/'cool-m06-1000-birth/matched-depth-evolution/constant-depth-to-900/direct-envelope-opacity-check'
        if (direct/'below-925-review.json').is_file():
            review = json.loads((direct/'below-925-review.json').read_bytes())
            segment = Path(review['segments'][-1]['command'][4])
            row = latest_model(segment/'history.jsonl')
            if segment.parent != direct or row != review['last']:
                raise ValueError('Reviewed depth-10 endpoint differs from saved history')
            candidates.append(dict(row, source=str(segment.relative_to(root))))
            extension = direct.parent/'constant-depth-to-850'
            for stage in ('evolution', 'continuation'):
                history = extension/stage/'history.jsonl'
                attempts = history.with_name('attempts.jsonl')
                if not history.exists() or not attempts.exists():
                    break
                def completed_records(path):
                    raw = path.read_bytes()
                    return list(map(json.loads, raw[:raw.rfind(b'\n')+1].splitlines()))
                states = completed_records(history)
                if len(states) < 2:
                    break
                if not same_restart(row, states[0]):
                    raise ValueError(f'Cooling restart differs from saved model in {stage}')
                trials = [t for t in completed_records(attempts) if t['start_years'] < states[-1]['years']]
                if not audited_steps(states, trials, {'accepted_this_invocation':len(states)-1}):
                    break
                row = states[-1]
                candidates.append(dict(row, source=str(history.relative_to(root))))
            # Advertise completed, audited extensions from the same birth clock.
            from cooling_interval import BASE, checked_history
            for colder in [root/BASE, *(root/f'm06/cooling-intervals/band-{i}/native' for i in (1,2))]:
                if not (colder/'native-review.json').is_file():
                    break
                states, _ = checked_history(colder)
                if not same_restart(row, states[0]):
                    raise ValueError(f'Cooling interval ancestry differs at {colder}')
                row = states[-1]
                candidates.append(dict(row, source=str((colder/'native-review.json').relative_to(root)),
                    qualification='Native birth clock with measured atmosphere knots and provisional gravity and composition transfer.'))
    review_file = case/'main-sequence-join-review.json'
    track = case/'birth-clock-main-sequence.csv'
    if review_file.is_file() and track.is_file():
        review = json.loads(review_file.read_bytes())
        with track.open() as stream:
            last = None
            for last in csv.DictReader(stream):
                pass
        if last is None or float(last['birth_years']) != review['through_birth_years']:
            raise ValueError(f'Joined track differs from its review in {case}')
        if mass == 8:
            for stage in ('checked-near-edge', 'checked-continuation'):
                extended = case/f'birth-clock-{stage}.csv'
                if extended.is_file():
                    with extended.open() as stream:
                        for last in csv.DictReader(stream):
                            pass
                    checked = json.loads((root/'helium-control-m08'/f'{stage}.json').read_bytes())
                    expected = checked['last']['years'] + review['offset_years']
                    if float(last['birth_years']) != expected or float(last['Teff_K']) != checked['last']['Teff_K']:
                        raise ValueError(f'Joined extension differs from reviewed native endpoint in {extended}')
                    track = extended
                    review['limits'] = 'Conditional on the early atmosphere, main-sequence join and checked surface-composition interval.'
                    break
        candidates.append(dict(years=float(last['birth_years']), Teff_K=float(last['Teff_K']), joined=True, source=str(track.relative_to(root)), qualification=review["limits"]))

    if candidates:
        row = max(candidates, key=lambda item: item['years'])
        return {k: row[k] for k in ('years', 'Teff_K', 'source', 'qualification', 'nuclear_fraction', 'central_X') if k in row}
    return None
