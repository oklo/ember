"""Check reviewed track evidence for every mass in the solar campaign."""
import argparse
import csv
import hashlib
import json
from itertools import zip_longest
import math
from pathlib import Path
import shlex
import sys

SOLAR_Z = 0.016482259798651158
# The campaign's z1/composition.dat tabulates the mixture after deuterium burning.
SOLAR_POST_D = (0.70567931775974979, 7.8560545197511764e-05, 0.27775986189640156,
                0.0033193055090442805, 3.729556751735148e-05, 0.00091762595953452165,
                0.00759300001893268, 0.0046150327436223241, 0.0)
SOLAR_BIRTH_D = 2.9639153768138622e-05
SOLAR_MASS_G = 1.3271244e26 / 6.67430e-8


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify_solar_mixture(mixture):
    require(len(mixture) == 9, 'incorrect number of initial isotope fractions')
    fuel = mixture[-1]
    require(any(abs(fuel-value) < 1e-14 for value in (0.0, SOLAR_BIRTH_D)),
            'unrecognized initial deuterium abundance')
    post = list(mixture)
    post[0] -= fuel / 2
    post[1] += 1.5 * fuel
    post[-1] = 0.0
    require(all(abs(a-b) < 1e-14 for a,b in zip(post,SOLAR_POST_D)), 'wrong solar isotope mixture')


def verify_birth_state(path, mass, mixture, row):
    """Tie the first source model to the mass and composition in its saved interior."""
    lines = path.read_text().splitlines()
    require(len(lines) > 10 and lines[0] == 'EMBER_EVOLUTION_CHECKPOINT 7'
            and lines[-1] == 'END', 'invalid birth checkpoint')
    header = 9 + int(lines[8])
    require(9 <= header < len(lines)-1, 'invalid checkpoint input count')
    cells, grams, seconds, *_ = lines[header].split()
    n = int(cells)
    require(n > 0 and header+n < len(lines)-1, 'truncated birth model')
    require(float(seconds) == 0 and math.isclose(float(grams), mass*SOLAR_MASS_G, rel_tol=1e-10),
            'birth checkpoint age or mass mismatch')
    require(len(mixture) == 9 and all(0 <= x <= 1 for x in mixture)
            and abs(sum(mixture)-1) < 1e-10, 'invalid initial isotope fractions')
    center = lines[header+1].split()
    surface = lines[header+n].split()
    for line in lines[header+1:header+n+1]:
        zone = line.split()
        require(len(zone) >= 16 and all(abs(float(x)-y) < 1e-10
                for x,y in zip(zone[7:16], mixture)), 'birth checkpoint mixture mismatch')
    temperature = float(row['central_T_K'])
    if mixture[-1] or 'central_density_g_cm3' in row:
        density = float(row['central_density_g_cm3'])
        require(math.isfinite(density) and density > 0
                and abs(math.log(density)-float(center[2])) < 1e-10,
                'initial central density differs from birth checkpoint')
    # Older compact histories expose fewer diagnostics than the deuterium-bearing runs.
    for field, fraction in (('H_mass_g', mixture[0]), ('He3_mass_g', mixture[1]), ('D_mass_g', mixture[-1])):
        if mixture[-1] or field in row:
            require(field in row and math.isclose(float(row[field]), float(grams)*fraction, rel_tol=1e-9),
                    f'initial history {field} differs from the uniform birth mixture')
    require(math.isfinite(temperature) and temperature > 0
            and abs(math.log(temperature)-float(center[3])) < 1e-10
            and abs(float(row['surface_X'])-float(surface[7])) < 1e-10,
            'first history model differs from birth checkpoint')


def source_rows(segments, files):
    """Check the join on the original grids before yielding each reviewed interval."""
    previous = following = None
    previous_native = None
    cached_stage = None
    stage_models = {}
    for segment in segments:
        path = files[segment['source']]
        start, stop = segment['start'], segment['stop']
        require(isinstance(start, int) and isinstance(stop, int) and 0 <= start < stop,
                'invalid source row interval')
        count = 0
        preceding = preceding_native = None
        prior_csv_step = None
        with path.open() as stream:
            records = (map(json.loads, stream) if path.suffix == '.jsonl' else csv.DictReader(stream))
            for index, row in enumerate(records):
                point = float(row[segment['age_column']]), float(row['Teff_K'])
                if index == stop:
                    following = point
                    break
                if index < start:
                    preceding = point
                    preceding_native = row if path.suffix == '.jsonl' else None
                    continue
                if index == start and previous is not None:
                    # Restarts commonly omit their duplicate initial row. A splice inside
                    # a previous step can instead be checked against its two bounding rows.
                    shared = preceding is not None and preceding[0] == previous[0] and math.isclose(
                        preceding[1], previous[1], rel_tol=3e-6)
                    if previous_native is not None and path.suffix == '.jsonl':
                        require(shared and preceding_native is not None,
                                f'native restart lacks a shared checkpoint state at {path}')
                        fields = ['central_X', 'H_mass_g', 'He3_mass_g', 'surface_X', 'surface_Z',
                                  'central_T_K', 'luminosity_Lsun']
                        if 'D_mass_g' in previous_native or 'D_mass_g' in preceding_native:
                            require('D_mass_g' in previous_native and 'D_mass_g' in preceding_native,
                                    f'native restart lost deuterium inventory at {path}')
                            fields.append('D_mass_g')
                        for field in fields:
                            require(math.isclose(float(previous_native[field]), float(preceding_native[field]),
                                                 rel_tol=1e-9, abs_tol=0),
                                    f'native restart changed {field} at {path}')
                    if not shared and following is not None and previous[0] < point[0] <= following[0]:
                        fraction = (point[0] - previous[0]) / (following[0] - previous[0])
                        expected = previous[1] + fraction * (following[1] - previous[1])
                        shared = math.isclose(point[1], expected, rel_tol=3e-6)
                    require(shared, f'no matching source state at join to {path.name}')
                method = segment['method']
                if 'stage_column' in segment:
                    method += row[segment['stage_column']]
                if method.startswith('Ember:') and path.suffix != '.jsonl':
                    # An assembled CSV can contain both native and reduced stages.
                    # Native rows must reproduce the corresponding solver history.
                    stage = row.get(segment.get('stage_column')) or segment.get('native_source')
                    require(stage in files and files[stage].suffix == '.jsonl',
                            f'missing native solver history for {method}')
                    if stage != cached_stage:
                        with files[stage].open() as native_stream:
                            stage_models = {}
                            for line in native_stream:
                                model = json.loads(line)
                                require(model['years'] not in stage_models, f'duplicate native age in {stage}')
                                stage_models[model['years']] = model
                        cached_stage = stage
                    model = stage_models.get(point[0])
                    require(model is not None, f'age missing from native stage {stage}')
                    if prior_csv_step is not None:
                        require(model['accepted'] == prior_csv_step + 1,
                                f'missing accepted native step in {path.name} at row {index}')
                    prior_csv_step = model['accepted']
                    for key in ('Teff_K', 'central_T_K', 'radius_Rsun', 'central_X', 'surface_X',
                                'surface_Z', 'luminosity_Lsun', 'nuclear_fraction'):
                        require(math.isclose(float(row[key]), model[key], rel_tol=1e-10),
                                f'{key} differs from native stage {stage} at age {point[0]}')
                if path.suffix == '.jsonl' and method.startswith('Ember:'):
                    require(point[0] == row['years'], f'native model uses a different age at {path}:{index}')
                    prior = preceding_native if index == start else previous_native
                    if prior is not None:
                        require(row['accepted'] == prior['accepted'] + 1,
                                f'missing accepted native step in {path.name} at row {index}')
                if method.startswith('Ember:'):
                    native_model = row if path.suffix == '.jsonl' else model
                    luminosity = float(native_model['luminosity_Lsun']) * 3.828e33
                    nuclear = luminosity * float(native_model['nuclear_fraction'])
                    for key, expected in [('L_erg_s', luminosity), ('Lnuc_erg_s', nuclear)]:
                        if key in row:
                            require(math.isclose(float(row[key]), expected, rel_tol=1e-10),
                                    f'{key} differs from native luminosity at {path.name}:{index}')
                elif 'L_erg_s' in row:
                    luminosity = float(row['L_erg_s'])
                    nuclear = float(row['Lnuc_erg_s'])
                else:
                    luminosity = float(row['luminosity_Lsun']) * 3.828e33
                    nuclear = luminosity * float(row['nuclear_fraction'])
                require(math.isfinite(luminosity) and luminosity > 0,
                        f'invalid surface luminosity in {path.name} at row {index}')
                yield point[0], point[1], method, luminosity, nuclear
                previous = point
                previous_native = row if path.suffix == '.jsonl' else None
                following = None
                count += 1
        require(count == stop - start, f'truncated source interval in {path.name}')


def contraction_review(review, files, mixture):
    """Confirm that the reviewed birth precedes a measured interval of contraction."""
    name = review.get('initial_contraction_review')
    if not name:
        return 'initial contraction has not been reviewed'
    require(name in files, 'missing initial contraction review file')
    record = json.loads(files[name].read_text())
    if record['decision'] == 'pending':
        require(bool(record['reason'].strip()), 'missing birth qualification')
        return record['reason']
    require(record['decision'] == 'accepted' and record['reviewer'].strip(), 'unreviewed start')
    require(math.isclose(mixture[-1], SOLAR_BIRTH_D, rel_tol=0, abs_tol=1e-14),
            'end-to-end solar track must start with the birth deuterium abundance')
    segment = review['segments'][0]
    require(record['source'] == segment['source'], 'birth review refers to a different history')
    require(record['initial_state_justification'].strip() and record['deuterium_treatment'].strip(),
            'missing initial model qualifications')
    index = record['contracted_row']
    require(segment['start'] < index < segment['stop'], 'contraction interval outside first source')
    path = files[segment['source']]
    with path.open() as stream:
        rows = map(json.loads, stream) if path.suffix == '.jsonl' else csv.DictReader(stream)
        for i, row in enumerate(rows):
            if i == segment['start']:
                birth = row
            if i == index:
                require(float(row[segment['age_column']]) > float(birth[segment['age_column']])
                        and 0 < float(row['radius_Rsun']) < float(birth['radius_Rsun']),
                        'recorded interval does not show contraction')
                for key in ('radius_Rsun', 'central_T_K', 'Teff_K'):
                    require(math.isclose(float(birth[key]), record['birth'][key], rel_tol=1e-10),
                            f'birth review differs from source {key}')
                return None
    raise ValueError('truncated birth evidence')


def verify_mass(review_file, mass):
    review = json.loads(review_file.read_text())
    require(review['mass_Msun'] == mass and review['initial_Z'] == SOLAR_Z, 'wrong mass or mixture')
    files = {}
    for name, evidence in review['files'].items():
        path = review_file.parent / evidence['path']
        with path.open('rb') as stream:
            digest = hashlib.file_digest(stream, 'sha256').hexdigest()
        require(digest == evidence['sha256'], f'{name} has changed since review')
        files[name] = path
    cfg = dict(shlex.split(line) for line in files['birth_config'].read_text().splitlines() if line.strip())
    require(float(cfg['mass_Msun']) == mass, 'birth mass mismatch')
    require(Path(cfg['composition']).resolve() == files['composition'].resolve(), 'wrong composition source')
    mixture = list(map(float, files['composition'].read_text().split()))
    verify_solar_mixture(mixture)
    require(abs(sum(mixture[3:8]) - SOLAR_Z) < 1e-14, 'birth metallicity mismatch')
    require(bool(files['physical_review'].read_text().strip()), 'missing physical review')
    previous = None
    crossing = None
    methods = set()
    rows = 0
    with files['history'].open() as stream:
        segments = review['segments']
        require(bool(segments) and all(files[s['source']].resolve() != files['history'].resolve()
                                       for s in segments), 'missing independent source history')
        first = segments[0]
        require(first['start'] == 0, 'native birth must start at first source model')
        path = files[first['source']]
        with path.open() as birth_stream:
            birth_model = (json.loads(birth_stream.readline()) if path.suffix == '.jsonl'
                           else next(csv.DictReader(birth_stream), {}))
        verify_birth_state(files['birth_checkpoint'], mass, mixture, birth_model)
        for row, source in zip_longest(csv.DictReader(stream), source_rows(segments, files)):
            require(row is not None and source is not None, 'history length differs from source')
            age, temp = float(row['age_yr']), float(row['Teff_K'])
            require(math.isfinite(age) and math.isfinite(temp) and temp > 0, 'invalid age or temperature')
            require(row['method'].strip(), 'missing method')
            require(source[:3] == (age, temp, row['method']), 'history differs from reviewed source interval')
            for column, expected in zip(('L_erg_s', 'Lnuc_erg_s'), source[3:]):
                value = float(row[column])
                require(math.isfinite(value) and math.isclose(value, expected, rel_tol=1e-10),
                        f'{column} differs from reviewed source interval')
            if previous is None:
                require(age == 0 and temp > 1 and row['method'].startswith('Ember:'), 'missing native birth')
            else:
                require(age > previous[0], 'unordered history')
                if crossing is None and previous[1] > 1 >= temp:
                    crossing = [previous[0], age]
            methods.add(row['method'])
            previous = age, temp
            rows += 1
    require(crossing is not None, 'no downward crossing of 1 K')
    reduced = any(not method.startswith('Ember:') for method in methods)
    require(set(review['reviewed_methods']) == methods, 'unreviewed method')
    for key in ['source_provenance', 'continuity_checks', 'reproduction'] + (['sensitivity_checks'] if reduced else []):
        require(review[key] and all(name in files for name in review[key]), f'missing {key}')
    if reduced:
        require(review.get('endpoint_controls') and 'control_first_1K_age_range_years' in review,
                'missing 1 K crossing evidence for reduced-model sensitivity checks')
    distinct_controls = 0
    if 'control_first_1K_age_range_years' in review:
        controls = review['endpoint_controls']
        require(controls and set(controls) <= set(review['sensitivity_checks']),
                'unreviewed endpoint control')
        nominal = {review['files'][s['source']]['sha256'] for s in segments}
        nominal.add(review['files']['history']['sha256'])
        distinct_controls = len({review['files'][k]['sha256'] for k in controls} - nominal)
        require(distinct_controls > 0, 'endpoint controls only repeat the nominal history')
        crossings = []
        for name in controls:
            prior = crossing_age = None
            with files[name].open() as stream:
                for row in csv.DictReader(stream):
                    age, temp = float(row['age_yr']), float(row['Teff_K'])
                    require(math.isfinite(age) and math.isfinite(temp) and age >= 0 and temp > 0,
                            f'invalid control model in {name}')
                    if 'L_erg_s' in row:
                        luminosity = float(row['L_erg_s'])
                        require(math.isfinite(luminosity) and luminosity > 0,
                                f'invalid surface luminosity in control {name}')
                    if 'Lnuc_erg_s' in row:
                        require(math.isfinite(float(row['Lnuc_erg_s'])),
                                f'invalid nuclear luminosity in control {name}')
                    if prior is None:
                        require(temp > 1, f'control {name} starts at or below the endpoint')
                    else:
                        a, t = prior
                        require(age > a, f'unordered control {name}')
                        if crossing_age is None and t > 1 >= temp:
                            crossing_age = a+(age-a)*(t-1)/(t-temp)
                    prior = age, temp
            require(crossing_age is not None, f'no 1 K crossing in control {name}')
            crossings.append(crossing_age)
        reported = review['control_first_1K_age_range_years']
        require(len(reported) == 2 and all(math.isclose(a,b,rel_tol=1e-12)
                for a,b in zip(reported,[min(crossings),max(crossings)])), 'control age range mismatch')
    birth_pending = contraction_review(review, files, mixture)
    if review['decision'] == 'pending' and review.get('pending_reason'):
        raise ValueError(f"pending scientific review — {review['pending_reason']}")
    require(review['decision'] == 'accepted' and review['reviewer'].strip(), 'pending scientific review')
    return dict(method='reduced' if reduced else 'native', crossing_age_bracket_yr=crossing, rows=rows,
                distinct_non_nominal_control_files=distinct_controls,
                initial_contraction_pending=birth_pending)


def audit(root):
    campaign = json.loads((root / 'campaign.json').read_text())
    masses = [m['mass_Msun'] for m in campaign['tracks']]
    require(sorted(masses) == [i / 100 for i in range(1, 31)], 'incorrect mass manifest')
    require(campaign['initial_Z'] == SOLAR_Z and campaign['endpoint_Teff_K'] == 1, 'incorrect campaign target')
    tracks = {}
    for mass in masses:
        review = root / f'm{round(100*mass):02d}' / 'endpoint-review.json'
        try:
            if not review.exists():
                tracks[f'{mass:.2f}'] = {'pending': 'no reviewed endpoint'}
                continue
            tracks[f'{mass:.2f}'] = verify_mass(review, mass)
        except (OSError, KeyError, TypeError, ValueError) as error:
            tracks[f'{mass:.2f}'] = {'pending': str(error)}
    return tracks


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign', type=Path)
    parser.add_argument('--json', action='store_true')
    args = parser.parse_args()
    tracks = audit(args.campaign)
    passed = sum('method' in row for row in tracks.values())
    complete = sum('method' in row and row['initial_contraction_pending'] is None for row in tracks.values())
    print(json.dumps(tracks, indent=2) if args.json else
          f'{passed}/30 reviewed 1 K endpoints; {complete}/30 reviewed end-to-end tracks')
    sys.exit(0 if complete == 30 else 1)
