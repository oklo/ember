"""Choose and evolve measured cooling intervals on the 0.06 solar mass birth clock."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

from cold_restart_check import audited_steps, same_restart

# Spacing criteria for a cooling pilot. These do not qualify a tail for release.
LIMITS = dict(elapsed_relative=.10, core_relative=.03, radius_relative=.01, nuclear_fraction_difference=.01)
MATCH_LIMITS = (.003, .01, .02)
BASE = Path('cool-m06-1000-birth/matched-depth-evolution/constant-depth-to-900/constant-depth-to-850/constant-depth-to-800')
FORECAST = Path('thermal-holdout/m06/birth-depth-10-900K-source/forecast-800')
FROZEN = '54264c139934a38da01dd9706054c1e688e7d21080ad64eb79bc1f07672b4e08'


def read(path):
    return json.loads(path.read_bytes())


def save(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    return [json.loads(line) for line in path.read_bytes().splitlines()]


def checked_history(case):
    review = read(case/'native-review.json')
    seed = Path(review['anchor']['checkpoint'])
    assert digest(seed) == review['anchor']['checkpoint_sha256']
    previous = rows(seed.parent/'history.jsonl')[-1]
    samples, hashes = [], {}
    for segment in map(Path, review['segments']):
        assert segment.parent == case
        h, a, r = (segment/name for name in ('history.jsonl', 'attempts.jsonl', 'report.json'))
        states = rows(h)
        assert same_restart(previous, states[0]) and audited_steps(states, rows(a), read(r))
        samples.extend(states[1:] if samples else states)
        previous = states[-1]
        hashes.update({str(p): digest(p) for p in (h, a, r)})
    assert samples and samples[-1] == review['last']
    return samples, hashes


def crossing(samples, temp, fields):
    for a, b in zip(samples, samples[1:]):
        if a['Teff_K'] >= temp >= b['Teff_K'] and a['Teff_K'] > b['Teff_K']:
            w = (a['Teff_K']-temp)/(a['Teff_K']-b['Teff_K'])
            return [a[k] + w*(b[k]-a[k]) for k in fields]
    raise ValueError(f'Missing downward crossing of {temp} K')


def cooling_csv(path):
    with path.open() as stream:
        return [{k:float(v) for k,v in row.items()} for row in csv.DictReader(stream)]


def freeze(root):
    """Extend the already registered forecast with identical inputs and verify its prefix."""
    old = root/FORECAST/'cooling.csv'
    assert digest(old) == FROZEN
    job = next(j for j in read(root/'recovery-plan.json') if j['name']=='m06-900K-anchor-forecast-to-800')
    cmd = job['command'].copy()
    assert cmd[2] == str(old) and cmd[-1] == 'Teff_stop=798'
    dest = root/'m06/cooling-intervals'
    dest.mkdir(exist_ok=True)
    path = dest/'forecast.csv'
    cmd[2], cmd[-1] = str(path), 'Teff_stop=698'
    record = dest/'forecast.json'
    if path.is_file() and record.exists() and read(record)['command'] == cmd and read(record)['sha256'] == digest(path):
        return
    subprocess.run(cmd, check=True, capture_output=True)
    prefix = old.read_bytes().splitlines()
    assert path.read_bytes().splitlines()[:len(prefix)] == prefix
    assert cooling_csv(path)[-1]['Teff_K'] <= 700
    save(dest/'forecast.json', dict(command=cmd, sha256=digest(path), prefix_sha256=FROZEN))


def compare(native, reduced, warm, cold):
    nf = ('years','central_T_K','radius_Rsun','nuclear_fraction')
    rf = ('age_yr','Tc_K','R_ph_cm','Lnuc_over_L')
    # Use the actual restart temperature at an edge slightly above/below its nominal value.
    start = min(native[0]['Teff_K'], reduced[0]['Teff_K'])
    assert abs(start-warm) <= 1
    na = crossing(native, start, nf)[0]
    ra = crossing(reduced, start, rf)[0]
    rb, tc, radius, nuc = crossing(reduced, cold, rf)
    age, central, rsol, frac = crossing(native, cold, nf)
    dt = age-na
    assert dt > 0 and rb > ra
    errors = dict(elapsed_relative=(rb-ra)/dt-1, core_relative=tc/central-1,
                  radius_relative=radius/(rsol*6.957e10)-1, nuclear_fraction_difference=nuc-frac)
    assert all(math.isfinite(v) for v in errors.values())
    # Bound the elapsed-time ratio with the timesteps bracketing each crossing.
    def age_bounds(series, temperature, key):
        for left,right in zip(series,series[1:]):
            if left['Teff_K'] >= temperature >= right['Teff_K']:
                if left['Teff_K']==temperature:
                    return left[key],left[key]
                if right['Teff_K']==temperature:
                    return right[key],right[key]
                return left[key],right[key]
        raise ValueError(f'Missing crossing time bounds at {temperature} K')

    def elapsed_bounds(series, key):
        a,b = age_bounds(series,start,key)
        c,d = age_bounds(series,cold,key)
        return c-b,d-a

    nmin,nmax = elapsed_bounds(native,'years')
    rmin,rmax = elapsed_bounds(reduced,'age_yr')
    bounded = nmin>0 and rmin>0
    bounds = [rmin/nmax-1,rmax/nmin-1] if bounded else None
    return dict(warm_K=start,crossing_K=cold,native_elapsed_yr=dt,reduced_elapsed_yr=rb-ra,
                errors=errors,elapsed_error_bracket=bounds)


def decide(root, index, case, stage):
    previous = read(case/'native-review.json')
    native, hashes = checked_history(case)
    if index == 1:
        warm, edge, comparison = 850, 800, 825
        check = read(root/'cool-m06-850-birth-depth10/boundary-comparison.json')
        atmosphere_ok = check['within_campaign_tolerances']
        assert read(case/'clock-comparison.json')['points']['825']['native_age_yr'] > native[0]['years']
    else:
        last = read(stage.parent/f'band-{index-1}'/'decision.json')
        warm, edge, comparison = last['edge_K'], last['target_K'], last['midpoint_K']
        atmosphere_ok = read(case/'reader-check.json')['single_interval_accurate']
    assert edge <= native[-1]['Teff_K']+1e-5 <= edge+1.00001
    forecast = stage.parent/'forecast.csv'
    assert digest(forecast) == read(stage.parent/'forecast.json')['sha256']
    score = compare(native, cooling_csv(forecast), warm, comparison)
    # A resolved surface reheating or nuclear support calls for closer spacing.
    monotone = all(b['Teff_K'] <= a['Teff_K'] + .01 for a,b in zip(native,native[1:]))
    weak_burning = max(r['nuclear_fraction'] for r in native) < .1
    passed = all(abs(v) <= LIMITS[k] for k,v in score['errors'].items())
    bounds = score['elapsed_error_bracket']
    resolved = bounds is not None and max(map(abs,bounds)) <= LIMITS['elapsed_relative']
    spacing = 50 if passed and resolved and monotone and weak_burning and atmosphere_ok else 25
    interpolate=False
    validation={}
    if index==2:
        for study in ('m05','m06','band1'):
            path=root/'atmosphere-validation'/study/'comparison.json'
            if path.is_file():
                validation[study]=dict(passes=read(path)['passes'],sha256=digest(path))
        interpolate=(len(validation)==3 and all(r['passes'] for r in validation.values()) and
                     atmosphere_ok and monotone and weak_burning and spacing<=warm-edge)
    save(stage/'decision.json', dict(interpolate_midpoint=interpolate,validation=validation,edge_K=edge, target_K=edge-spacing, midpoint_K=edge-spacing//2,
        comparison=score, limits=LIMITS, atmosphere_ok=atmosphere_ok, matching_limits=MATCH_LIMITS,
        cooling_without_resolved_reversal=monotone, nuclear_fraction_below_one_tenth=weak_burning,
        reference='native birth clock', previous_case=str(case), seed=previous['segments'][-1]+'/final.checkpoint',
        source_sha256=hashes, forecast_sha256=digest(forecast)))


def prepare(root, index, stage):
    decision = read(stage/'decision.json')
    ck = Path(decision['seed'])
    native = rows(ck.parent/'history.jsonl')[-1]
    ck_rows = ck.read_bytes().splitlines()
    off = 9+int(ck_rows[8])
    count, mass, seconds = ck_rows[off].split()[:3]
    assert abs(float(seconds)/31557600 - native['years']) < 1
    surface = list(map(float, ck_rows[off+int(count)].split()))
    x,h3,h4 = surface[7:10]
    z = sum(surface[10:15])
    assert abs(x-native['surface_X']) < 1e-12 and abs(z-native['surface_Z']) < 1e-12
    template = root/'cool-m06-800-birth-depth10'
    seed = (read(template/'measured-800K.json')['receipt'] if index==1 else
            read(Path(decision['previous_case'])/'reader-check.json')['columns']['endpoint']['receipt'])
    assert read(Path(seed))['quality']['usable']
    knots = [('verify', 'edge_K')] if index == 3 else [('midpoint','midpoint_K'), ('endpoint','target_K')]
    for name,key in knots:
        directory = stage/name
        (directory/'grid_hhe').mkdir(parents=True, exist_ok=True)
        spec = read(template/'spec.json')
        spec.update(teff=decision[key], years=native['years'], tag=f'_m06_band_{index}_{name}',
            checkpoint=str(ck), checkpoint_sha256=digest(ck), seed_receipt=seed, Z_over_X=z/x,
            arms=dict(anchor=dict(X=x, He3=h3, He4=h4, Z=z)),
            logg=math.log10(6.67430e-8*float(mass)/(native['radius_Rsun']*6.957e10)**2),
            purpose='Measured cooling interval at the native restart gravity and composition.')
        from ember_convection_seed import neighbor_regions
        regions=neighbor_regions(read(Path(seed)),spec)
        driver=template/'solve.py'
        if regions:
            spec['initial_regions']=regions
            driver=root/'cool-m06-850/detached-region/solve.py'
        shutil.copyfile(driver, directory/'solve.py')
        save(directory/'spec.json', spec)


def measurements(root, index, stage, d, knots):
    from archive_matching_points import matching_points
    columns = {}
    for name,key in knots:
        spec = read(stage/name/'spec.json')
        chosen = read(root/'recovery'/f'm06-band-{index}-{name}'/'selected-column.json')
        path = Path(chosen['receipt'])
        receipt = read(path)
        actual = receipt['spec']
        numerical = {'resume_guess','resume_region_top','initial_regions'}
        assert actual == read(path.parent/'spec.json') and path.parent == Path(chosen['directory'])
        assert {k:v for k,v in actual.items() if k not in numerical} == {k:v for k,v in spec.items() if k not in numerical}
        assert receipt['quality']['usable'] and abs(receipt['teff_used']-d[key]) < .01
        assert spec['checkpoint_sha256'] == digest(Path(d['seed'])) and spec['checkpoint']==d['seed']
        layer = matching_points(receipt)['10']
        measured = [layer['T'], 1e6*layer['P_bar'], layer['depth_cm_mu_He3corr']]
        assert all(math.isfinite(v) and v>0 for v in measured)
        columns[name] = dict(spec=spec, receipt=str(path), sha256=digest(path), measured=measured)
    return columns


def probe(root, family, temp, spec):
    args = read(root.parent/'highz-native-opacity-oct12-v1/matching_depth_controls.json')['command'][1:]
    mix = spec['arms']['anchor']
    query = f"thin {temp} {spec['logg']} {mix['X']} {mix['He3']} {mix['Z']}\n"
    output = subprocess.run([root.parent/'boundary-corrected-item1-oct13/probe',*args,family],
                            input=query, text=True, capture_output=True, check=True)
    return list(map(float,output.stdout.split()[6:9]))


def build_family(root, index, stage):
    from extend_measured_boundary import extend
    d = read(stage/'decision.json')
    source = Path(d['previous_case'])
    if d.get('interpolate_midpoint'):
        assert index==2
        for study,check in d['validation'].items():
            path=root/'atmosphere-validation'/study/'comparison.json'
            assert digest(path)==check['sha256'] and read(path)['passes']
        columns=measurements(root,index,stage,d,[('endpoint','target_K')])
        end=columns['endpoint']
        case=stage/'native'
        extend(source,case,end['spec'],d['target_K'],d['edge_K'],end['measured'])
        save(case/'reader-check.json',dict(columns=columns,interpolated_midpoint=True,
            single_interval_accurate=False,validation=d['validation'],
            note='Bracketed by measured endpoints. Final evolved-surface column remains required.'))
        return
    columns = measurements(root, index, stage, d, [('midpoint','midpoint_K'), ('endpoint','target_K')])
    mid, end = columns['midpoint'], columns['endpoint']
    # Probe the single-interval prediction at the actual integer midpoint.
    temp_table = stage/'endpoint-only'
    extend(source,temp_table,end['spec'],d['target_K'],d['edge_K'],end['measured'])
    predicted = probe(root,temp_table/'family.dat',d['midpoint_K'],mid['spec'])
    errors = [v/m-1 for v,m in zip(predicted,mid['measured'],strict=True)]
    middle = stage/'midpoint-family'
    extend(source,middle,mid['spec'],d['midpoint_K'],d['edge_K'],mid['measured'])
    case = stage/'native'
    extend(middle,case,end['spec'],d['target_K'],d['midpoint_K'],end['measured'])
    save(case/'reader-check.json',dict(columns=columns, relative_T_P_depth=errors,
        single_interval_accurate=all(abs(e)<=l for e,l in zip(errors,MATCH_LIMITS,strict=True)),
        note='Both measured knots are used. Gravity and composition dependence away from the restart is transferred from the warmer table.'))


def evolve(root, stage):
    from atmosphere_repair import members
    from coordinate import solver
    d = read(stage/'decision.json')
    case, source, ck = stage/'native', Path(d['previous_case']), Path(d['seed'])
    spec = read(case/'reader-check.json')['columns']['endpoint']['spec']
    assert digest(ck) == spec['checkpoint_sha256']
    old, new = source/'family.dat', case/'family.dat'
    cfg = (source/'config.txt').read_text()
    assert cfg.count(str(old)) == 1
    (case/'config.txt').write_text(cfg.replace(str(old), str(new)))
    raw = ck.read_bytes().splitlines(keepends=True)
    off = 9+int(raw[8])
    ids = dict(shlex.split(row.decode()) for row in raw[9:off])
    for key,value in members(old).items():
        assert ids.pop(key) == value
    ids.update(members(new))
    start = case/'start.checkpoint'
    entries = ''.join(f'{json.dumps(k)} {json.dumps(v)}\n' for k,v in sorted(ids.items()))
    start.write_bytes(b''.join(raw[:8])+f'{len(ids)}\n'.encode()+entries.encode()+b''.join(raw[off:]))
    shutil.copyfile(str(ck)+'.lithium', str(start)+'.lithium')
    previous = rows(ck.parent/'history.jsonl')[-1]
    segments = []
    accepted = 0
    for n in range(4):
        seg = case/f'segment-{n}'
        cmd = [str(solver),'--lifetime',str(case/'config.txt'),str(max(4e10,3*spec['years'])),str(seg),
               '--restart',str(start),'--cpu-seconds','900']
        if not n:
            cmd.extend(['--restart-step-years','10000'])
        if not (seg/'report.json').exists():
            with (case/f'segment-{n}.log').open('w') as log:
                result = subprocess.run(cmd, stdout=log, stderr=log)
            if not (seg/'report.json').exists():
                result.check_returncode()
        report = read(seg/'report.json')
        states = rows(seg/'history.jsonl')
        assert same_restart(previous,states[0]) and audited_steps(states,rows(seg/'attempts.jsonl'),report)
        previous = states[-1]
        accepted += len(states)-1
        segments.append(str(seg))
        start = seg/'final.checkpoint'
        stop = report['stop_reason'].strip().lower()
        reached_edge = 'outside' in stop and 'teff' in stop and 'table' in stop
        if reached_edge and d['target_K']-1e-5 <= previous['Teff_K'] <= d['target_K']+1:
            assert start.is_file() and Path(str(start)+'.lithium').is_file()
            break
        if report['stop_reason'].strip().lower() != 'cpu budget':
            raise ValueError(f'Native evolution stopped inside the measured interval: {report["stop_reason"]}')
    else:
        raise ValueError('Native pilot used its four CPU budgets before reaching the measured edge')
    save(case/'native-review.json',dict(anchor=dict(checkpoint=str(ck),checkpoint_sha256=digest(ck)),
        accepted=accepted,segments=segments,last=previous,stop=report['stop_reason'],target_K=d['target_K']))


def verify(root, stage):
    d = read(stage/'decision.json')
    column = measurements(root,3,stage,d,[('verify','edge_K')])['verify']
    predicted = probe(root,Path(d['previous_case'])/'family.dat',d['edge_K'],column['spec'])
    errors = [p/m-1 for p,m in zip(predicted,column['measured'],strict=True)]
    passed = all(abs(e)<=l for e,l in zip(errors,MATCH_LIMITS,strict=True))
    save(stage/'surface-review.json',dict(edge_K=d['edge_K'],column=column,
        relative_T_P_depth=errors,tolerances=MATCH_LIMITS,passes=passed,clock_comparison=d['comparison'],
        next_action=(f"Measure the next {d['edge_K']-d['target_K']} K native interval." if passed else
                     "Refit the atmosphere response at both measured surface compositions and rerun the last native interval.")))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('root', type=Path)
    parser.add_argument('index', type=int, choices=[1,2,3])
    parser.add_argument('action', choices=['freeze','decide','prepare','family','native','verify'])
    args = parser.parse_args()
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
    root = args.root.resolve()
    os.environ['EMBER_GRID_ROOT'] = str(root)
    sys.path.append(str(root/'thermal-holdout/m06'))
    parent = root/'m06/cooling-intervals'
    stage = parent/f'band-{args.index}'
    if args.action == 'freeze':
        freeze(root)
    elif args.action == 'decide':
        case = root/BASE if args.index == 1 else parent/f'band-{args.index-1}'/'native'
        decide(root,args.index,case,stage)
    elif args.action == 'prepare':
        prepare(root,args.index,stage)
    elif args.action == 'family':
        build_family(root,args.index,stage)
    elif args.action == 'verify':
        verify(root,stage)
    else:
        evolve(root,stage)
