"""Campaign release policy, bounded work orders, and evidence-based progress."""
from functools import lru_cache
import gzip
import hashlib
import json
import os
from pathlib import Path
import tempfile
import time

HERE = Path(__file__).resolve().parent
SPEC = HERE / 'solar-release-v1.json'


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, record):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile('w', dir=path.parent, delete=False) as stream:
        json.dump(record, stream, indent=2, allow_nan=False)
        stream.write('\n')
    os.replace(stream.name, path)


def orders(root):
    local = root / 'work-orders.json'
    return read(local if local.exists() else HERE / 'work-orders.json')


def prerequisites(targets, jobs):
    seen = set()
    todo = list(targets)
    while todo:
        name = todo.pop()
        if name not in seen:
            seen.add(name)
            todo.extend(jobs[name].get('after', []))
    return seen


def decision(root, order):
    """A dated judgement only applies to the exact question, release and evidence reviewed."""
    path = root / 'decisions' / (order['id'] + '.json')
    if not path.exists():
        return 'open'
    record = read(path)
    if (record.get('release_sha256') != digest(SPEC)
            or record.get('order_sha256') != order_hash(order)
            or not record.get('reviewer') or not record.get('rationale')
            or not record.get('evidence')):
        return 'review stale'
    for name, expected in record['evidence'].items():
        source = root / name
        if not source.resolve().is_relative_to(root.resolve()) or not source.is_file() or digest(source) != expected:
            return 'review stale'
    value = record.get('decision')
    return value if value in ('adopted', 'sensitivity', 'rejected') else 'review stale'


def order_hash(order):
    return hashlib.sha256(json.dumps(order, sort_keys=True).encode()).hexdigest()


def routing(root, jobs):
    priorities, gates = {}, {}
    definitions = orders(root)
    preferred = root / 'shepherd' / 'recommendation.json'
    selected = None
    if preferred.exists():
        vote = read(preferred)
        state_path = root/'status.json'
        state = read(state_path).get('recovery', {}) if state_path.exists() else {}
        if (vote.get('brief_version') == 6 and vote.get('model') == 'gpt-6-luna' and vote.get('expires', 0) > time.time()
                and vote.get('release_sha256') == digest(SPEC) and state.get(vote.get('job')) == 'ready'):
            definitions.sort(key=lambda item: item['id'] != vote.get('order'))
            selected = vote.get('job')
    for index, order in enumerate(definitions):
        gates.update(order.get('gates', {}))
        # An explicit release handoff names the preconfigured restart job.
        # Only a current adoption can release that job; pilots use evidence gates.
        for name in order.get('restart_on_adoption', []):
            if decision(root, order) != 'adopted':
                gates[name] = None
        if decision(root, order) == 'open':
            for name in prerequisites(order['targets'], jobs):
                priorities.setdefault(name, index)
                if name == selected and index == 0:
                    priorities[name] = -1
    return priorities, gates


def thermal_turns(samples, tolerance_K=.1):
    """Resolved extrema in time order, suppressing sub-tolerance oscillations."""
    high=low=samples[0]
    direction=None
    events=[]
    for row in samples[1:]:
        temp=row['Teff_K']
        if direction!='up' and temp<low['Teff_K']:low=row
        if direction!='down' and temp>high['Teff_K']:high=row
        if direction!='up' and temp-low['Teff_K']>tolerance_K:
            if direction=='down':events.append(dict(kind='minimum',years=low['years'],Teff_K=low['Teff_K']))
            direction='up';high=row
        elif direction!='down' and high['Teff_K']-temp>tolerance_K:
            if direction=='up':events.append(dict(kind='maximum',years=high['years'],Teff_K=high['Teff_K']))
            direction='down';low=row
    return events


@lru_cache(maxsize=2)
def comparisons(path, size, mtime_ns):
    spec = read(SPEC)
    if digest(path) != spec['comparison_archive_sha256']:
        raise ValueError('boundary archive differs from release specification')
    rows = {}
    for track in json.loads(gzip.decompress(path.read_bytes()))['tracks']:
        key = track['key']
        if key in spec['comparisons']:
            if track['mass_Msun'] != spec['comparisons'][key] or track['Z_over_solar'] != 1:
                raise ValueError(f'wrong comparison identity: {key}')
            rows[key] = dict(mass_Msun=track['mass_Msun'], last=track['rows'][-1],
                             thermal_turns=thermal_turns(track['rows']),rewarming_fraction=track['rewarming_fraction'], scope=track['notes'])
    if rows.keys() != spec['comparisons'].keys():
        raise ValueError('incomplete boundary comparison')
    return rows


def is_column(job):
    from ember_campaign_recovery import column_job
    return column_job(job)


def board(root, status):
    import sys
    if str(HERE.parent) not in sys.path:
        sys.path.insert(0, str(HERE.parent))
    from birth_history import frontier
    from verify_endpoints import audit
    audit_rows = audit(root)
    states = status['recovery']
    jobs = {j['name']: j for j in expand_plan(root, read(root / 'recovery-plan.json'))}
    work = {}
    for order in orders(root):
        names = prerequisites(order['targets'], jobs)
        target_states = {n: states.get(n, 'unobserved') for n in order['targets']}
        disposition = decision(root, order)
        if disposition != 'open':
            next_action = disposition
        elif all(s == 'passed' for s in target_states.values()):
            next_action = 'adjudicate evidence'
        elif any(states.get(n) in ('ready','running') for n in names):
            next_action = 'advance available prerequisites'
        else:
            next_action = 'diagnose blocked prerequisites'
        work[order['id']] = dict(order, disposition=disposition, next_action=next_action,
            owners={n: 'recovery/'+n for n in sorted(names) if states.get(n) in ('failed','interrupted','invalid column record')},
            unfinished_dependencies={n: [k for k in jobs[n].get('after', []) if states.get(k) not in ('passed','interpolated')]
                                     for n in sorted(names) if states.get(n) in ('ready','waiting','running','dependency failed','failed','interrupted','await scientific assessment')},
            ready_columns=sorted(n for n in names if states.get(n)=='ready' and is_column(jobs[n])),
            order_sha256=order_hash(order), target_states=target_states,
            ready=sorted(n for n in names if states.get(n) == 'ready'),
            running=sorted(n for n in names if states.get(n) == 'running'),
            failed_prerequisites=sorted(n for n in names if states.get(n) in ('failed','interrupted','invalid column record')))
    tracks = {}
    for key, result in audit_rows.items():
        mass = round(float(key)*100)
        location = root / f'm{mass:02d}' / 'endpoint-review.json'
        record = read(location) if location.exists() else {}
        complete = ('method' in result and result.get('initial_contraction_pending') is None
                    and record.get('release_sha256') == digest(SPEC))
        try:
            head = frontier(root, mass)
        except (ValueError, OSError, IndexError, KeyError) as exc:
            head = {'observation_error': str(exc)}
        tracks[key] = dict(birth_history=head, dispatch=status['tracks'].get(key),
            endpoint=result, complete=complete, release_bound=record.get('release_sha256') == digest(SPEC),
            work_orders=[k for k,v in work.items() if mass in v['masses']])
    path = HERE.parents[1] / 'research/record/boundary/tracks.json.gz'
    stat = path.stat()
    dispatch_path = root/'dispatch-review.json'
    dispatch = read(dispatch_path) if dispatch_path.exists() else None
    return dict(warm_routes={k:v['superseded_by'] for k,v in jobs.items() if v.get('superseded_by')}, dispatch=dispatch, release=read(SPEC)['release'], release_sha256=digest(SPEC), initial_Z=status['initial_Z'],
                updated_unix=time.time(), scheduler_heartbeat=status['updated_unix'],
                hold=status['hold'], phase=status['phase'], work_orders=work, tracks=tracks,
                lanes=read(root/'campaign.json')['lanes'],
                occupied_columns=[n for n,s in states.items() if s=='running' and n in jobs and is_column(jobs[n])],
                occupied_native=[k for k,row in status['tracks'].items() if row['status'] in ('running','atmosphere running')],
                occupied_recovery=[n for n,s in states.items() if s=='running'],
                complete=sum(v['complete'] for v in tracks.values()), comparisons=comparisons(path, stat.st_size, stat.st_mtime_ns))


def expand_plan(root, specs):
    """Extend the known H/He crossing workflow only for the registered stellar masses."""
    jobs = {j['name']: j for j in specs}
    stellar = next((o for o in orders(root) if o['id']=='stellar-evolution'), {})
    for mass in stellar.get('masses', []):
        prefix = f'm{mass:02}'
        check = prefix+'-release-edge-review'
        if check in jobs:
            continue
        template = jobs['m12-hydrogen-edge-column']
        directory = root/'warm-marcs'/f'extended-start-{prefix}'/'hydrogen-edge-column'
        trial = directory.parent/'later-evolution-hydrogen-trial'
        column = prefix+'-hydrogen-edge-column'
        assessment = prefix+'-hydrogen-edge-assessment'
        import sys
        driver = [sys.executable,str(HERE/'composition_edge.py'),str(root),str(mass)]
        base = dict(trigger='review birth atmosphere')
        if column not in jobs:
            prep = prefix+'-release-edge-prepare'
            deps = [prefix+'-later-evolution-hydrogen-trial']
            if mass != 11:
                deps.append('m11-release-edge-review')
            jobs[prep] = dict(base,name=prep,after=deps,command=driver+['prepare'],wall_seconds=120,
                outputs=[str(directory/k) for k in ('spec.json','prediction.json','assess.py','solve.py','candidate-family.dat')])
            python = template['resume_column']['python']
            jobs[column] = dict(base,name=column,after=[prep],cwd=str(directory),
                command=[python,str(directory/'solve.py'),'anchor','spec.json','3600'],
                outputs=[str(directory/f'hhe_arm_anchor_solargrid_{prefix}_hydrogen_edge.json')],
                wall_seconds=template['wall_seconds'],resume_column=dict(template['resume_column'],directory=str(directory)))
            jobs[assessment] = dict(base,name=assessment,after=[column],command=[sys.executable,str(directory/'assess.py')],
                wall_seconds=60, outputs=[str(directory/'assessment.json')])
        if assessment not in jobs:
            continue
        jobs[check] = dict(base,name=check,after=[assessment,prefix+'-later-evolution-hydrogen-trial'],command=driver+['review'],
            wall_seconds=120, outputs=[str(trial/'review.json')])
    if 'm06-800-birth-depth10-column' in jobs:
        from cooling_work import jobs_for_cooling
        for job in jobs_for_cooling(root,jobs):
            jobs.setdefault(job['name'],job)
    from atmosphere_work import apply
    apply(root,jobs)
    return list(jobs.values())
