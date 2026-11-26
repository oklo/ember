"""Prepare and assess a measured H/He edge using the existing campaign column pipeline."""
import argparse
import json
from pathlib import Path
import shlex
import shutil
import sys

from cold_restart_check import audited_steps, same_restart
from release import digest, read, save


def process(root, mass, action):
    run = root/'warm-marcs'/f'extended-start-m{mass:02}'
    trial = run/'later-evolution-hydrogen-trial'
    folder = run/'hydrogen-edge-column'
    name = f'm{mass:02}-hydrogen-edge-column'
    outcome = read(trial/'outcome.json')
    source = Path(outcome['source'])
    ck = trial/'evolution/final.checkpoint'
    history = ck.parent/'history.jsonl'
    rows = [json.loads(line) for line in history.read_bytes().splitlines()]
    trials = [json.loads(line) for line in (ck.parent/'attempts.jsonl').read_bytes().splitlines()]
    if (outcome['crossed_original_domain'] is not True or not source.resolve().is_relative_to(run.resolve())
            or digest(source/'final.checkpoint') != outcome['checkpoint_sha256']
            or not same_restart(json.loads((source/'history.jsonl').read_bytes().splitlines()[-1]), rows[0])
            or not audited_steps(rows, trials, read(ck.parent/'report.json'))
            or rows[-1]['years'] != outcome['age_yr']):
        raise ValueError('saved composition crossing failed ancestry or step checks')
    sys.path.insert(0,str(root))
    sys.path.insert(0,str(root/'warm-marcs/extended-start-m30/low-gravity-pilot'))
    from test_hydrogen_edge import surface_state, probe
    surface = surface_state(ck.read_bytes().splitlines(), rows[-1])
    if surface != outcome['final_surface']:
        raise ValueError('checkpoint composition differs from crossing evidence')
    temp, gravity, x, he3, z = surface
    surface[0] = round(temp)
    cfg = dict(shlex.split(line) for line in (run/'config.txt').read_text().splitlines() if line.strip())
    original = Path(cfg['atmosphere_pressure_matched'])
    family = trial/'family.dat'
    predicted = probe(family, surface)
    if probe(original,surface) is not None or not predicted:
        raise ValueError('retained table does not span the measured edge')
    if action == 'review':
        assessment = read(folder/'assessment.json')
        spec = read(folder/'spec.json')
        if (assessment['passes'] is not True or digest(Path(assessment['source'])) != assessment['receipt_sha256']
                or spec['checkpoint'] != str(ck) or spec['mass_Msun'] != mass/100
                or spec['teff'] != surface[0] or spec['logg'] != gravity
                or spec['arms']['anchor'] != dict(X=x,He3=he3,He4=1-x-he3-z,Z=z)
                or read(folder/'prediction.json') != predicted
                or digest(folder/'candidate-family.dat') != digest(family)):
            raise ValueError('independent assessment differs from saved boundary')
        save(trial/'review.json',dict(accepted=len(rows)-1,last=rows[-1],history_sha256=digest(history),
            atmosphere_receipt_sha256=assessment['receipt_sha256'],checkpoint_sha256=digest(ck),
            qualification='Crossing between the original edge and the measured endpoint only. Birth atmosphere remains conditional.'))
        return
    template = root/'warm-marcs/extended-start-m20/hydrogen-edge-column'
    folder.mkdir(exist_ok=True)
    (folder/'grid_hhe').mkdir(exist_ok=True)
    spec = read(template/'spec.json')
    spec.update(tag=f'_solargrid_m{mass:02}_hydrogen_edge',teff=surface[0],logg=gravity,
                mass_Msun=mass/100, years=outcome['age_yr'], checkpoint=str(ck), Z_over_X=z/x,
                purpose=f'Independent atmosphere at the measured {mass/100:.2f} solar-mass crossing.',
                arms={'anchor':dict(X=x,He3=he3,He4=1-x-he3-z,Z=z)})
    if (folder/'spec.json').exists() and read(folder/'spec.json') != spec:
        raise ValueError('existing column has a different scientific identity')
    save(folder/'spec.json',spec)
    save(folder/'prediction.json',predicted)
    shutil.copyfile(family,folder/'candidate-family.dat')
    shutil.copyfile(template/'solve.py',folder/'solve.py')
    assessor = (template/'assess.py').read_text()
    if assessor.count('m20-hydrogen-edge-column') != 1:
        raise ValueError('unexpected campaign assessor')
    (folder/'assess.py').write_text(assessor.replace('m20-hydrogen-edge-column',name))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign',type=Path)
    parser.add_argument('mass',type=int,choices=range(11,30))
    parser.add_argument('action',choices=['prepare','review'])
    args = parser.parse_args()
    process(args.campaign.resolve(),args.mass,args.action)
