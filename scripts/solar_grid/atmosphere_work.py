"""Plan atmosphere work from material coverage and paired native crossing tests."""
from pathlib import Path
import re
import sys

from release import digest, read

HERE=Path(__file__).resolve().parent


def apply(root,jobs):
    anchor=jobs.get('m06-800-birth-depth10-column')
    if anchor is None:return
    python=jobs['m06-depth-10-from-850-score-800']['command'][0]
    climate=anchor['command'][0]
    check=root/'atmosphere-validation'
    trigger=anchor['trigger']

    def add(name,after,cmd,output,wall,**opts):
        jobs[name]=dict(name=name,after=after,trigger=trigger,command=list(map(str,cmd)),outputs=[str(output)],
                        wall_seconds=wall,**opts)
        return name

    for mass in (2,5):
        key=f'm{mass:02}-measured-1000-endpoint'
        if key in jobs:
            jobs[key]['command']=[python,str(HERE/'measured_endpoint.py'),str(root),str(mass)]

    fixed=add('m06-850-fixed-profile-check',['m06-850-birth-depth10-score'],
              [climate,HERE/'fixed_profile.py',root],check/'early850/fixed-profile.json',1800,resource='atmosphere')
    dependencies=dict(m05=['m05-extended-below-1100K','m05-1050-score'],
                      m06=['m06-depth-10-from-900-native-850','m06-900-birth-depth10-column'],
                      early850=[fixed,'m06-depth-10-from-900-native-850'],band1=['m06-band-1-native'])
    for study,after in dependencies.items():
        prefix=f'{study}-boundary-replay-'
        folder=check/study
        cmd=[python,HERE/'boundary_replay.py',root,study]
        prep=add(prefix+'prepare',after,cmd+['prepare'],folder/'setup.json',120)
        arms=[add(prefix+arm,[prep],cmd+[arm],folder/arm/'finished.json',18000) for arm in ('measured','trial')]
        add(prefix+'score',arms,cmd+['score'],folder/'comparison.json',90)
    for study in ('m05','m06','band1'):
        dest=jobs['m06-band-2-decide'].setdefault('after_attempt',[])
        name=f'{study}-boundary-replay-score'
        if name not in dest:dest.append(name)
    jobs['m06-band-2-midpoint']['interpolate_when']=['m06/cooling-intervals/band-2/decision.json','interpolate_midpoint']

    baseline=root/'cool-m06-800-birth-depth10/solve.py'
    detached=root/'cool-m06-850/detached-region/solve.py'
    # The other audited driver differs only in its first numerical pressure guess.
    compatible={'2fa02c3dc345a94378a709d4f941c3d9d137bddfe22e00f0b6a6539e8c2b506e',
                'c29bcdfaf295266ef1633c5ef0b4021e0b1dddfb357be8babe8395b17e0753de'}
    for name,job in list(jobs.items()):
        policy=job.get('resume_column')
        if not policy:continue
        driver=Path(policy['directory'])/'solve.py'
        if (policy['python']==climate and driver.is_file() and baseline.is_file() and
                detached.is_file() and digest(driver) in compatible):
            policy.update(detached_solver=str(detached),retry_detached_convection=True)
        m=re.fullmatch(r'm(\d\d)-hydrogen-edge-column',name)
        if not m:continue
        mass=int(m[1])
        folder=Path(policy['directory'])
        work=folder/'warm-source'
        cmd=[python,HERE/'warm_boundary.py',root,str(mass)]
        prefix=f'm{mass:02}-warm-source-'
        deps=[d for d in job.get('after',[]) if d!=prefix+'triage']
        triage=add(prefix+'triage',deps,cmd+['triage'],work/'triage.json',120)
        job['after']=[*deps,triage]
        f=work/'triage.json'
        if not f.exists() or not (folder/'spec.json').exists():continue
        route=read(f)
        if not route['warm'] or route['spec_sha256']!=digest(folder/'spec.json'):continue
        assessment=f'm{mass:02}-hydrogen-edge-assessment'
        if assessment not in jobs:continue
        if (root/'recovery'/assessment/'result.json').exists() and (folder/'assessment.json').is_file() and read(folder/'assessment.json')['passes']:continue
        prep=add(prefix+'prepare',[triage],cmd+['prepare'],work/'route.json',120)
        ref=add(prefix+'reference',[prep],cmd+['reference'],work/'reference/boundary.json',21600,resource='atmosphere')
        opacity=add(prefix+'opacity',[prep],cmd+['opacity'],work/'opacity.json',43200,resource='atmosphere')
        target=add(prefix+'target',[opacity],cmd+['target'],work/'target/boundary.json',21600,resource='atmosphere')
        job['superseded_by']=target
        # Preserve the canonical assessment target used by the mass-specific evolution gate.
        jobs[assessment]=dict(jobs[assessment],after=[ref,target],command=list(map(str,cmd+['assess'])))
