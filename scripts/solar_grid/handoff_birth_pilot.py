"""Return the independently checked pilot to its own native mass queue."""
import json,math,shlex,shutil,tempfile
from pathlib import Path
from coordinate import root,read,write,dirs,last,solar,solver

edge=root/'birth-atmosphere-2850'
def handoff():
    review=read(edge/'reviewed-column.json')
    native=Path(review['directory'])
    assert review['passes'] and read(native/'spec.json')['arms']['anchor']['Z']==solar
    assert read(native/'prediction.json')==read(edge/'prediction.json')==review['prediction']
    e=edge/'pilot';r=read(e/'report.json')
    assert r['stop_reason']=='accepted-step budget' and r['accepted_this_invocation']==12
    attempts=[json.loads(l) for l in (e/'attempts.jsonl').read_text().splitlines()]
    assert sum(bool(a['accepted'] and a['audit_pass'] and a['converged']) for a in attempts)==12
    cfg=dict(shlex.split(l) for l in (edge/'trial.cfg').read_text().splitlines() if l.strip())
    assert float(cfg['mass_Msun'])==.09 and Path(cfg['atmosphere_pressure_matched'])==edge/'trial-family.dat'
    history=[json.loads(l) for l in (e/'history.jsonl').read_text().splitlines()]
    assert history[0]['years']==0 and math.isclose(history[0]['surface_Z'],solar,abs_tol=1e-14)
    assert r['age_years']==history[-1]['years'] and all((e/f).is_file() for f in ['final.checkpoint','final.checkpoint.lithium'])
    cases=dirs(9)
    for case in cases:
        if read(case/'launch.json')['action']=='validated birth pilot':return case
        assert not (case/'evolution/history.jsonl').exists()
    dest=cases[-1].parent/f's{len(cases):04}'
    staging=Path(tempfile.mkdtemp(dir=dest.parent,prefix='.validated-'))
    shutil.copyfile(edge/'trial.cfg',staging/'config.txt')
    (staging/'evolution').symlink_to(e,target_is_directory=True)
    write(staging/'launch.json',dict(action='validated birth pilot',mass_Msun=.09,initial_Z=solar,
        start_age_years=0,short_step=False,source=str(e),qualification='Conditional local atmosphere extension',
        column_review=str(edge/'reviewed-column.json')))
    staging.rename(dest)
    return dest
write(edge/'handoff.json',dict(case=str(handoff())))
