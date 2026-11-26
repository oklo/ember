"""Register the two measured intervals and the independent final surface check."""
from pathlib import Path
import sys


def jobs_for_cooling(root, existing):
    template = existing['m06-800-birth-depth10-column']
    driver = [sys.executable, str(Path(__file__).with_name('cooling_interval.py').resolve()), str(root)]
    parent = root/'m06/cooling-intervals'
    jobs = []

    def add(name, after, command, outputs, wall, **kwargs):
        jobs.append(dict(name=name, after=after, trigger=template['trigger'], command=command,
                         outputs=list(map(str,outputs)),wall_seconds=wall,**kwargs))
        return name

    freeze = add('m06-freeze-cooling-interval-forecast',['m06-register-800K-forecast'],
                 driver+['1','freeze'],[parent/'forecast.json'],360)
    previous = 'm06-depth-10-from-850-score-800'
    for index in (1,2,3):
        cmd = driver+[str(index)]
        stage = parent/f'band-{index}'
        prefix = f'm06-band-{index}-'
        deps = [previous,freeze]
        if index==1:
            deps.append('m06-850-birth-depth10-score')
        choose = add(prefix+'decide',deps,cmd+['decide'],[stage/'decision.json'],90)
        names = ['verify'] if index==3 else ['midpoint','endpoint']
        prepared = add(prefix+'prepare',[choose],cmd+['prepare'],
                       [stage/name/file for name in names for file in ('solve.py','spec.json')],90)
        columns = []
        for name in names:
            directory = stage/name
            policy = dict(directory=str(directory),python=template['command'][0],cpu_seconds=3600,max_restarts=2)
            columns.append(add(prefix+name,[prepared],
                [policy['python'],str(directory/'solve.py'),'anchor','spec.json',str(policy['cpu_seconds'])],
                [directory/f'hhe_arm_anchor_m06_band_{index}_{name}.json'],template['wall_seconds'],
                cwd=str(directory),resume_column=policy))
        if index==3:
            add(prefix+'surface-review',columns,cmd+['verify'],[stage/'surface-review.json'],90)
        else:
            family = add(prefix+'family',columns,cmd+['family'],[stage/'native/reader-check.json'],120)
            previous = add(prefix+'native',[family],cmd+['native'],[stage/'native/native-review.json'],18000)
    return jobs
