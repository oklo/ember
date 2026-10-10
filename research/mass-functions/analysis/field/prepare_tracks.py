"""Independent field-age extensions of the cluster physics, without source edits."""
from pathlib import Path
import json, re

W=Path(__file__).resolve().parent
ROOT=W.parents[1]
C=ROOT/'out/cluster-mass-functions-oct9-v1'
inventory=json.loads((C/'track_inventory.json').read_text())
chosen={}
for r in inventory:
    if r.get('cluster')!='pleiades' or not r.get('completed_target'):continue
    cl='f2' if r['cloud']=='f2' else 'nc'
    key=(cl,r['mass_Msun']);name=r['case']
    score=sum(v for tag,v in [('grid_v1',10),('shallow',20),('v2',5),('earlier_ramp',30),
            ('earlier_bounded',40),('cool500',50),('hybrid',50),('primordial',60)] if tag in name)
    if key not in chosen or score>chosen[key][0]:chosen[key]=(score,r)
cases=[]
for (cloud,mass),(score,r) in sorted(chosen.items()):
    if cloud=='f2' and mass>.075:continue
    case=f'field_{cloud}_m{mass:.3f}_v1'
    dest=W/'tracks'/case;dest.mkdir(parents=True,exist_ok=True)
    src=C/'tracks'/r['case'];cfg=(src/'configuration.txt').read_text()
    def setkey(k,v):
        global cfg
        line=f'{k} "{v}"'
        cfg=re.sub(r'^'+re.escape(k)+r'\s+.*$',lambda _:line,cfg,flags=re.M) if re.search(r'^'+re.escape(k)+r'\s+',cfg,re.M) else cfg+line+'\n'
    setkey('maximum_step_years','3e8')
    # Keep the lower-gravity cool table for primordial-D starts.
    boundary=('sonora' if mass<=.014 else 'cluster')
    setkey('atmosphere_pressure_matched',C/f'physics/{boundary}_pleiades_{cloud}_cool500.dat')
    (dest/'configuration.txt').write_text(cfg)
    scope=json.loads((src/'scope.json').read_text()) if (src/'scope.json').exists() else {}
    scope.update(target_years=10e9,maximum_steps=1500,maximum_cpu_seconds=1200,
        cluster='field',mass_Msun=mass,cloud=cloud,parent_case=r['case'],
        initial_contraction_scale_years=r['initial_contraction_scale_years'],
        goal='Field population kernel to 10 Gyr; stop at a real table limit, never extrapolate an unfinished state.')
    (dest/'scope.json').write_text(json.dumps(scope,indent=2)+'\n');cases.append(case)
(W/'cases.json').write_text(json.dumps(cases,indent=2)+'\n')
print(len(cases),'fresh field-age cases')
