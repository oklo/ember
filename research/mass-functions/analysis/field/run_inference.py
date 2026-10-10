from pathlib import Path
import json
import numpy as np
import population as p
W=Path(__file__).resolve().parent
records={};cases=[]
for cloud in ['nc','f2']:
    for obs in ['component_luminosity','luminosity','temperature']:
        cases.append((f'{cloud}_{obs}_mixed',dict(cloud=cloud,observable=obs,metallicity=None)))
    for label,kwargs in [('young100',dict(young_age=.1)),('young300',dict(young_age=.3)),
        ('coldbound',dict(cold_factor=.5)),('binary05',dict(binary_fraction=.05)),
        ('age0',dict(age_scale=0.)),('age2',dict(age_scale=2.)),('solar',dict(metallicity=0.))]:
        cases.append((f'{cloud}_component_{label}',dict(cloud=cloud,observable='component_luminosity',**{'metallicity':None,**kwargs})))
for key,kw in cases:
    x=p.fit(**kw);records[key]=x['record']
    np.savez_compressed(W/f'fit_{key}.npz',alpha=p.ALPHA,beta=p.BETA,**{k:v for k,v in x.items() if k!='record'})
    print(key,x['record']['alpha_quantiles'],x['record']['beta_quantiles'],x['record']['best']['chi2_distribution'],flush=True)
(W/'population_final_fits.json').write_text(json.dumps(records,indent=2)+'\n')
