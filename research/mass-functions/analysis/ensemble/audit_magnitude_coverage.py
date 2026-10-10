"""Compare observed magnitudes without requiring an assigned luminosity/mass.

This checks one selection lost when only objects with luminosities are counted.
It is not a new membership or bolometric-correction fit.
"""
from pathlib import Path
from contextlib import closing
import json,sqlite3
import numpy as np
import pandas as pd
W=Path(__file__).resolve().parent;C=W.parent/'cluster-mass-functions-oct9-v1'
d=pd.read_csv(C/'catalogues/objects_for_mass_inference.csv',low_memory=False)
with closing(sqlite3.connect(C/'cluster_catalogues.sqlite')) as con:
    ids=pd.read_sql("select distinct object_id from source_records where source='olivares2018_bouy2015' and selected=1",con).object_id
meta=json.loads((W/'runs/pleiades_r3_assort_s41/initial.json').read_text());r=np.pi/180
sep=np.arccos(np.clip(np.sin(d.dec_2016_deg*r)*np.sin(meta['dec']*r)+np.cos(d.dec_2016_deg*r)*np.cos(meta['dec']*r)*np.cos((d.ra_2016_deg-meta['ra'])*r),-1,1))/r
q=d[d.object_id.isin(ids)&(d.Kmag<17)&(sep<3)]
rows=[];edges=np.arange(8,17.01,.5)
for lo,hi in zip(edges[:-1],edges[1:]):
    z=q[(q.Kmag>=lo)&(q.Kmag<hi)]
    rows.append(dict(K_lo=lo,K_hi=hi,N=len(z),N_missing_luminosity=int(z.logL_observed.isna().sum())))
coverage=pd.DataFrame(rows);coverage.to_csv(W/'dance_luminosity_coverage.csv',index=False)
q[q.logL_observed.isna()].to_csv(W/'dance_missing_luminosities.csv',index=False)
m=pd.read_csv(W/'model_magnitude_draws.csv')
e=pd.read_csv(W/'structural_endpoints.csv')
valid=e[e.structural_screen_pass&~((e.binary_mode=='random')&(e.binary_fraction>0))]
means=m.groupby(['run','cloud','K_lo','K_hi']).N.mean().reset_index()
means=means.merge(valid[['run','expected_birth_ratio']],on='run',validate='many_to_one')
means.N*=.90/(1-.043)
comp=means.groupby(['expected_birth_ratio','K_lo','K_hi']).N.agg(['min','median','max']).reset_index()
comp=comp.merge(coverage,on=['K_lo','K_hi'],validate='many_to_one')
comp.to_csv(W/'dance_magnitude_comparison.csv',index=False)
receipt=dict(selected_core_objects=len(q),without_luminosity=int(q.logL_observed.isna().sum()),
             faint_range='13.5 <= Ks < 17',
             faint_objects=int(coverage[coverage.K_lo>=13.5].N.sum()),
             faint_without_luminosity=int(coverage[coverage.K_lo>=13.5].N_missing_luminosity.sum()),
             interpretation='Counts missing luminosities are included in this magnitude check, but do not receive invented masses. The empirical magnitude-luminosity conversion and single-cluster dynamics remain assumptions.')
(W/'dance_magnitude_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
with closing(sqlite3.connect(W/'cluster_followup.sqlite')) as con:
    coverage.to_sql('dance_luminosity_coverage',con,if_exists='replace',index=False)
    comp.to_sql('dance_magnitude_comparison',con,if_exists='replace',index=False)
print(json.dumps(receipt,indent=2));print(comp.to_string(index=False,float_format=lambda x:f'{x:.4g}'))
