"""Save queryable model provenance and final fits alongside the observations."""
from pathlib import Path
import json,sqlite3
import pandas as pd
import population as p
W=Path(__file__).resolve().parent
selected={}
for cl in ['nc','f2']:
    for z in [-.3,0,.3]:
        for r in p.select(cl,max_mass=.251,metallicity=z):selected[r['case']]=r
tracks=[];states=[]
for name,r in sorted(selected.items()):
    launch=json.loads((W/'tracks'/name/'launch.json').read_text())
    tracks.append(dict(case=name,mass_Msun=r['mass'],cloud=r['cloud'],metallicity=r['scope'].get('metallicity',0),
        end_age_Gyr=r['age'],end_Teff_K=r['Teff'],reached_requested_age=r['complete'],
        initial_age_correction_years=r['scope']['initial_contraction_scale_years'],
        configuration=str(W/'tracks'/name/'configuration.txt'),
        history=str(W/'tracks'/name/'evolution/history.jsonl'),
        binary_sha256=launch.get('binary_sha256'),configuration_sha256=launch.get('configuration_sha256'),
        scope_json=json.dumps(r['scope'])))
    h=r['history'];keep=[c for c in ['years','Teff_K','luminosity_Lsun','radius_Rsun','X_center','Y3_center'] if c in h]
    df=h[keep].copy();df.insert(0,'case',name);df.insert(1,'saved_state',range(len(df)));states.append(df)
t=pd.DataFrame(tracks);s=pd.concat(states,ignore_index=True)
t.to_csv(W/'selected_tracks.csv',index=False)
fits=json.loads((W/'population_final_fits.json').read_text())
f=pd.DataFrame([dict(case=k,cloud=v['cloud'],observable=v['observable'],alpha_median=v['alpha_quantiles'][1],
    beta_median=v['beta_quantiles'][1],chi2_distribution=v['best']['chi2_distribution'],record_json=json.dumps(v)) for k,v in fits.items()])
with sqlite3.connect(W/'field_census.sqlite') as con:
    t.to_sql('ember_tracks',con,index=False,if_exists='replace')
    s.to_sql('ember_track_states',con,index=False,if_exists='replace')
    f.to_sql('population_fits',con,index=False,if_exists='replace')
    pd.read_csv(W/'stellar_ember_mass_estimates.csv').to_sql('stellar_mass_estimates',con,index=False,if_exists='replace')
    con.execute('CREATE INDEX IF NOT EXISTS ember_track_state_case ON ember_track_states("case")')
    assert con.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
    schema=con.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name").fetchall()
rec=dict(selected_tracks=len(t),saved_track_states=len(s),completed_to_10Gyr=int(t.reached_requested_age.sum()),
    cold_edge_tracks=int((~t.reached_requested_age&(t.end_Teff_K<505)).sum()),
    other_partial_tracks=t[~t.reached_requested_age&(t.end_Teff_K>=505)][['case','end_age_Gyr']].to_dict('records'),
    database_tables=[r[0] for r in schema],sqlite_integrity='ok')
(W/'catalogue_receipt.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
