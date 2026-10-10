"""Small queryable ledger of the completed analysis; preserve raw runs separately."""
from pathlib import Path
import sqlite3,json,hashlib,subprocess
import pandas as pd,numpy as np
W=Path(__file__).resolve().parent
db=W/'cluster_dissolution.sqlite'
with sqlite3.connect(db) as con:
    for stem in ['dynamical_history','dynamical_endpoints','mock_survey_ratios','observed_aperture_ratios','observed_ratio_counting_uncertainty','observed_source_contributions']:
        pd.read_csv(W/(stem+'.csv')).to_sql(stem,con,if_exists='replace',index=False)
    inputs=[];components=[];ends=[];cost=[]
    for p in sorted((W/'runs').glob('*_r*_s*')):
        if not (p/'initial.json').exists():continue
        meta=json.loads((p/'initial.json').read_text());inputs.append(dict(run=p.name,configuration_json=json.dumps(meta)))
        z=np.load(p/'initial.npz');d=pd.DataFrame(dict(run=p.name,id=np.arange(1,len(z['mass'])+1),initial_mass_Msun=z['mass'],initial_x_pc=z['pos'][:,0],initial_y_pc=z['pos'][:,1],initial_z_pc=z['pos'][:,2]));components.append(d)
        if (p/'final_descendants.csv').exists():
            d=pd.read_csv(p/'final_descendants.csv');d.insert(0,'run',p.name);ends.append(d)
        for rec in [p/'receipt.json',p/'continuation/receipt.json']:
            if rec.exists():cost.append(dict(run=p.name,stage=str(rec.parent.relative_to(p)),**json.loads(rec.read_text())))
    pd.DataFrame(inputs).to_sql('run_inputs',con,if_exists='replace',index=False)
    pd.concat(components).to_sql('initial_components',con,if_exists='replace',index=False)
    if ends:pd.concat(ends).to_sql('final_descendants',con,if_exists='replace',index=False)
    pd.DataFrame(cost).to_sql('execution_cost',con,if_exists='replace',index=False)
    selection=[]
    for cl in ['hyades','pleiades']:
        d=pd.read_csv(W/f'{cl}_spatial_selection.csv')
        selection.append(d[['object_id','cluster','ra_2016_deg','dec_2016_deg','distance_for_luminosity_pc','source_names','angle_to_center_deg','distance_to_center_pc','p_mass_BD','p_mass_reference_star']])
    pd.concat(selection).to_sql('observed_spatial_selection',con,if_exists='replace',index=False)
    for table in ['initial_components','final_descendants']:
        con.execute(f'CREATE INDEX IF NOT EXISTS {table}_run_id ON {table}(run,id)')
    con.execute('CREATE INDEX IF NOT EXISTS observed_selection_id ON observed_spatial_selection(object_id)')
    scope={p.name:p.read_text() for p in W.glob('*scope.json')}
    pd.DataFrame([dict(name=k,content=v) for k,v in scope.items()]).to_sql('scope',con,if_exists='replace',index=False)

files=list(W.glob('*.py'))+[W/'install/bin/petar.omp.bse.galpy.d',W/'vendor_receipt.json',W/'build_command.json']
receipt=dict(files={str(p.relative_to(W)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},cpu_seconds=sum(r['cpu_seconds'] for r in cost),sum_job_wall_seconds=sum(r['wall_seconds'] for r in cost),completed_cluster_runs=len(ends),sqlite=str(db),scope='Evolution cost includes the deliberately stopped compact-run stage and its continuation, excludes software compilation and short validation probes. Sum of job wall times is not elapsed campaign wall time.')
(W/'analysis_receipt.json').write_text(json.dumps(receipt,indent=2))
with (W/'petar_m4_compatibility.patch').open('w') as f:
    subprocess.run(['git','diff','--','galpy-interface/galpy_interface.h'],cwd=W/'vendor/PeTar',stdout=f,check=True)
print(json.dumps({k:v for k,v in receipt.items() if k!='files'},indent=2))
