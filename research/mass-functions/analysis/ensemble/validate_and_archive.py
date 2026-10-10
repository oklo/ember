"""Check population accounting and save a compact, queryable final ledger."""
from pathlib import Path
from contextlib import closing
import argparse,hashlib,json,sqlite3,time
import numpy as np
import pandas as pd
W=Path(__file__).resolve().parent;V1=W.parent/'cluster-dissolution-oct9-v1'
parser=argparse.ArgumentParser()
parser.add_argument('--allow-incomplete',action='store_true',help='Validate finished endpoints and explicitly record pending runs.')
args=parser.parse_args()
e=pd.read_csv(W/'structural_endpoints.csv');checks={}
planned={p.parent.name for base in [V1,W] for p in (base/'runs').glob('*/initial.json') if p.parent.name.startswith(('hyades_r','pleiades_r'))}
assert len(planned)==17
assert set(e.run)<=planned
pending=sorted(planned-set(e.run))
assert args.allow_incomplete or not pending,('Unfinished planned runs',pending)
assert e.run.is_unique
assert np.allclose(e.time_Myr,np.where(e.cluster=='hyades',650,125),rtol=0,atol=1e-6)
assert np.all(e.total_mass_Msun<=e.birth_mass_Msun*(1+1e-8))
assert np.all((e.aperture_halfmass_pc>0)&(e.aperture_halfmass_pc<e.matched_aperture_pc))
checks['completed_full_age_models']=len(e)
checks['pending_runs']=pending
checks['campaign_complete']=not pending
checks['all_endpoint_ages_and_total_mass_checks_pass']=True
for which in ['model','observed']:
    counts=pd.read_csv(W/f'{which}_count_draws.csv')
    keys=(['run','cloud','faint_extension'] if which=='model' else ['cluster'])+['scenario','rep']
    assert not counts.duplicated(keys).any()
    for geometry in ['profile','angular_profile']:
        profiles=pd.read_csv(W/f'{which}_{geometry}_draws.csv')
        summed=profiles.groupby(keys+['population']).N.sum().unstack('population').reset_index()
        joined=counts.merge(summed,on=keys,validate='one_to_one')
        assert len(joined)==len(counts)
        assert (joined.N_BD==joined.BD).all() and (joined.N_star==joined.stars).all()
    checks[which+'_disjoint_bins_preserve_counts']=True
    mf=pd.read_csv(W/f'{which}_mass_function_draws.csv')
    for name,lo,hi in [('N_BD',.03,.08),('N_star',.1,.2)]:
        selected=mf[(mf.mass_lo_Msun>=lo)&(mf.mass_hi_Msun<=hi)]
        summed=selected.groupby(keys).N.sum().rename('mass_bin_sum').reset_index()
        joined=counts.merge(summed,on=keys,validate='one_to_one')
        assert len(joined)==len(counts)
        assert (joined[name]==joined.mass_bin_sum).all()
    checks[which+'_mass_bins_preserve_comparison_counts']=True
cost=[];initial=[];final=[]
for r in e.itertuples():
    base=W if r.workspace==W.name else V1;p=base/'runs'/r.run
    m=np.load(p/'initial.npz')['mass'];initial.append(pd.DataFrame(dict(run=r.run,id=np.arange(1,len(m)+1),mass_Msun=m)))
    f=pd.read_csv(W/f'{r.run}_descendants.csv');assert f.id.is_unique;assert len(f)==r.N
    assert np.all((f.id>=1)&(f.id<=len(m)))
    f.insert(0,'run',r.run);final.append(f)
    endpoint=p/'continuation/receipt.json' if (p/'continuation.json').exists() else p/'receipt.json'
    assert json.loads(endpoint.read_text())['exit_code']==0,(r.run,endpoint)
    for path in [p/'receipt.json',p/'continuation/receipt.json']:
        if path.exists():cost.append(dict(run=r.run,stage=str(path.parent.relative_to(p)),**json.loads(path.read_text())))
cost=pd.DataFrame(cost);cost.to_csv(W/'execution_cost.csv',index=False)
checks['retained_descendants_have_unique_initial_ids']=True
checks['cpu_hours_completed']=float(cost.cpu_seconds.sum()/3600)
checks['cpu_hours_new_completed']=float(cost[cost.run.isin(e[e.workspace==W.name].run)].cpu_seconds.sum()/3600)
checks['cost_excludes_pending_runs']=bool(pending)
checks['cost_scope']='Includes discarded original seed43 work before its saved-state continuation. Software installation and small force-split probes excluded; sum of job wall times is not elapsed campaign time.'
checks['files']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in W.glob('*.py')}
checks['endpoint_sha256']={r.run:hashlib.sha256(Path(r.snapshot).read_bytes()).hexdigest() for r in e.itertuples()}
checks['validation_unix']=time.time()
(W/'validation_receipt.json').write_text(json.dumps(checks,indent=2)+'\n')
with closing(sqlite3.connect(W/'cluster_followup.sqlite')) as con:
    pd.concat(initial,ignore_index=True).to_sql('initial_components',con,if_exists='replace',index=False)
    pd.concat(final,ignore_index=True).to_sql('final_descendants',con,if_exists='replace',index=False)
    cost.to_sql('execution_cost',con,if_exists='replace',index=False)
    for name in ['dance_binary_sensitivity','stellar_radial_checks']:
        pd.read_csv(W/f'{name}.csv').to_sql(name,con,if_exists='replace',index=False)
    assert con.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
print(json.dumps({k:v for k,v in checks.items() if k!='files'},indent=2))
