"""Export complete selected member lists alongside provisional Ember masses."""
from pathlib import Path
import json,sqlite3,shutil,ast
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

W=Path(__file__).resolve().parent
doc=W.parents[1]/'docs/research/independent'
objects=pd.read_csv(W/'catalogues/objects_for_mass_inference.csv',dtype={'gaia_dr3_id':str},low_memory=False)
masses=pd.read_csv(W/'catalogues/ember_candidate_masses.csv',dtype={'gaia_dr3_id':str},low_memory=False)
assert not masses.object_id.duplicated().any()
assert not masses.loc[masses.gaia_dr3_id.notna(),['cluster','gaia_dr3_id']].duplicated().any()
good=masses.mass_median_Msun.notna()
assert (masses.loc[good,'mass_sensitivity_low_Msun']>0).all()
assert (masses.loc[good,'mass_median_Msun']>=masses.loc[good,'mass_sensitivity_low_Msun']).all()
assert (masses.loc[good,'mass_median_Msun']<=masses.loc[good,'mass_sensitivity_high_Msun']).all()
assert (masses.loc[good,'mass_median_Msun']<=.25*(1+1e-12)).all()
columns=['object_id']+[c for c in masses.columns if c not in objects.columns]
merged=objects.merge(masses[columns],on='object_id',how='left',validate='one_to_one')
counts={}
for cl in ('hyades','pleiades'):
    selected=merged[(merged.cluster==cl)&merged.selected_by_any].copy()
    selected.to_csv(W/f'catalogues/{cl}_members_and_ember_masses.csv',index=False)
    counts[cl]=len(selected)

receipt=json.loads((W/'observed_mass_function_receipt.json').read_text());rows=[]
for cl,c in receipt['clusters'].items():
    for i,(lo,hi) in enumerate(zip(c['bins_Msun'][:-1],c['bins_Msun'][1:])):
        rows.append(dict(cluster=cl,mass_lower_Msun=lo,mass_upper_Msun=hi,
            candidates_low=c['count_quantiles'][0][i],candidates_median=c['count_quantiles'][1][i],
            candidates_high=c['count_quantiles'][2][i],tail_candidates_median=c['tail_count_quantiles'][1][i],
            dN_dlog10M=c['count_quantiles'][1][i]/np.log10(hi/lo)))
pd.DataFrame(rows).to_csv(W/'observed_mass_function_bins.csv',index=False)

# Scientific edge check: retain both allowed masses at a folded relation,
# preserve an analytic monotonic inverse, and never clamp out-of-range data.
tree=ast.parse((W/'infer_observed_mass_functions.py').read_text())
definition=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='invert_all_branches')
namespace=dict(np=np,PchipInterpolator=PchipInterpolator)
exec(compile(ast.Module(body=[definition],type_ignores=[]),'<branch check>','exec'),namespace)
invert=namespace['invert_all_branches'];rng=np.random.default_rng(118)
mass=np.array([.01,.02,.04]);lum=2*np.log10(mass)
answer,amb=invert(mass,lum,np.array([-3.6,-2.8,-9]),rng)
assert np.allclose(answer[:2],10**(np.array([-3.6,-2.8])/2),rtol=1e-12)
assert np.isnan(answer[-1]) and not amb.any()
answer,amb=invert(mass,np.array([-3.,-2.,-3.]),np.full(2000,-2.5),rng)
assert amb.all() and (answer<.02).sum()>800 and (answer>.02).sum()>800
with sqlite3.connect(W/'cluster_catalogues.sqlite') as con:
    integrity=con.execute('pragma integrity_check').fetchone()[0]
assert integrity=='ok'
checks=dict(database_integrity=integrity,rows_with_luminosity=len(masses),selected_catalogue_counts=counts,
    unique_object_ids=True,unique_nonnull_gaia_ids_per_cluster=True,ordered_mass_quantiles=True,
    analytic_mass_inverse=True,multiple_mass_branches_preserved=True,out_of_range_values_not_clamped=True)
(W/'mass_inference_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
for ext in ('png','pdf'):
    shutil.copy2(W/f'observed_mass_functions.{ext}',doc/f'CLUSTER_MASS_FUNCTIONS_20261009.{ext}')
print(json.dumps(checks,indent=2))
