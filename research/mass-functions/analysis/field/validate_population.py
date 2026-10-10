"""Check quadrature, native metallicity nodes and the assumed lower mass cutoff."""
import json
from pathlib import Path
import numpy as np
import population as p

W=Path(__file__).resolve().parent
records=json.loads((W/'population_final_fits.json').read_text())
baseline=records['nc_component_luminosity_mixed']
result={}

# The population quadrature is independent of the saved evolutionary mesh.
p.MASS=np.geomspace(.003,.09,801)
p.LOGM=np.log10(p.MASS);p.DM=np.gradient(p.LOGM);p.AGE_INTERVALS=2400
refined=p.fit('nc','component_luminosity',metallicity=None)
r=refined['record']
result['double_mass_and_age_quadrature']={
    'baseline_alpha_quantiles':baseline['alpha_quantiles'],
    'refined_alpha_quantiles':r['alpha_quantiles'],
    'baseline_beta_quantiles':baseline['beta_quantiles'],
    'refined_beta_quantiles':r['beta_quantiles'],
    'largest_LF_fraction_change':float(np.max(np.abs(np.array(r['best']['predicted'])/baseline['best']['predicted']-1))),
    'refined_chi2_LF':r['best']['chi2_distribution']}
print('quadrature',result['double_mass_and_age_quadrature'],flush=True)

# Native non-solar tracks must be reproduced by interpolation of differences.
node_checks=[]
for z in [-.3,.3]:
    tracks=p.select('nc',metallicity=z,max_mass=.251)
    masses=np.array([t['mass'] for t in tracks]);ages=np.array([.05,.1,.5,1,5,9.])
    ll,teff,valid,_=p.state_grid('nc',ages=ages,masses=masses,metallicity=z)
    for i,t in enumerate(tracks):
        h=t['history'];ta=h.years.to_numpy()/1e9+t['scope']['initial_contraction_scale_years']/1e9
        covered=(ages>=ta[0])&(ages<=ta[-1])&valid[i]
        expected=np.interp(np.log(ages),np.log(ta),np.log10(h.luminosity_Lsun))
        if covered.any():node_checks.append(float(np.max(abs(ll[i,covered]-expected[covered]))))
result['native_metallicity_luminosity_max_error_dex']=max(node_checks)
assert max(node_checks)<1e-10

p.MASS=np.geomspace(.003,.09,401)
p.LOGM=np.log10(p.MASS);p.DM=np.gradient(p.LOGM);p.AGE_INTERVALS=1200
cut=p.fit('nc','component_luminosity',metallicity=None,mass_min=.014)
result['lower_mass_cutoff_0014']=cut['record']
np.savez_compressed(W/'fit_nc_component_cutoff0014.npz',alpha=p.ALPHA,beta=p.BETA,
    **{k:v for k,v in cut.items() if k!='record'})

lf=np.load(W/'component_luminosity_function.npz')
fraction=lf['imputed_draws'].mean(axis=0)/lf['draws'].mean(axis=0)
result['imputed_luminosity_fraction_per_fitted_bin']=fraction[1:7].tolist()
result['imputed_luminosity_fraction_fitted_total']=float(lf['imputed_draws'].mean(axis=0)[1:7].sum()/lf['draws'].mean(axis=0)[1:7].sum())
result['field_metallicity_nodes']=[-.3,0,.3]
result['field_metallicity_weights']=[.235,.63,.135]
assert np.all(np.isfinite(refined['probability']))
assert np.isclose(refined['probability'].sum(),1)

# Repair an obsolete description, without changing numerical results.
for key,rec in records.items():
    rec['limits']='Three-node local metallicity mixture when metallicity is null; otherwise the stated fixed metallicity. Young-fraction calibration and clouds are sensitivity cases. Not a reconstructed birth IMF.'
(W/'population_final_fits.json').write_text(json.dumps(records,indent=2)+'\n')
(W/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
