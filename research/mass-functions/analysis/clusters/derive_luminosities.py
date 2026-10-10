"""Observed luminosities and an empirical Gaia bolometric-correction control.

Uses Cifuentes 2020 SED luminosities/photometry, not that catalogue's masses.
No cluster age or evolutionary model enters the bolometric correction.
The DR2/DR3 bandpass and gravity/metallicity differences are assessed against
the independent Pleiades DANCe SED sample; this is a provisional conversion.
"""
from pathlib import Path
import json,io,sqlite3
import numpy as np
import pandas as pd
from numpy.polynomial import Polynomial

W=Path(__file__).resolve().parent; C=W/'catalogues'
lines=[l for l in (W/'raw/cifuentes2020.tsv').read_text().splitlines() if l.strip() and not l.startswith('#')]
ref=pd.read_csv(io.StringIO('\n'.join([lines[0]]+lines[3:])),sep='\t',dtype=str)
ref=ref.map(lambda x:x.strip() if isinstance(x,str) else x)
for col in ['Gmag','BPmag','RPmag','Lbol','e_Lbol','Dist','e_Dist','e_Gmag','e_RPmag']:
    ref[col]=pd.to_numeric(ref[col],errors='coerce')
good=np.ones(len(ref),bool)
for flag in ['Multiple','Young','RUWE','Excess','q_Gmag','q_RPmag']:
    good &= ref[flag].eq('false').to_numpy()
ref['G_RP']=ref.Gmag-ref.RPmag
ref['BCG']=4.74-2.5*np.log10(ref.Lbol)-(ref.Gmag-5*np.log10(ref.Dist/10))
good &= (ref.Lbol>0)&(ref.e_Lbol/ref.Lbol<.2)&np.isfinite(ref.BCG)&ref.G_RP.between(.7,1.85)
x=ref.loc[good,'G_RP'].to_numpy();y=ref.loc[good,'BCG'].to_numpy()
keep=np.ones(len(x),bool)
for _ in range(5):
    poly=Polynomial.fit(x[keep],y[keep],4)
    res=y-poly(x);scale=1.4826*np.median(abs(res[keep]-np.median(res[keep])))
    keep=abs(res-np.median(res[keep]))<max(3*scale,.06)
lo,hi=np.quantile(x[keep],[.01,.99])
# Deterministic held-out test: no cluster members are used to tune this fit.
rng=np.random.default_rng(20910);test=rng.random(len(x))<.2
cv=Polynomial.fit(x[keep & ~test],y[keep & ~test],4)
cvres=y[keep & test]-cv(x[keep & test])
receipt=dict(source='https://doi.org/10.1051/0004-6361/202038295',
    fitted_stars=int(keep.sum()),degree=4,coefficient_order='increasing powers of G-RP',
    coefficients=poly.convert().coef.tolist(),colour_valid_min=float(lo),colour_valid_max=float(hi),
    residual_scatter_mag=float(scale),held_out_stars=len(cvres),
    held_out_rms_mag=float(np.sqrt(np.mean(cvres**2))),
    note='Empirical DR2 colour-to-bolometric correction applied provisionally to DR3; cluster crosscheck below. No evolutionary masses used.')
obj=pd.read_csv(C/'combined_objects.csv',dtype={'gaia_dr3_id':str},low_memory=False)
gaia=pd.read_csv(C/'gaia_dr3_photometry.csv',dtype={'source_id':str})
gaia=gaia.rename(columns={'ruwe':'gaia_ruwe'})
assert not gaia.source_id.isna().any()
assert not obj.loc[obj.gaia_dr3_id.notna(), ['cluster','gaia_dr3_id']].duplicated().any()
# Unmatched faint objects have no Gaia identifier; they must remain separate.
obj=obj.merge(gaia,left_on='gaia_dr3_id',right_on='source_id',how='left',validate='many_to_one')
obj['G_RP_DR3']=obj.phot_g_mean_mag-obj.phot_rp_mean_mag
obj['BP_RP_DR3']=obj.phot_bp_mean_mag-obj.phot_rp_mean_mag
# Existing dereddened absolute magnitudes take precedence. For other nearby
# candidates use the parallax and leave extinction uncertainty explicit.
obj['distance_for_luminosity_pc']=obj.distance_pc
use_gaia=(obj.parallax>0)&(obj.parallax_error/obj.parallax<.2)
obj.loc[use_gaia,'distance_for_luminosity_pc']=1000/obj.loc[use_gaia,'parallax']
mg=obj.MG_dereddened.copy()
missing=mg.isna()&obj.distance_for_luminosity_pc.notna()
mg.loc[missing]=obj.loc[missing,'phot_g_mean_mag']-5*np.log10(obj.loc[missing,'distance_for_luminosity_pc']/10)
obj['MG_used']=mg
obj['extinction_treatment']=np.where(obj.MG_dereddened.notna(),'published dereddened M_G','unreddened provisional; correction pending')
gr=obj.G_RP_DR3
valid=gr.between(lo,hi)&mg.notna()&(obj.phot_rp_mean_flux_over_error>10)
obj['logL_from_G_control']=np.where(valid,(4.74-mg-poly(gr))/2.5,np.nan)
obj['logL_observed']=np.nan;obj['logL_method']='unavailable'
obj.loc[valid,'logL_observed']=obj.loc[valid,'logL_from_G_control']
obj.loc[valid,'logL_method']='Cifuentes2020_Gaia_colour_control'
# SED-derived luminosities are preferred and rescaled to an individual
# distance where available; retain the 136.2-pc source distance otherwise.
sed=(obj.Lbol_Lsun_136p2pc>0)&(obj.cluster=='pleiades')
distance=obj.distance_for_luminosity_pc.fillna(136.2)
sedlog=np.log10(obj.Lbol_Lsun_136p2pc)+2*np.log10(distance/136.2)
obj.loc[sed,'logL_observed']=sedlog[sed]
obj.loc[sed,'logL_method']='Bouy2015_SED_rescaled_distance'
faint=obj.logL_SED_2018.notna()
obj.loc[faint,'logL_observed']=obj.loc[faint,'logL_SED_2018']
obj.loc[faint,'logL_method']='Zapatero2018_SED_published_distance'
obj['e_logL_stat']=np.nan
distance_error_logL=pd.Series(np.nan,index=obj.index)
distance_error_logL.loc[use_gaia]=2/np.log(10)*obj.loc[use_gaia,'parallax_error']/obj.loc[use_gaia,'parallax']
obj['e_logL_distance']=distance_error_logL
dg=poly.deriv()(gr)
errG=2.5/np.log(10)/obj.phot_g_mean_flux_over_error
errRP=2.5/np.log(10)/obj.phot_rp_mean_flux_over_error
# G enters both M_G and G-RP; propagate that covariance algebraically.
photo_variance=((1+dg)*errG)**2+(dg*errRP)**2
obj.loc[valid,'e_logL_stat']=np.sqrt((max(scale,float(np.std(cvres)))/2.5)**2+
    photo_variance[valid]/2.5**2+distance_error_logL[valid]**2)
sed_flux_error=obj.e_Lbol_Lsun_136p2pc/obj.Lbol_Lsun_136p2pc/np.log(10)
obj.loc[sed,'e_logL_stat']=np.sqrt(sed_flux_error[sed]**2+distance_error_logL[sed].fillna(0)**2)
obj['luminosity_distance_uncertainty_missing']=obj.logL_observed.notna()&distance_error_logL.isna()
obj.loc[faint,'e_logL_stat']=obj.loc[faint,'e_logL_SED_2018']
# A statistical error is not a complete uncertainty budget: model SED fill-in,
# age-dependent colours, extinction and system binarity are recorded separately.
cross=sed&valid&obj.selected_by_any&obj.logL_observed.between(-3.5,-1)
delta=obj.loc[cross,'logL_from_G_control']-obj.loc[cross,'logL_observed']
receipt['pleiades_independent_SED_check']=dict(objects=int(cross.sum()),
    median_logL_difference=float(delta.median()),robust_scatter_dex=float(1.4826*np.median(abs(delta-delta.median()))),
    p90_abs_difference_dex=float(np.percentile(abs(delta),90)))
receipt['coverage']={cl:obj[(obj.cluster==cl)&obj.selected_by_any].logL_method.value_counts().to_dict()
                     for cl in ['hyades','pleiades']}
receipt['uncertainty_limits']=['SED distance errors remain missing where no individual parallax is available.',
    'Published SED uncertainties and parallax errors are combined conservatively without a covariance model.',
    'Extinction, unresolved binaries, atmosphere-model fill-in and metallicity/gravity systematics require separate treatment.',
    'The empirical calibration uses observed Gaia colours; DR2/DR3 differences are assessed only through the cluster crosscheck.']
(W/'luminosity_calibration.json').write_text(json.dumps(receipt,indent=2)+'\n')
obj.to_csv(C/'objects_with_luminosity.csv',index=False)
ref['fit_selected']=False;ref.loc[ref.index[good][keep],'fit_selected']=True
ref.to_csv(C/'bolometric_calibrators.csv',index=False)
with sqlite3.connect(W/'cluster_catalogues.sqlite') as con:
    obj.to_sql('objects_with_luminosity',con,if_exists='replace',index=False)
print(json.dumps(receipt,indent=2))
