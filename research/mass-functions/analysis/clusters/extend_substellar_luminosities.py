"""Add observed SED luminosities and empirical K-band corrections for faint members.

No published evolutionary mass, radius or temperature is used. Kinematic
distances assume membership and a common Hyades velocity, fitted to Gaia data.
They do not use luminosity. This is a provisional uncertainty-aware conversion.
"""
from pathlib import Path
import json, sqlite3, re, csv
import numpy as np
import pandas as pd
from numpy.polynomial import Polynomial
from astropy.coordinates import SkyCoord
from astropy import units as u

W=Path(__file__).resolve().parent
obj=pd.read_csv(W/'catalogues/objects_with_luminosity.csv',dtype={'gaia_dr3_id':str},low_memory=False)
ref=pd.read_csv(W/'raw/sanghi2023_properties.csv')
with sqlite3.connect(W/'cluster_catalogues.sqlite') as con:
    sources=pd.read_sql_query('select * from source_records',con)
receipt=dict(source='https://doi.org/10.5281/zenodo.10086810',scope=__doc__)

def basis(ra,dec):
    a=np.deg2rad(np.asarray(ra));d=np.deg2rad(np.asarray(dec))
    return np.stack([np.stack([-np.sin(a),np.cos(a),0*a],axis=-1),
        np.stack([-np.sin(d)*np.cos(a),-np.sin(d)*np.sin(a),np.cos(d)],axis=-1)],axis=-2)

core=obj[(obj.cluster=='hyades')&obj.selected_by_any&
    obj.associations.str.contains('inside_tidal_radius',na=False)&obj.parallax.between(15,30)&
    (obj.parallax_error/obj.parallax<.02)&(obj.gaia_ruwe<1.4)&obj.pmra.notna()&obj.pmdec.notna()]
mat=basis(core.ra_2016_deg,core.dec_2016_deg)
y=4.74047*core[['pmra','pmdec']].to_numpy()/core.parallax.to_numpy()[:,None]
keep=np.ones(len(core),bool)
for _ in range(6):
    velocity=np.linalg.lstsq(mat[keep].reshape(-1,3),y[keep].ravel(),rcond=None)[0]
    residual=np.linalg.norm(y-mat@velocity,axis=1);keep=residual<2
receipt['hyades_velocity_fit']=dict(input_stars=len(core),retained_stars=int(keep.sum()),
    velocity_ICRS_km_s=velocity.tolist(),median_tangential_residual_km_s=float(np.median(residual[keep])),
    rejection_limit_km_s=2,scope='Common velocity fitted to Gaia proper motions and parallaxes, conditional on core membership.')

obj['faint_pmra']=obj.pmra;obj['faint_pmdec']=obj.pmdec
obj['faint_e_pmra']=np.nan;obj['faint_e_pmdec']=np.nan
idx={v:i for i,v in enumerate(obj.object_id)}
with (W/'raw/hyades_schneider2024_table1.csv').open() as f:phot=list(csv.reader(f))[2:]
pm24={re.sub(r'\s+','',r[0].split(' a a ')[0].split(' b b ')[0].split(' c c ')[0].split(' d d ')[0].split(' e e ')[0]).replace('++','+'):
      [float(re.findall(r'[-+]?\d+\.\d+',v)[0]) for v in r[1:3]] for r in phot}
for _,r in sources.iterrows():
    raw=json.loads(r.original);i=idx[r.object_id]
    if r.source=='schneider2022':
        pairs=re.findall(r'([-+]?\d+\.\d+)±(\d+\.\d+)',raw['source_line'].replace('−','-'))[:2]
        if len(pairs)==2:
            obj.loc[i,['faint_pmra','faint_e_pmra','faint_pmdec','faint_e_pmdec']]=[float(x) for p in pairs for x in p]
    elif r.source=='schneider2024':
        matches=[(n,v) for n,v in pm24.items() if n.startswith(r.source_id)]
        assert len(matches)==1,(r.source_id,matches)
        obj.loc[i,['faint_pmra','faint_pmdec']]=matches[0][1]
        obj.loc[i,['e_Jmag','e_Kmag']]=[raw['e_Jmag'],raw['e_Kmag']]

obj['luminosity_distance_method']='catalogue_distance'
obj['e_luminosity_distance_pc']=obj.e_distance_pc
par=(obj.parallax>0)&(obj.parallax_error/obj.parallax<.2)
obj.loc[par,'luminosity_distance_method']='Gaia_parallax'
obj.loc[par,'e_luminosity_distance_pc']=1000*obj.loc[par,'parallax_error']/obj.loc[par,'parallax']**2
pred=basis(obj.ra_2016_deg,obj.dec_2016_deg)@velocity
mu=obj[['faint_pmra','faint_pmdec']].to_numpy()/1000
kin=np.sum(pred*mu,axis=1)/(4.74047*np.sum(mu**2,axis=1))
perp=np.linalg.norm(pred-4.74047*kin[:,None]*mu,axis=1)
obj['hyades_fitted_kinematic_distance_pc']=np.where(obj.cluster=='hyades',kin,np.nan)
obj['hyades_kinematic_perpendicular_km_s']=np.where(obj.cluster=='hyades',perp,np.nan)
faint_h=(obj.cluster=='hyades')&obj.source_names.str.contains('schneider',na=False)
use=faint_h&obj.distance_for_luminosity_pc.isna()&np.isfinite(kin)&(kin>0)&(perp<3)
obj.loc[use,'distance_for_luminosity_pc']=kin[use]
obj.loc[use,'luminosity_distance_method']='Gaia_core_velocity_conditional_membership'
# 20% is an explicit working dispersion/measurement allowance, not a quoted
# formal uncertainty. Retain comparisons with measured and published distances.
obj.loc[use,'e_luminosity_distance_pc']=.20*kin[use]
check=faint_h&obj.distance_for_luminosity_pc.notna()&~use&np.isfinite(kin)
ratio=kin[check]/obj.loc[check,'distance_for_luminosity_pc'].to_numpy()-1
receipt['hyades_distance_check']=dict(count=int(check.sum()),median_fractional_difference=float(np.median(ratio)),
    rms_fractional_difference=float(np.sqrt(np.mean(ratio**2))),new_kinematic_distances=int(use.sum()),
    assumed_fractional_uncertainty=.20,
    note='Some published candidate parallaxes disagree substantially with common motion; retain their measured distances and do not treat kinematic membership as confirmed.')
known=faint_h&obj.kinematic_distance_pc.notna()&np.isfinite(kin)
receipt['hyades_published_convergent_point_check']=dict(count=int(known.sum()),
    maximum_fractional_difference=float(np.max(abs(kin[known]/obj.loc[known,'kinematic_distance_pc']-1))))
obj['logL_conversion_systematic_dex']=0.04

# Match observed SEDs to catalogue positions. Correct to J2000 where proper
# motions exist. Require a unique candidate within 2 arcsec in both directions.
ra0=obj.ra_2016_deg-obj.faint_pmra.fillna(0)*16/(3.6e6*np.cos(np.deg2rad(obj.dec_2016_deg)))
dec0=obj.dec_2016_deg-obj.faint_pmdec.fillna(0)*16/3.6e6
cat=SkyCoord(ra0.to_numpy()*u.deg,dec0.to_numpy()*u.deg)
reference=SkyCoord(ref.ra_j2000_formula.to_numpy()*u.deg,ref.dec_j2000_formula.to_numpy()*u.deg)
near,sep,_=cat.match_to_catalog_sky(reference)
second=cat.match_to_catalog_sky(reference,nthneighbor=2)[1].arcsec
matches=[]
for i,j in enumerate(near):
    if sep.arcsec[i]>2 or second[i]<4:continue
    r=ref.iloc[j]
    if not np.isfinite(r.log_lbol_lsun):continue
    # Prefer independent observed distance; a new reference parallax can fill
    # a missing distance, while a spectral-type distance cannot.
    if not np.isfinite(obj.loc[i,'distance_for_luminosity_pc']) and r.dist_formula_source=='dist_plx_formula':
        obj.loc[i,'distance_for_luminosity_pc']=r.dist_formula
        obj.loc[i,'e_luminosity_distance_pc']=r.disterr_formula
        obj.loc[i,'luminosity_distance_method']='Sanghi2023_parallax'
    distance=obj.loc[i,'distance_for_luminosity_pc']
    if not np.isfinite(distance):continue
    obj.loc[i,'logL_observed']=r.log_lbol_lsun+2*np.log10(distance/r.dist_formula)
    obj.loc[i,'logL_method']='Sanghi2023_observed_SED_rescaled'
    ed=obj.loc[i,'e_luminosity_distance_pc'];extra=2/np.log(10)*ed/distance if np.isfinite(ed) else 0
    obj.loc[i,'e_logL_stat']=np.hypot(r.log_lbol_lsun_err,extra)
    obj.loc[i,'logL_conversion_systematic_dex']=.08 if r.lbol_flag else .04
    matches.append(dict(object_id=obj.loc[i,'object_id'],reference_name=r['name'],separation_arcsec=float(sep.arcsec[i]),flagged_SED=bool(r.lbol_flag)))
receipt['SED_matches']=matches

# Empirical K-band correction as a function of absolute K, from parallax
# calibrators. Use young and field relations separately and carry their
# difference as a model uncertainty. No extrapolation beyond fitted domains.
calibrations={}
for band,err in [('K_MKO','Kerr_MKO'),('Ks_2MASS','Kserr_2MASS')]:
    M=ref[band]-5*np.log10(ref.dist_formula/10)
    BC=4.74-2.5*ref.log_lbol_lsun-M
    quality=((ref.dist_formula_source=='dist_plx_formula')&~ref.lbol_flag&ref[err].between(0,.15)&
        ref.log_lbol_lsun_err.between(0,.15)&np.isfinite(BC)&M.between(8,18)&
        ~ref.flag.str.contains('blend|binary|lowZ|disk|over-L',case=False,na=False))
    field=ref.age_category.isin(['Field','FLD-G'])
    for young in (False,True):
        good=quality&(field if not young else ~field)
        x=M[good].to_numpy();y=BC[good].to_numpy();keep=np.ones(len(x),bool)
        assert len(x)>25,(band,young,len(x))
        for _ in range(4):
            fit=Polynomial.fit(x[keep],y[keep],3);res=y-fit(x)
            scatter=1.4826*np.median(abs(res[keep]-np.median(res[keep])))
            keep=abs(res-np.median(res[keep]))<max(3*scatter,.12)
        domain=np.quantile(x[keep],[.02,.98])
        calibrations[band,young]=(fit,domain,scatter)
        receipt[f'{band}_{"young" if young else "field"}']=dict(count=int(keep.sum()),
            coefficients=fit.convert().coef.tolist(),absolute_magnitude_domain=domain.tolist(),scatter_mag=float(scatter))

for i,r in obj.iterrows():
    if np.isfinite(r.logL_observed) or not np.isfinite(r.Kmag):continue
    is_mko=any(s in str(r.Kmag_source) for s in ('schneider','zapatero'))
    # Other catalogues' mixed K systems need explicit filter metadata first.
    if not is_mko:continue
    distance=r.distance_for_luminosity_pc
    if r.cluster=='pleiades' and not np.isfinite(distance):
        distance=136.2;obj.loc[i,'distance_for_luminosity_pc']=distance
        obj.loc[i,'e_luminosity_distance_pc']=5.
        obj.loc[i,'luminosity_distance_method']='Pleiades_member_distance_136p2pc_depth5pc'
    if not np.isfinite(distance):continue
    # Small near-IR extinction, explicit mean correction for Pleiades.
    AK=.012 if r.cluster=='pleiades' else 0.
    M=r.Kmag-5*np.log10(distance/10)-AK
    young=r.cluster=='pleiades';fit,domain,scatter=calibrations['K_MKO',young]
    if not domain[0]<=M<=domain[1]:continue
    obj.loc[i,'logL_observed']=(4.74-M-fit(M))/2.5
    obj.loc[i,'logL_method']='Sanghi2023_empirical_K_MKO_'+('young' if young else 'field')
    ed=obj.loc[i,'e_luminosity_distance_pc'];ed=.1*distance if not np.isfinite(ed) else ed
    ek=r.e_Kmag if np.isfinite(r.e_Kmag) else .05
    emag=np.hypot(ek,5/np.log(10)*ed/distance)
    obj.loc[i,'e_logL_stat']=np.hypot((1+fit.deriv()(M))*emag,scatter)/2.5
    other=calibrations['K_MKO',not young][0]
    obj.loc[i,'logL_conversion_systematic_dex']=max(.06,abs(fit(M)-other(M))/2.5)

# The published Calar SEDs used 133.5 pc. Rescale only when this input is the
# retained SED estimate, avoiding any duplicate adjustment.
use=obj.logL_method.eq('Zapatero2018_SED_published_distance')
obj.loc[use,'distance_for_luminosity_pc']=136.2
obj.loc[use,'e_luminosity_distance_pc']=5.
obj.loc[use,'luminosity_distance_method']='Pleiades_member_distance_136p2pc_depth5pc'
obj.loc[use,'logL_observed']+=2*np.log10(136.2/133.5)
obj.loc[use,'e_logL_stat']=np.hypot(obj.loc[use,'e_logL_stat'],2/np.log(10)*5/136.2)
receipt['coverage']={cl:obj[(obj.cluster==cl)&obj.selected_by_any].logL_method.value_counts().to_dict() for cl in ('hyades','pleiades')}
(W/'substellar_luminosity_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
obj.to_csv(W/'catalogues/objects_for_mass_inference.csv',index=False)
with sqlite3.connect(W/'cluster_catalogues.sqlite') as con:
    obj.to_sql('objects_for_mass_inference',con,if_exists='replace',index=False)
print(json.dumps({k:receipt[k] for k in ('hyades_velocity_fit','hyades_distance_check','coverage')},indent=2))
