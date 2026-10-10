"""Empirical bolometric LF of K24 individual components with volume corrections.

Prefer measured SED integrals and published super-magnitudes to temperatures.
Temperatures specify the published completeness radius, not the Ember mass.
"""
from pathlib import Path
import re,json,sqlite3
import numpy as np
import pandas as pd
from numpy.polynomial import Polynomial
from scipy.special import ndtr,ndtri
W=Path(__file__).resolve().parent;C=W.parent/'cluster-mass-functions-oct9-v1'
rng=np.random.default_rng(991024)
d=pd.read_csv(W/'field_individual_observations.csv');s=pd.read_csv(C/'raw/sanghi2023_properties.csv');best=pd.read_csv(W/'best2024_objects.csv')
def key(v):
    if not isinstance(v,str) or v in ['Null','nan']:return ''
    v=v.replace('−','-').lstrip('- ').upper();v=re.sub(r'\s+','',v)
    v=re.sub(r'^2MASS[WI]?J','2MASSJ',v);v=re.sub(r'^WISE[AP]*J','WISEJ',v)
    return v
aliases=['DefaultName','CommonName','2MASS','WISE','SDSS','ULAS','PSO','CWISE','Gl/GJ','LHS','LSPM','SIPS','SCR']
aliases=[a for a in aliases if a in d]
index={}
for i,row in d.iterrows():
    for col in aliases:
        k=key(row[col])
        if k:index.setdefault(k,set()).add(i)
def match(values):
    ids=set()
    for v in values:ids|=index.get(key(v),set())
    return next(iter(ids)) if len(ids)==1 else None
matches=[]
for j,r in s.iterrows():
    if 'specblend' in str(r.flag) or r.lbol_flag or not np.isfinite(r.log_lbol_lsun):continue
    i=match([r['name'],r.name_simbadable,r.designation_2mass,r.designation_WISE])
    if i is None:continue
    p=d.at[i,'Parallax (mas)'];ep=d.at[i,'Parallax_unc (mas)'];sp=r.plx_formula
    if not p>0 or not sp>0:continue
    value=r.log_lbol_lsun+2*np.log10(sp/p)
    intrinsic=max(.02**2,r.log_lbol_lsun_err**2-(2/np.log(10)*r.plxerr_formula/sp)**2)
    error=np.sqrt(intrinsic+(2/np.log(10)*ep/p)**2)
    matches.append(dict(object_id=d.at[i,'object_id'],source='Sanghi2023',source_row=int(j),logL=value))
    d.loc[i,['logL','e_logL','luminosity_method']]=[value,error,'Sanghi2023_SED_measured_flux']
for j,r in best.iterrows():
    if 'B' in str(r['flags']):continue # unresolved combined light is not a component
    i=match([r['name']])
    if i is None or d.at[i,'luminosity_method']=='Sanghi2023_SED_measured_flux':continue
    p=d.at[i,'Parallax (mas)'];sp=r.parallax_mas
    if not p>0 or not sp>0:continue
    value=(4.74-r.Mbol)/2.5+2*np.log10(sp/p)
    error=max(.025,r.e_Mbol/2.5)
    d.loc[i,['logL','e_logL','luminosity_method']]=[value,error,'Best2024_bolometric_or_supermag']
    matches.append(dict(object_id=d.at[i,'object_id'],source='Best2024',source_row=int(j),logL=value))

# An empirical W2 correction fills objects without enough near-IR flux. The
# calibration uses measured bolometric fluxes; evolutionary quantities ignored.
mw=s.W2+5*np.log10(s.plx_formula)-10;mb=4.74-2.5*s.log_lbol_lsun;bc=mb-mw
good=(s.age_category=='Field')&~s.flag.fillna('').str.contains('specblend|lowg|over-L|v-red|lowZ')&~s.lbol_flag&(s.plx_formula/s.plxerr_formula>10)&(s.W2err<.1)&np.isfinite(bc)
x=mw[good].to_numpy();y=bc[good].to_numpy()
# Extend the measured calibration through late T/Y dwarfs using direct
# bolometric estimates in the matched 25-pc sample. Do not extrapolate W2.
wm=pd.to_numeric(d['W2 (mag)'],errors='coerce')+5*np.log10(d['Parallax (mas)'])-10
direct=d.luminosity_method.str.startswith(('Best2024','Sanghi2023_SED'))&d.single_system&np.isfinite(wm)&(wm>13)
x=np.r_[x,wm[direct]];y=np.r_[y,4.74-2.5*d.logL[direct]-wm[direct]]
keep=np.ones(len(x),bool)
for _ in range(5):
    poly=Polynomial.fit(x[keep],y[keep],4);res=y-poly(x);scatter=1.4826*np.median(abs(res[keep]-np.median(res[keep])));keep=abs(res-np.median(res[keep]))<3*scatter
domain=np.quantile(x[keep],[.005,.995])
pred=poly.convert();w2=pd.to_numeric(d['W2 (mag)'],errors='coerce');p=d['Parallax (mas)'];absolute=w2+5*np.log10(p)-10
use=d.logL.isna()&absolute.between(*domain)&(p>0)&~d.remnant
use &= ~d['WISE_contam?'].fillna('').str.contains('yes|true',case=False)
d.loc[use,'logL']=(4.74-absolute-pred(absolute))[use]/2.5
we=pd.to_numeric(d['W2_unc (mag)'],errors='coerce').fillna(.05) if 'W2_unc (mag)' in d else pd.Series(.05,index=d.index)
d.loc[use,'e_logL']=np.sqrt((scatter/2.5)**2+.04**2+((1+pred.deriv()(absolute))*we/2.5)[use]**2+(2/np.log(10)*d['Parallax_unc (mas)']/p)[use]**2)
d.loc[use,'luminosity_method']='Sanghi2023_empirical_W2_BC'
# Missing deblended luminosities: a separately flagged conditional estimate
# from measured field L versus catalogue T, with substantial intrinsic scatter.
# This is not a direct flux measurement and is excluded in a saved control LF.
tcat=pd.to_numeric(d['Teff (K)'],errors='coerce')
train=d.logL.notna()&tcat.between(450,2300)&d.single_system
tx=np.log10(tcat[train]);ty=d.logL[train]
tf=Polynomial.fit(tx,ty,3);tr=ty-tf(tx)
ts=max(.12,float(1.4826*np.median(abs(tr-np.median(tr)))))
impute=d.logL.isna()&tcat.between(450,2300)
d.loc[impute,'logL']=tf(np.log10(tcat[impute]));d.loc[impute,'e_logL']=max(.18,ts)
d.loc[impute,'luminosity_method']='temperature_conditioned_missing_flux'
d.to_csv(W/'field_luminosity_observations.csv',index=False)
pd.DataFrame(matches).to_csv(W/'luminosity_crossmatches.csv',index=False)

sp=d.SpTOpt_indx.fillna(d.SpTNIR_indx)
temp=pd.to_numeric(d['Teff (K)'],errors='coerce').to_numpy();et=pd.to_numeric(d['Teff_unc (K)'],errors='coerce').fillna(100).to_numpy()
group=np.where(np.isfinite(temp),100+np.floor(temp/300),np.floor((sp-40)/5));multiple=~d.single_system
groups={}
for a in sorted(set(group[np.isfinite(group)])):
    for m in [False,True]:
        mask=d.hard_20pc&~d.remnant&(group==a)&(multiple==m)
        if mask.sum():groups[(a,m)]=(int(d.logL[mask].notna().sum()),int(mask.sum()))
par=d['Parallax (mas)'].to_numpy();ep=d['Parallax_unc (mas)'].to_numpy();ep=np.where(np.isfinite(ep)&(ep>0),ep,.01*par)
logl=d.logL.to_numpy();el=d.e_logL.fillna(.08).to_numpy()
intrinsic=np.sqrt(np.maximum(.02**2,el**2-(2/np.log(10)*ep/par)**2))
finite=np.isfinite(logl)&(par/ep>=7)&~d.remnant.to_numpy()
edges=np.arange(13,21.1,1);draws=[];measured_draws=[];imputed_draws=[];binary_draws=[]
selection_sum=np.zeros(len(d));weights_sum=np.zeros(len(d))
for it in range(2048):
    # Importance-resampled Gaussian parallax likelihood, p(parallax)~parallax^-4
    # (local uniform volume). Finite input catalogue boundary remains a caveat.
    z=ndtri((np.arange(128)+.5)/128)
    if it==0:
        pdist=par[:,None]+ep[:,None]*z;wt=np.where(pdist>20,np.maximum(pdist,20)**-4,0)
        cdf=np.cumsum(wt,axis=1)/np.maximum(wt.sum(axis=1)[:,None],1e-300)
    idx=np.minimum((cdf<rng.random(len(d))[:,None]).sum(axis=1),127)
    pp=pdist[np.arange(len(d)),idx]
    tt=temp+et*rng.normal(size=len(d))
    limit=np.where(np.isfinite(tt)&(tt<600),15.,20.)
    valid=finite&(pp>=1000/limit)&(~np.isfinite(tt)|(tt>=450))
    correction=np.where(np.isfinite(tt)&(tt<1200),1.06,np.where(np.isfinite(tt)&(tt<1350),1.03,1.))
    comp=np.ones(len(d))
    for k,(n,total) in groups.items():
        if n:comp[(group==k[0])&(multiple==k[1])]=max(.1,rng.beta(n+.5,total-n+.5))
    luminosity=logl+2*np.log10(par/pp)+intrinsic*rng.normal(size=len(d))+rng.normal(0,.02)
    luminosity[impute]=tf(np.log10(np.maximum(tt[impute],200)))+rng.normal(0,max(.18,ts),impute.sum())+rng.normal(0,.05)
    mbol=4.74-2.5*luminosity
    weight=correction/(4*np.pi/3*limit**3)/comp
    valid &= np.isfinite(mbol)
    # Exp(1) bootstrap includes Poisson counting variance in the weighted LF.
    ww=weight*rng.exponential(size=len(d))
    draws.append(np.histogram(mbol[valid],edges,weights=ww[valid])[0])
    measured_draws.append(np.histogram(mbol[valid&~impute],edges,weights=ww[valid&~impute])[0])
    imputed_draws.append(np.histogram(mbol[valid&impute],edges,weights=ww[valid&impute])[0])
    # Deliberately strong residual-binary sensitivity: 5% of apparent singles
    # are equal-light pairs. Resolved/known multiples are not corrected twice.
    sf=d.single_system.to_numpy();f=.05*sf
    bd=np.histogram(mbol[valid],edges,weights=(ww*(1-f))[valid])[0]
    bd+=np.histogram((mbol+2.5*np.log10(2))[valid],edges,weights=(ww*2*f)[valid])[0]
    binary_draws.append(bd)
    selection_sum+=valid;weights_sum+=np.where(valid,weight,0)
draws=np.array(draws);mean=draws.mean(axis=0);cov=np.cov(draws,rowvar=False)
np.savez_compressed(W/'component_luminosity_function.npz',edges=edges,draws=draws,mean=mean,covariance=cov,measured_draws=measured_draws,imputed_draws=imputed_draws,binary_draws=binary_draws)
out=pd.DataFrame(dict(Mbol_low=edges[:-1],Mbol_high=edges[1:],density_pc3=mean,error_pc3=np.sqrt(np.diag(cov))))
out.to_csv(W/'component_luminosity_function.csv',index=False)
d['selection_probability_LF']=selection_sum/len(draws);d['mean_volume_weight_LF']=weights_sum/len(draws)
d.to_csv(W/'field_luminosity_observations.csv',index=False)
receipt=dict(individual_rows=len(d),direct_luminosity_matches=len(matches),
 imputed_luminosities=int(impute.sum()),imputation_scatter_dex=max(.18,ts),imputed_density_by_bin=np.mean(imputed_draws,axis=0).tolist(),
 W2_calibration=dict(count=int(keep.sum()),domain=domain.tolist(),scatter_mag=float(scatter),coefficients=pred.coef.tolist()),
 volume='20 pc for Teff >=600 K and warmer objects without a BD temperature; 15 pc for 450--600 K; colder excluded.',
 correction='K24 sky factors; Gaussian parallax likelihood with local-volume prior; missing photometry correction within spectral-type/multiplicity classes; measurement and Poisson bootstrap covariance.',
 luminosity_methods=d[d.hard_20pc].luminosity_method.value_counts().to_dict(),
 temperature_rows=int(np.isfinite(temp).sum()),temperature_rows_with_luminosity=int((np.isfinite(temp)&np.isfinite(logl)).sum()),
 remaining='Unknown binaries, age-metallicity correlations, and unknown objects outside the parent catalogue are not fully recoverable. No progenitor/birth-population correction.',
 input_photometry_groups={str(k):v for k,v in groups.items()},LF=out.to_dict('records'))
(W/'component_luminosity_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
with sqlite3.connect(W/'field_census.sqlite') as con:
    d.to_sql('luminosity_observations',con,if_exists='replace',index=False)
    out.to_sql('component_luminosity_function',con,if_exists='replace',index=False)
print(json.dumps(receipt,indent=2))
