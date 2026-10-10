"""Observed field census and published completeness-corrected distributions.

Preserve all original rows in SQLite, but count only individual-component rows.
Published evolutionary masses are never inputs to the luminosity conversion.
"""
from pathlib import Path
import json,sqlite3,re
import numpy as np
import pandas as pd
from numpy.polynomial import Polynomial
from bs4 import BeautifulSoup

W=Path(__file__).resolve().parent; C=W.parent/'cluster-mass-functions-oct9-v1'
raw=pd.read_csv(W/'raw/Kirkpatrick2024_20pcCensus_Table4_Parameters.csv',dtype=str)
raw['row_id']=np.arange(len(raw))
d=raw[pd.to_numeric(raw['#CompsOnThisRow'])==1].copy().reset_index(drop=True)
assert len(d)==3582
for name in ['Parallax (mas)','Parallax_unc (mas)','G (mag)','G_unc (mag)',
             'G_RP (mag)','G_RP_unc (mag)','K (mag)','Kerr (mag)','Ks (mag)',
             'Kserr (mag)','RA (deg)','Dec (deg)','SpTOpt_indx','SpTNIR_indx']:
    d[name]=pd.to_numeric(d[name],errors='coerce')
d['object_id']='K24_'+d.row_id.astype(str)
d['single_system']=pd.to_numeric(d['#CompsInThisSystem'],errors='coerce').eq(1)
d['distance_pc']=1000/d['Parallax (mas)']
d['remnant']=d['AdoptedInitialMassNote'].str.startswith('wd',na=False)|d['SpecTypeOpt'].str.match(r'^D[A-Z]',na=False)
d['logL']=np.nan;d['e_logL']=np.nan;d['luminosity_method']='unavailable'

cal=json.loads((C/'luminosity_calibration.json').read_text());p=Polynomial(cal['coefficients'])
gr=d['G (mag)']-d['G_RP (mag)'];mg=d['G (mag)']+5*np.log10(d['Parallax (mas)'])-10
valid=gr.between(cal['colour_valid_min'],cal['colour_valid_max'])&(d['Parallax (mas)']>0)&~d.remnant
e_dist=2/np.log(10)*d['Parallax_unc (mas)']/d['Parallax (mas)']
d.loc[valid,'logL']=(4.74-mg-p(gr))[valid]/2.5
d.loc[valid,'e_logL']=np.sqrt(.04**2+e_dist[valid]**2+
    (((1+p.deriv()(gr))*d['G_unc (mag)'].fillna(.03))**2+
     (p.deriv()(gr)*d['G_RP_unc (mag)'].fillna(.03))**2)[valid]/2.5**2)
d.loc[valid,'luminosity_method']='Cifuentes2020_empirical_Gaia_BC'

# Empirical absolute-Ks bolometric correction for M dwarfs, including resolved
# components without usable Gaia colors. Fit measured SED luminosities only.
r=pd.read_csv(C/'catalogues/bolometric_calibrators.csv')
for name in ['Ksmag','Lbol','Dist','e_Ksmag']:r[name]=pd.to_numeric(r[name],errors='coerce')
mk=r.Ksmag-5*np.log10(r.Dist/10);bc=4.74-2.5*np.log10(r.Lbol)-mk
good=r.fit_selected&mk.between(4.5,10.5)&np.isfinite(bc)&(r.e_Ksmag<.10)
x=mk[good].to_numpy();y=bc[good].to_numpy();keep=np.ones(len(x),bool)
for _ in range(5):
    fit=Polynomial.fit(x[keep],y[keep],3);res=y-fit(x)
    scatter=1.4826*np.median(abs(res[keep]-np.median(res[keep])))
    keep=abs(res-np.median(res[keep]))<max(3*scatter,.12)
domain=np.quantile(x[keep],[.01,.99])
faint=json.loads((C/'substellar_luminosity_receipt.json').read_text())
calibrations=[('Ks (mag)','Kserr (mag)',fit,domain,scatter,'Cifuentes2020_empirical_Ks_BC')]
for band,err,key in [('K (mag)','Kerr (mag)','K_MKO_field'),('Ks (mag)','Kserr (mag)','Ks_2MASS_field')]:
    z=faint[key];calibrations.append((band,err,Polynomial(z['coefficients']),z['absolute_magnitude_domain'],z['scatter_mag'],'Sanghi2023_'+key))
for band,err,poly,domain_here,scat,label in calibrations:
    absolute=d[band]+5*np.log10(d['Parallax (mas)'])-10
    ok=d.logL.isna()&absolute.between(*domain_here)&(d['Parallax (mas)']>0)&~d.remnant
    if band.startswith('Ks'):ok &= ~d['2MASS_contam?'].fillna('').str.contains('yes|true',case=False)
    d.loc[ok,'logL']=(4.74-absolute-poly(absolute))[ok]/2.5
    magerr=np.hypot(d[err].fillna(.05),5/np.log(10)*d['Parallax_unc (mas)']/d['Parallax (mas)'])
    d.loc[ok,'e_logL']=np.sqrt(((1+poly.deriv()(absolute))*magerr)[ok]**2+scat**2)/2.5
    d.loc[ok,'e_logL']=np.hypot(d.loc[ok,'e_logL'],.04)
    d.loc[ok,'luminosity_method']=label

# Export Best's object table and corrected luminosity function independently.
# It is a single-object sample; never add it to K24 (the objects overlap).
soup=BeautifulSoup((W/'literature/best2024.html').read_text(),'html.parser')
for annotation in soup.find_all('annotation'):annotation.decompose()
best=[]
for row in soup.find(id='S2.T1').find_all('tr')[2:]:
    cells=[x.get_text(' ',strip=True) for x in row.find_all(['td','th'])]
    if len(cells)!=11:continue
    def pair(t):
        n=re.findall(r'[-+]?\d+(?:\.\d+)?',t.replace('−','-'))
        return (float(n[0]),float(n[1]) if len(n)>1 else np.nan) if n else (np.nan,np.nan)
    par,ep=pair(cells[1]);mb,emb=pair(cells[8])
    best.append(dict(name=cells[0],parallax_mas=par,e_parallax_mas=ep,
        spectral_type=cells[2],flags=cells[3],Mbol=mb,e_Mbol=emb,method=cells[9],
        photometric_distance='{' in cells[1],original=json.dumps(cells)))
best=pd.DataFrame(best);assert len(best)==504
best.to_csv(W/'best2024_objects.csv',index=False)
lf=[]
for row in soup.find(id='S3.T4').find_all('tr'):
    cells=[x.get_text(' ',strip=True) for x in row.find_all(['td','th'])]
    if len(cells)!=8 or '≤' not in cells[0]:continue
    lo=float(re.match(r'\d+\.\d+',cells[0])[0])
    def errors(t):
        nums=re.findall(r'\d+(?:\.\d+)?',t)
        return [float(nums[0]),float(nums[-1])] if len(nums)>1 else [float(nums[0])]*2
    density=float(cells[5])*1e-3;err=np.array(errors(cells[7]))*1e-3
    lf.append(dict(Mbol_low=lo,Mbol_high=lo+1,density_pc3=density,
        error_low_pc3=err[0],error_high_pc3=err[1],fit_use=14<=lo<21,
        original=json.dumps(cells)))
lf=pd.DataFrame(lf);assert len(lf)==11
lf.to_csv(W/'best2024_corrected_single_luminosity_function.csv',index=False)

# K24 Table 17: measured temperature distribution, including components.
# This is a separate inference/control, not a second likelihood for the same
# objects. Rows with lower limits are excluded from the main temperature fit.
temperature=[]
for lo,raw_count,correction,dlim,density,error in [
    (1950,27,1,20,.81,.19),(1800,16,1,20,.48,.16),(1650,24,1,20,.72,.19),
    (1500,26,1,20,.78,.19),(1350,35,1,20,1.04,.23),(1200,69,1.03,20,2.12,.30),
    (1050,45,1.03,20,1.38,.26),(900,50,1.06,20,1.58,.27),(750,63,1.06,20,1.99,.31),
    (600,98,1.06,20,3.10,.37),(450,56,1.06,15,4.20,.67)]:
    assert abs(raw_count*correction/(4*np.pi/3*dlim**3)-density*1e-3)<.006e-3
    temperature.append(dict(Teff_low_K=lo,Teff_high_K=lo+150,count=raw_count,
        sky_correction=correction,complete_radius_pc=dlim,density_pc3=density*1e-3,
        error_pc3=error*1e-3,fit_use=lo>=600))
pd.DataFrame(temperature).to_csv(W/'kirkpatrick2024_corrected_temperature_function.csv',index=False)

# Probabilistic 20-pc selection with a local uniform-volume distance prior.
# The finite source list cannot fully recover objects omitted just outside
# 20 pc; compare to a hard cut and an inner 15-pc sample in the inference.
from scipy.special import ndtr,ndtri
z=ndtri((np.arange(512)+.5)/512)
par=d['Parallax (mas)'].to_numpy();ep=d['Parallax_unc (mas)'].to_numpy()
ep=np.where(np.isfinite(ep)&(ep>0),ep,.01*np.abs(par))
draw=par[:,None]+ep[:,None]*z
weight=np.where(draw>20,np.maximum(draw,20)**-4,0)
d['probability_within_20pc']=np.divide((weight*(draw>=50)).sum(axis=1),weight.sum(axis=1),out=np.zeros(len(d)),where=weight.sum(axis=1)>0)
d['hard_20pc']=par>=50
d['parallax_snr']=par/ep
d['volume_20pc_pc3']=4*np.pi/3*20**3
d.to_csv(W/'field_individual_observations.csv',index=False)
with sqlite3.connect(W/'field_census.sqlite') as con:
    raw.to_sql('source_kirkpatrick2024',con,if_exists='replace',index=False)
    d.to_sql('individual_observations',con,if_exists='replace',index=False)
    best.to_sql('source_best2024',con,if_exists='replace',index=False)
    lf.to_sql('best_corrected_single_LF',con,if_exists='replace',index=False)
receipt=dict(component_rows=len(d),best_objects=len(best),
    Cifuentes_Ks_calibration=dict(count=int(keep.sum()),domain=domain.tolist(),scatter_mag=float(scatter),coefficients=fit.convert().coef.tolist()),
    hard_20pc_components=int(d.hard_20pc.sum()),luminosity_methods=d[d.hard_20pc].luminosity_method.value_counts().to_dict(),
    remnant_rows=int(d.remnant.sum()),
    population_warning='K24 individual components and Best singles must not be conflated or counted twice.',
    missing_luminosity_warning='A component without deblended photometry retains missing luminosity; no parent flux is assigned to every child.',
    young_single_L0_L7_fraction=dict(value=.092,error_low=.026,error_high=.030,age_limit_Myr=200,age_sensitivity_Myr=[100,300]),
    sources=['https://doi.org/10.3847/1538-4365/ad24e2','https://arxiv.org/abs/2401.09535'])
(W/'observations_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
