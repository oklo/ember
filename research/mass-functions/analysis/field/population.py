"""Forward model a local field population using Ember dwell times.

Masses from the observed catalogue are deliberately never read. Independent
fits use either the corrected component temperature function or the corrected
single-object luminosity function, never their product.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.special import logsumexp

W=Path(__file__).resolve().parent
BETA=np.linspace(-1.5,.5,101)  # p(present age) proportional to exp(beta*age/Gyr)
ALPHA=np.linspace(-1.5,2,141) # dN/dM proportional to M**(-alpha)
MASS=np.geomspace(.003,.09,401)
LOGM=np.log10(MASS)
DM=np.gradient(LOGM)
T_EDGES=np.arange(600.,2101.,150.)
MBOL_EDGES=np.arange(14.,20.1,1.) # colder bins exceed current table coverage
AGE_INTERVALS=1200
STEPHENS=np.polynomial.Polynomial([2265.9,347.82,-60.558,3.151,-.060481,.00024506])
L_TEMP=np.sort(STEPHENS(np.array([9.75,17.25])))

def inventory():
    rows=[]
    for d in sorted((W/'tracks').iterdir()):
        p=d/'evolution/history.jsonl'
        if not p.exists():continue
        try:h=pd.read_json(p,lines=True)
        except ValueError:continue
        if h.empty:continue
        s=json.loads((d/'scope.json').read_text());last=h.iloc[-1]
        # Recompute the approximate omitted contraction time from each run's
        # own seed; changing mass or metallicity changes that seed.
        s['initial_contraction_scale_years']=6.67430e-8*(s['mass_Msun']*1.98847e33)**2/(7*h.iloc[0].radius_Rsun*6.957e10*h.iloc[0].luminosity_Lsun*3.828e33)/31557600
        report=json.loads((d/'evolution/report.json').read_text()) if (d/'evolution/report.json').exists() else {}
        rows.append(dict(case=d.name,mass=s['mass_Msun'],cloud=s['cloud'],
            history=h,scope=s,age=float(last.years/1e9),Teff=float(last.Teff_K),
            Mbol=float(4.74-2.5*np.log10(last.luminosity_Lsun)),
            complete=bool(report.get('requested_age_reached')),report=report))
    return rows

def select(cloud,thin=False,max_mass=.091,metallicity=0.):
    candidates={}
    for r in inventory():
        if abs(r['scope'].get('metallicity',0.)-metallicity)>1e-6:continue
        if (r['cloud']!=cloud and not (r['mass']>.091 and r['cloud']=='nc')) or r['mass']>max_mass:continue
        if '_thin_' in r['case'] and not thin:continue
        # Older-than-cold-limit states cannot re-enter these selections for
        # monotonically cooling, sub-hydrogen-burning objects.
        if not r['complete'] and not (r['Teff']<505 and r['Mbol']>20) and r['age']<9:continue
        score=1000*int(r['complete'])+r['age']
        if '_v5' in r['case']:score+=20
        if '_v7' in r['case']:score+=30
        if '_v8' in r['case']:score+=35
        if '_thin_' in r['case']:score+=40
        m=r['mass']
        if m not in candidates or score>candidates[m][0]:candidates[m]=(score,r)
    return [candidates[m][1] for m in sorted(candidates)]

def age_integral(lo,hi,beta=BETA):
    lo=np.asarray(lo);hi=np.asarray(hi)
    b=np.asarray(beta).reshape((-1,)+(1,)*lo.ndim)
    # beta=0 is a uniform distribution in age, not in log(age).
    b_safe=np.where(abs(b)<1e-12,1.,b)
    return np.where(abs(b)<1e-12,hi-lo,(np.exp(b*hi)-np.exp(b*lo))/b_safe)

def dwell(h,offset,values,edges,age_max=10.):
    """Integrate age PDF over each saved segment's crossings of an observable bin.

    Interpolate observables linearly with log(age), with the contraction-age
    offset. D-burning loops remain separate segments, not forced monotonic.
    """
    age=h.years.to_numpy()/1e9+offset
    age=np.maximum(age,1e-8)
    logage=np.log(age);v=np.asarray(values)
    out=np.zeros((len(BETA),len(edges)-1))
    for i in range(len(age)-1):
        if age[i]>=age_max:break
        if not age[i+1]>age[i]:continue
        for j,(vlo,vhi) in enumerate(zip(edges[:-1],edges[1:])):
            if abs(v[i+1]-v[i])<1e-14:
                if not vlo<=v[i]<vhi:continue
                low,high=age[i],min(age[i+1],age_max)
            else:
                fraction=sorted([(vlo-v[i])/(v[i+1]-v[i]),(vhi-v[i])/(v[i+1]-v[i])])
                f0=max(0,fraction[0]);f1=min(1,fraction[1])
                if f1<=f0:continue
                low=np.exp(logage[i]+f0*(logage[i+1]-logage[i]))
                high=min(np.exp(logage[i]+f1*(logage[i+1]-logage[i])),age_max)
                if high<=low:continue
            out[:,j]+=age_integral(low,high)
    return out/age_integral(0.,10.)[:,None]

def dwell_kernel(cloud,age_scale=1.,young_age=.2,thin=False):
    tracks=select(cloud,thin=thin)
    masses=np.array([r['mass'] for r in tracks]);rows=[]
    for r in tracks:
        h=r['history'];offset=r['scope']['initial_contraction_scale_years']/1e9*age_scale
        kt=dwell(h,offset,np.log10(h.Teff_K),np.log10(T_EDGES))
        kl=dwell(h,offset,4.74-2.5*np.log10(h.luminosity_Lsun),MBOL_EDGES)
        ky=dwell(h,offset,np.log10(h.Teff_K),np.log10(L_TEMP),young_age)
        ka=dwell(h,offset,np.log10(h.Teff_K),np.log10(L_TEMP))
        rows.append(np.concatenate([kt,kl,ky,ka],axis=1))
    native=np.stack(rows)
    dense=np.stack([np.stack([np.interp(LOGM,np.log10(masses),native[:,ib,k],left=np.nan,right=np.nan)
          for k in range(native.shape[-1])],axis=1) for ib in range(len(BETA))])
    # Non-covered low masses are not implicitly filled with zero.
    supported=(MASS>=masses[0])&(MASS<=masses[-1])
    return dense,supported,tracks

def state_grid(cloud,age_scale=1.,thin=False,cold_factor=1.,ages=None,masses=None,metallicity=0.):
    tracks=select(cloud,thin=thin,max_mass=.251 if masses is not None and np.max(masses)>.091 else .091,metallicity=metallicity)
    tm=np.array([r['mass'] for r in tracks])
    if ages is None:ages=np.geomspace(.001,10.,1200)
    if masses is None:masses=MASS
    ll=[];lt=[];valid=[]
    for r in tracks:
        h=r['history'];offset=r['scope']['initial_contraction_scale_years']/1e9*age_scale
        age=np.maximum(h.years.to_numpy()/1e9+offset,1e-8)
        logt=np.interp(np.log(ages),np.log(age),np.log10(h.Teff_K))
        logl=np.interp(np.log(ages),np.log(age),np.log10(h.luminosity_Lsun))
        if not r['complete'] and r['Teff']<505:
            # Values below the table edge are bounds used solely to assess the
            # interpolation of a censored mass bracket. No cold model is claimed.
            logt[ages>age[-1]]+=np.log10(cold_factor)
            logl[ages>age[-1]]+=4*np.log10(cold_factor)
        elif not r['complete']:
            # A rare unfinished late-age row is used only over its computed
            # interval. At later ages interpolate between completed mass rows.
            logt[ages>age[-1]]=np.nan;logl[ages>age[-1]]=np.nan
        ll.append(logl);lt.append(logt);valid.append(ages>=age[0])
    ll=np.array(ll);lt=np.array(lt);valid=np.array(valid)
    out=[]
    for a in [ll,lt,valid.astype(float)]:
        out.append(np.array([np.interp(np.log10(masses),np.log10(tm[np.isfinite(a[:,j])]),a[np.isfinite(a[:,j]),j],left=np.nan,right=np.nan)
                    for j in range(len(ages))]).T)
    if metallicity!=0:
        # Interpolate the metallicity DIFFERENCE to the finer solar grid. This
        # preserves its resolved D-burning and hydrogen-burning transitions,
        # while reproducing each calculated non-solar mass node exactly.
        sf=state_grid(cloud,age_scale,thin,cold_factor,ages,masses,0.)
        sc=state_grid(cloud,age_scale,thin,cold_factor,ages,tm,0.)
        for k,coarse,fine in [(0,sc[0],sf[0]),(1,np.log10(sc[1]),np.log10(sf[1]))]:
            interp=np.array([np.interp(np.log10(masses),np.log10(tm),coarse[:,j],left=np.nan,right=np.nan)
                             for j in range(len(ages))]).T
            out[k]=fine+out[k]-interp
        out[2]*=sf[2]
    return out[0],10**out[1],out[2]>.999,tracks

def kernel(cloud,age_scale=1.,young_age=.2,thin=False,cold_factor=1.,temperature_shift=0,metallicity=None):
    if metallicity is None:
        pieces=[kernel(cloud,age_scale,young_age,thin,cold_factor,temperature_shift,z) for z in [-.3,0.,.3]]
        return sum(w*p[0] for w,p in zip([.235,.63,.135],pieces)),np.logical_and.reduce([p[1] for p in pieces]),[r for p in pieces for r in p[2]]
    ae=np.geomspace(.001,10.,AGE_INTERVALS+1);ages=np.sqrt(ae[:-1]*ae[1:])
    ll,t,valid,tracks=state_grid(cloud,age_scale,thin,cold_factor,ages,metallicity=metallicity)
    t=t+temperature_shift;mb=4.74-2.5*ll
    masks=[]
    for values,edges in [(t,T_EDGES),(mb,MBOL_EDGES)]:
        masks.extend([(values>=lo)&(values<hi)&valid for lo,hi in zip(edges[:-1],edges[1:])])
    ltype=(t>=L_TEMP[0])&(t<L_TEMP[1])&valid
    masks.extend([ltype&(ages[None,:]<young_age),ltype])
    aw=age_integral(ae[:-1],ae[1:])/age_integral(0.,10.)[:,None]
    dense=np.stack([aw@mask.T for mask in masks],axis=-1)
    tm=np.array([r['mass'] for r in tracks]);supported=(MASS>=tm.min())&(MASS<=tm.max())
    return dense,supported,tracks

def fit(cloud,observable='temperature',age_scale=1.,young_age=.2,thin=False,mass_min=None,temperature_shift=0,cold_factor=1.,metallicity=0.,binary_fraction=0.):
    dense,supported,tracks=kernel(cloud,age_scale,young_age,thin,cold_factor,temperature_shift,metallicity)
    if mass_min is not None:supported &= MASS>=mass_min
    covariance=None
    if observable=='temperature':
        data=pd.read_csv(W/'kirkpatrick2024_corrected_temperature_function.csv').query('fit_use').sort_values('Teff_low_K')
        y=data.density_pc3.to_numpy();em=ep=data.error_pc3.to_numpy();which=slice(0,10)
    elif observable=='luminosity':
        data=pd.read_csv(W/'best2024_corrected_single_luminosity_function.csv').query('14<=Mbol_low<20').sort_values('Mbol_low')
        y=data.density_pc3.to_numpy();em=data.error_low_pc3.to_numpy();ep=data.error_high_pc3.to_numpy();which=slice(10,16)
    elif observable=='component_luminosity':
        data=np.load(W/'component_luminosity_function.npz')
        draws=(1-binary_fraction/.05)*data['draws']+binary_fraction/.05*data['binary_draws']
        y=draws.mean(axis=0)[1:7];covariance=np.cov(draws,rowvar=False)[1:7,1:7]
        em=ep=np.sqrt(np.diag(covariance));which=slice(10,16)
    else:raise ValueError(observable)
    mw=(MASS[None,:]/.05)**(1-ALPHA[:,None])*DM[None,:]
    allpred=np.einsum('am,bmk->abk',mw[:,supported],dense[:,supported,:],optimize=True)
    pred=allpred[:,:,which]
    if covariance is None:
        err=(em+ep)/2
        for _ in range(5):
            numerator=np.sum(pred*y/err**2,axis=-1)
            denominator=np.sum(pred**2/err**2,axis=-1)
            amp=numerator/denominator
            err=np.where(amp[:,:,None]*pred<y,em,ep)
        chi=np.sum(((amp[:,:,None]*pred-y)/err)**2,axis=-1)
    else:
        inverse=np.linalg.inv(covariance)
        numerator=np.einsum('abi,ij,j->ab',pred,inverse,y)
        denominator=np.einsum('abi,ij,abj->ab',pred,inverse,pred)
        amp=numerator/denominator
        residual=amp[:,:,None]*pred-y
        chi=np.einsum('abi,ij,abj->ab',residual,inverse,residual)
    sigma_amp=1/np.sqrt(denominator)
    fy=np.divide(allpred[:,:,-2],allpred[:,:,-1],out=np.zeros_like(amp),where=allpred[:,:,-1]>0)
    chi_age=((fy-.092)/np.where(fy<.092,.026,.030))**2
    # Approximate analytical marginalization over a positive amplitude with
    # a locally log-uniform prior. Bins use published asymmetric errors.
    logp=-.5*(chi+chi_age)+np.log(sigma_amp/amp)
    prob=np.exp(logp-logsumexp(logp))
    flat_idx=np.argmax(logp);ia,ib=np.unravel_index(flat_idx,logp.shape)
    def quant(x,w):
        return np.interp([.16,.5,.84],np.cumsum(w)/w.sum(),x).tolist()
    rec=dict(cloud=cloud,observable=observable,age_scale=age_scale,young_age_Gyr=young_age,
        mass_min=float(MASS[supported][0]),mass_max=float(MASS[supported][-1]),metallicity=metallicity,binary_fraction=binary_fraction,
        alpha_quantiles=quant(ALPHA,prob.sum(axis=1)),beta_quantiles=quant(BETA,prob.sum(axis=0)),
        best=dict(alpha=float(ALPHA[ia]),beta=float(BETA[ib]),amplitude_pc3_dex=float(amp[ia,ib]),
                  chi2_distribution=float(chi[ia,ib]),chi2_young=float(chi_age[ia,ib]),
                  young_fraction=float(fy[ia,ib]),predicted= (amp[ia,ib]*pred[ia,ib]).tolist()),
        data=y.tolist(),error_low=em.tolist(),error_high=ep.tolist(),
        track_cases=[r['case'] for r in tracks],prior='uniform alpha [-1.5,2], beta [-1.5,.5], local log-uniform positive normalization',
        limits='Three-node local metallicity mixture when metallicity is null; otherwise the stated fixed metallicity. Young-fraction calibration and clouds are sensitivity cases. Not a reconstructed birth IMF.')
    return dict(record=rec,probability=prob,amplitude=amp,amplitude_error=sigma_amp,
                prediction=pred,young_fraction=fy)

if __name__=='__main__':
    result={}
    for cl in ['f2','nc']:
        for obs in ['temperature','luminosity']:
            x=fit(cl,obs);key=cl+'_'+obs;result[key]=x['record']
            np.savez_compressed(W/f'fit_{key}.npz',alpha=ALPHA,beta=BETA,**{k:v for k,v in x.items() if k!='record'})
            print(key,json.dumps(x['record']['best']),flush=True)
    (W/'population_fits.json').write_text(json.dumps(result,indent=2)+'\n')
