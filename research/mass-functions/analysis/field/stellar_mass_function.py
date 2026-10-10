"""Nearby M-dwarf anchor: empirical luminosities mapped through Ember.

The substellar population is inferred separately by population.py. This file
counts individual-component rows only, corrects missing component photometry
within spectral/multiplicity classes, and propagates age and counting errors.
"""
from pathlib import Path
import json,sys
import numpy as np
import pandas as pd
from scipy.special import ndtr
import population as pop

W=Path(__file__).resolve().parent
RNG=np.random.default_rng(94102026)
BINS=np.array([.08,.10,.125,.16,.20,.25])
N=512

def run():
    obs=pd.read_csv(W/'field_luminosity_observations.csv')
    # Match missing-photometry corrections within coarse spectral type and
    # multiplicity, rather than assuming binary components have single-star
    # photometric coverage. No unknown component gets its parent's flux.
    sp=obs.SpTOpt_indx.fillna(obs.SpTNIR_indx)
    obs['sptype_group']=np.floor((sp-40)/2)
    obs['known_multiple']=~obs.single_system
    groups={}
    for (a,m),g in obs[obs.hard_20pc&~obs.remnant].groupby(['sptype_group','known_multiple']):
        if 0<=a<=4:groups[(a,m)]=(int(g.logL.notna().sum()),len(g))
    sel=obs[~obs.remnant & obs.logL.notna() & (obs.probability_within_20pc>1e-4)
            & obs.logL.between(-4.5,-1.6) & (obs.parallax_snr>=7)].copy()
    # Log-luminosity interpolation over the completed field-age grid. A small
    # extrapolation above .25 is an overflow buffer only; never a plotted track.
    one_age=np.geomspace(.015,10,120);mass=np.geomspace(.075,.25,241)
    grids=[pop.state_grid('nc',ages=one_age,masses=mass,metallicity=z) for z in [-.3,0,.3]]
    ll=np.concatenate([g[0] for g in grids],axis=1);valid=np.concatenate([g[2] for g in grids],axis=1)
    tracks=[r for g in grids for r in g[3]];ages=np.tile(one_age,3)
    if max(r['mass'] for r in tracks)<.25:raise RuntimeError('The stellar field grid is not complete to .25 Msun.')
    ms=mass[:,None]*np.ones((1,len(ages)))
    logmass=np.log10(mass)
    fit=np.load(W/'fit_nc_component_luminosity_mixed.npz')
    bprob=fit['probability'].sum(axis=0);beta=fit['beta']
    choice=RNG.choice(len(beta),N,p=bprob/bprob.sum())
    age_width=np.tile(np.gradient(one_age),3)*np.repeat([.235,.63,.135],len(one_age))
    results=[];single_results=[];reference=[];refs_single=[];inner=[];hard=[]
    estimates=np.full((N,len(sel)),np.nan);sub=sel.reset_index(drop=True)
    overflow=0
    for ir in range(N):
        target=sub.logL.to_numpy()+RNG.normal(size=len(sub))*sub.e_logL.fillna(.05).to_numpy()
        roots=np.full((len(sub),len(ages)),np.nan);jac=np.zeros_like(roots)
        for ia in range(len(ages)):
            ok=valid[:,ia]&np.isfinite(ll[:,ia]);x=ll[ok,ia];m=logmass[ok]
            if len(x)<2:continue
            if np.any(np.diff(x)<=0):raise RuntimeError('Nonmonotonic stellar mass-luminosity relation; preserve multiple roots.')
            # The overflow extrapolation prevents an artificial pile-up at the
            # upper grid edge. Results above .25 are excluded from all bins.
            slope=(x[-1]-x[-2])/(m[-1]-m[-2]);x=np.r_[x,x[-1]+slope*(np.log10(.35)-m[-1])];m=np.r_[m,np.log10(.35)]
            roots[:,ia]=10**np.interp(target,x,m,left=np.nan,right=np.nan)
            deriv=np.gradient(x,m)
            jac[:,ia]=np.interp(target,x,1/deriv,left=0,right=0)
        # dN/dM ~ M^(-0.3) is a broad local stellar prior, not the BD fit.
        weights=np.nan_to_num(roots**.7)*jac*(age_width*np.exp(beta[choice[ir]]*ages))[None,:]
        sums=weights.sum(axis=1);cdf=np.cumsum(weights,axis=1)/np.maximum(sums[:,None],1e-300)
        picks=(cdf<RNG.random(len(sub))[:,None]).sum(axis=1);picks=np.minimum(picks,len(ages)-1)
        draws=roots[np.arange(len(sub)),picks];draws[sums==0]=np.nan
        estimates[ir]=draws
        # Bayesian bootstrap supplies counting uncertainty. Within-volume
        # probability is retained, with an inner-volume control reported below.
        base=RNG.exponential(size=len(sub))
        completeness=np.ones(len(sub))
        for key,(available,total) in groups.items():
            mask=(sub.sptype_group==key[0])&(sub.known_multiple==key[1])
            completeness[mask]=max(.1,RNG.beta(available+.5,total-available+.5))
        weight=base*sub.probability_within_20pc.to_numpy()/completeness
        # Residual unresolved companions are not mistaken for the whole known
        # multiplicity fraction. A separate sensitivity budget is retained.
        f=RNG.uniform(0,.05)*sub.single_system.to_numpy()
        split=draws*2**(-jac[np.arange(len(sub)),picks])
        component_weight=weight*(1-f);extra_weight=weight*2*f
        density=(np.histogram(draws,BINS,weights=component_weight)[0]+np.histogram(split,BINS,weights=extra_weight)[0])/(4*np.pi/3*20**3)/np.diff(np.log10(BINS))
        results.append(density)
        single_weight=weight*sub.single_system.to_numpy()
        single_results.append(np.histogram(draws,BINS,weights=single_weight)[0]/(4*np.pi/3*20**3)/np.diff(np.log10(BINS)))
        ref=(draws>=.1)&(draws<.2)
        split_ref=(split>=.1)&(split<.2)
        reference.append((component_weight[ref].sum()+extra_weight[split_ref].sum())/(4*np.pi/3*20**3))
        refs_single.append(single_weight[ref].sum()/(4*np.pi/3*20**3))
        inside=sub.distance_pc.to_numpy()<15
        inner.append((base/completeness)[ref&inside].sum()/(4*np.pi/3*15**3))
        hard.append((base/completeness)[ref&sub.hard_20pc.to_numpy()].sum()/(4*np.pi/3*20**3))
    np.savez_compressed(W/'stellar_mass_draws.npz',bins=BINS,density=results,single_density=single_results,
        reference_density=reference,single_reference_density=refs_single,
        inner15_reference_density=inner,hard20_reference_density=hard,masses=estimates,
        object_id=sub.object_id.to_numpy(dtype='U'),beta=beta[choice])
    q=np.nanquantile(estimates,[.16,.5,.84],axis=0)
    for name,a in zip(['ember_mass_lo','ember_mass_median','ember_mass_hi'],q):sub[name]=a
    sub.to_csv(W/'stellar_ember_mass_estimates.csv',index=False)
    rec=dict(objects_with_luminosities=len(sub),reference_Msun=[.1,.2],
        reference_density_quantiles=np.quantile(reference,[.16,.5,.84]).tolist(),
        single_reference_density_quantiles=np.quantile(refs_single,[.16,.5,.84]).tolist(),
        inner15_reference_density_quantiles=np.quantile(inner,[.16,.5,.84]).tolist(),
        hard20_reference_density_quantiles=np.quantile(hard,[.16,.5,.84]).tolist(),
        photometry_completeness_groups={str(k):v for k,v in groups.items()},
        tracks=[r['case'] for r in tracks],
        age_assumption='Local age distribution inferred from brown-dwarf population; stellar prior dN/dM proportional to M^-0.3 with luminosity Jacobian.',
        metallicity_nodes=[-.3,0,.3],metallicity_weights=[.235,.63,.135],
        residual_binaries='Uniform 0--5% of apparent singles treated statistically as equal-light pairs; known components not counted twice.',
        caution='Metallicity and age-metallicity correlation remain approximate. Missing-photometry correction assumes exchangeability within spectral-type/multiplicity classes. The 15-pc check omits the small residual-binary adjustment.')
    (W/'stellar_mass_receipt.json').write_text(json.dumps(rec,indent=2)+'\n')
    print(json.dumps(rec,indent=2))
if __name__=='__main__':run()
