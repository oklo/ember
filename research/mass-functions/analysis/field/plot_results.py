from pathlib import Path
import json,csv,shutil,sqlite3
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,FuncFormatter,NullFormatter
from scipy.stats import chi2
W=Path(__file__).resolve().parent;C=W.parent/'cluster-mass-functions-oct9-v1';DOC=W.parents[1]/'docs/research/independent'
RNG=np.random.default_rng(10092026)
MASS=np.geomspace(.01,.08,160)
records=json.loads((W/'population_final_fits.json').read_text())
stellar=np.load(W/'stellar_mass_draws.npz');BINS=stellar['bins'];ref=stellar['reference_density']

def sample(key,n=512):
    f=np.load(W/f'fit_{key}.npz');prob=f['probability'];idx=RNG.choice(prob.size,n,p=prob.ravel()/prob.sum())
    ia,ib=np.unravel_index(idx,prob.shape)
    amp=np.maximum(1e-6,f['amplitude'][ia,ib]+RNG.normal(size=n)*f['amplitude_error'][ia,ib])
    alpha=f['alpha'][ia];beta=f['beta'][ib]
    density=amp[:,None]*(MASS[None,:]/.05)**(1-alpha[:,None])
    return density,alpha,beta,amp

# A sensitivity ensemble, not a posterior probability assigned to different
# cloud/age-calibration assumptions. Each case carries its statistical draws.
keys=[]
for cloud in ['nc','f2']:
    keys.extend([cloud+'_component_luminosity_mixed']+
        [cloud+'_component_'+x for x in ['young100','young300','coldbound','binary05','age0','age2','solar']])
ensembles=[sample(k) for k in keys]
bd=np.concatenate([x[0] for x in ensembles]);alphas=np.concatenate([x[1] for x in ensembles]);betas=np.concatenate([x[2] for x in ensembles]);amps=np.concatenate([x[3] for x in ensembles])
refdraw=RNG.choice(ref,len(bd));normalized=bd/refdraw[:,None]
nominal=sample('nc_component_luminosity_mixed',4096)
nomref=RNG.choice(ref,len(nominal[0]));single=sample('nc_luminosity_mixed',4096)
def band(a):return np.quantile(a,[.16,.5,.84],axis=0)

def integral(alpha,amp,lo=.03,hi=.08):
    power=1-alpha
    return amp/np.log(10)*np.where(abs(power)<1e-8,np.log(hi/lo),((hi/.05)**power-(lo/.05)**power)/np.where(abs(power)<1e-8,1,power))
ratio=integral(alphas,amps)/refdraw;nomratio=integral(nominal[1],nominal[3])/nomref
missing_age=np.divide(np.exp(betas*10)-np.exp(betas*9.383),np.expm1(betas*10),
    out=np.full_like(betas,(10-9.383)/10),where=abs(betas)>=1e-8)
missing_bound=.135*integral(alphas,amps,.025,.065)*missing_age
receipt=dict(mass_range_Msun=[.01,.25],brown_dwarf_powerlaw='dN/dM proportional to M^(-alpha); normalization per cubic parsec per dex',
    component_population='Individual components in K24; temperatures choose published completeness radius; luminosities and Ember determine masses.',
    slope_nominal_quantiles=np.quantile(nominal[1],[.16,.5,.84]).tolist(),
    slope_sensitivity_quantiles=np.quantile(alphas,[.16,.5,.84]).tolist(),
    beta_nominal_quantiles=np.quantile(nominal[2],[.16,.5,.84]).tolist(),
    BD_003_008_to_stars_010_020_nominal=np.quantile(nomratio,[.16,.5,.84]).tolist(),
    BD_003_008_to_stars_010_020_sensitivity=np.quantile(ratio,[.16,.5,.84]).tolist(),
    reference_density_quantiles=np.quantile(ref,[.16,.5,.84]).tolist(),
    missing_late_metalrich_cell_density_bound=np.quantile(missing_bound,[.5,.84,.975]).tolist(),
    uncertainty='16th--84th percentiles across statistical draws and declared sensitivity cases. This is not a full systematic-error posterior.',
    sensitivity_cases=keys,limits=['Adopted 0.01 Msun lower population cutoff; no free low-mass cutoff inference.',
       'Field metallicity represented by three nodes, Y=.27, scaled GS98 mixture; no age-metallicity correlation.',
       'Approximate missing-component luminosities are flagged and propagated; single-object LF is a separate cross-check.',
       'Residual unresolved binaries represented by a 0--5% equal-light-pair sensitivity, not a measured individual classification.',
       'Temperature-function mismatch shows remaining atmosphere/calibration limitations.',
       'Cluster curves remain incomplete; local field comparison is not a measurement of cluster birth differences.'])

plt.rcParams.update({'font.size':13,'axes.labelsize':14,'legend.fontsize':11,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':220})
fig,axes=plt.subplots(1,2,figsize=(13.7,6.2))
blue='#284d78';teal='#267783';orange='#b96535'
for ax in axes:
    ax.set_xscale('log');ax.set_yscale('log');ax.set_xlim(.01,.25)
    ax.xaxis.set_major_locator(FixedLocator([.01,.02,.04,.08,.1,.2,.25]));ax.xaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:g}'));ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r'Mass ($M_\odot$)',labelpad=9);ax.grid(axis='y',which='major',alpha=.13)
ax=axes[0];lo,med,hi=band(bd);nmed=band(nominal[0])[1]
ax.fill_between(MASS,lo,hi,color=blue,alpha=.18,lw=0)
ax.plot(MASS,nmed,color=blue,lw=2.5,label='20 pc census: individual components')
sl,sm,sh=band(single[0]);ax.plot(MASS,sm,color='#758292',lw=1.6,ls='--',label='25 pc: single-object check')
q=band(stellar['density']);ax.fill_between(BINS,np.r_[q[0],q[0,-1]],np.r_[q[2],q[2,-1]],step='post',color=blue,alpha=.18,lw=0)
ax.stairs(q[1],BINS,color=blue,lw=2.5,baseline=None)
ax.set_ylabel(r'Objects pc$^{-3}$ per dex in mass');ax.set_ylim(.003,.11)
ax.legend(loc='upper left',frameon=False,handlelength=2.7)
ax.text(.04,.76,'Ember field mass function',transform=ax.transAxes,fontsize=13)
ax.text(.04,.70,'Brown-dwarf fit; stellar census above 0.08 '+r'$M_\odot$',transform=ax.transAxes,fontsize=10.5,color='#3e4145')

ax=axes[1];ax.axvspan(.1,.2,color='#6b6e73',alpha=.05)
lo,med,hi=band(normalized);ax.fill_between(MASS,lo,hi,color=blue,alpha=.16,lw=0)
ax.plot(MASS,band(nominal[0]/nomref[:,None])[1],color=blue,lw=2.6,label='Field: corrected')
star_normalized=stellar['density']/ref[:,None];q=band(star_normalized)
assert np.allclose(np.sum(star_normalized[:,1:4]*np.diff(np.log10(BINS))[1:4],axis=1),1,atol=1e-12)
ax.fill_between(BINS,np.r_[q[0],q[0,-1]],np.r_[q[2],q[2,-1]],step='post',color=blue,alpha=.16,lw=0)
ax.stairs(q[1],BINS,baseline=None,color=blue,lw=2.6)
cluster_ratios={}
cb=np.array([.008,.01,.014,.02,.028,.04,.056,.08,.1,np.sqrt(.02),.2,.25])
for cl,color,label in [('hyades',teal,'Hyades: observed'),('pleiades',orange,'Pleiades: observed')]:
    m=np.load(C/f'{cl}_mass_sensitivity_draws.npz')['mass_Msun'];nr=((m>=.1)&(m<.2)).sum(axis=1)
    h=np.array([np.histogram(x,cb)[0] for x in m])/nr[:,None]/np.diff(np.log10(cb))
    low,mean,high=band(h);ax.stairs(np.where(mean>0,mean,np.nan),cb,baseline=None,color=color,lw=1.8,label=label)
    ax.fill_between(cb,np.r_[np.maximum(low,.0001),max(low[-1],.0001)],np.r_[high,high[-1]],step='post',color=color,alpha=.09,lw=0)
    rr=((m>=.03)&(m<.08)).sum(axis=1)/nr;cluster_ratios[cl]=np.quantile(rr,[.16,.5,.84]).tolist()
ax.set_ylim(.025,9);ax.set_ylabel('Number per dex / number at 0.10–0.20 '+r'$M_\odot$')
ax.legend(loc='upper left',frameon=False,handlelength=2.7)
ax.text(.97,.055,'Same stellar normalization',ha='right',va='bottom',transform=ax.transAxes,fontsize=10.5,color='#414449')
fig.text(.075,.055,'Field shading includes statistical uncertainty and age, cloud, metallicity and multiplicity sensitivity.',fontsize=11)
fig.text(.075,.021,'Cluster curves include identified tails; missing members and unresolved binaries remain uncorrected.',fontsize=11)
fig.subplots_adjust(left=.075,right=.987,bottom=.18,top=.97,wspace=.26)
for ext in ['png','pdf']:
    out=W/f'field_mass_function.{ext}';fig.savefig(out);shutil.copy2(out,DOC/f'FIELD_MASS_FUNCTION_20261009.{ext}')
plt.close(fig)

# Show the observables that the population model must reproduce. The
# temperature mismatch is retained instead of selecting only successful plots.
fig,axs=plt.subplots(1,3,figsize=(13.8,4.7));xx=np.arange(14.5,20,1)
for ax,obs,label in zip(axs,['component_luminosity','luminosity','temperature'],['Individual components','Single-object check','Temperature check']):
    r=records['nc_'+obs+'_mixed'];x=np.arange(675,2101,150) if obs=='temperature' else xx
    ax.errorbar(x,1000*np.array(r['data']),yerr=1000*np.array([r['error_low'],r['error_high']]),fmt='o',ms=4.5,color='#24272b',capsize=2.5,label='Corrected census')
    for cl,col,lab in [('nc',blue,'Clear atmosphere'),('f2',orange,'Cloudy atmosphere')]:
        rr=records[cl+'_'+obs+'_mixed'];ax.plot(x,1000*np.array(rr['best']['predicted']),color=col,lw=2,label=lab)
    ax.set_xlabel(r'$T_{\rm eff}$ (K)' if obs=='temperature' else r'$M_{\rm bol}$')
    ax.set_ylabel(r'Objects ($10^{-3}$ pc$^{-3}$ per bin)');ax.set_ylim(bottom=0)
    ax.text(.04,.96,label,transform=ax.transAxes,va='top',fontsize=11.5)
    ax.grid(axis='y',alpha=.12)
    if obs=='temperature':ax.set_xticks([600,1000,1400,1800,2100])
axs[0].legend(loc='lower left',fontsize=9,frameon=False)
fig.text(.07,.035,'Fits are separate; the overlapping 20 pc and 25 pc samples are never counted twice in one likelihood.',fontsize=10.5)
fig.subplots_adjust(left=.07,right=.985,bottom=.18,top=.96,wspace=.32)
for ext in ['png','pdf']:
    out=W/f'field_mass_function_checks.{ext}';fig.savefig(out);shutil.copy2(out,DOC/f'FIELD_MASS_FUNCTION_CHECKS_20261009.{ext}')
plt.close(fig)
receipt['cluster_observed_ratios']=cluster_ratios
receipt['goodness_of_fit']={k:dict(chi2_LF=v['best']['chi2_distribution'],bins=len(v['data']),
     nominal_dof=len(v['data'])-3,nominal_tail_probability=float(chi2.sf(v['best']['chi2_distribution'],len(v['data'])-3)))
     for k,v in records.items() if k.endswith('_mixed')}
(W/'field_mass_function_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
rows=[]
for m,a,b,c,n in zip(MASS,*band(bd),nmed):rows.append(dict(mass_Msun=m,density_low=a,density_median=b,density_high=c,density_nominal=n,segment='brown_dwarf_population_fit'))
for m,a,b,c in zip(np.sqrt(BINS[:-1]*BINS[1:]),*band(stellar['density'])):rows.append(dict(mass_Msun=m,density_low=a,density_median=b,density_high=c,density_nominal=b,segment='stellar_census'))
pd.DataFrame(rows).to_csv(W/'field_mass_function.csv',index=False)
with sqlite3.connect(W/'field_census.sqlite') as con:
    pd.DataFrame(rows).to_sql('field_mass_function',con,if_exists='replace',index=False)
print(json.dumps({k:v for k,v in receipt.items() if k not in ['sensitivity_cases','limits','goodness_of_fit']},indent=2))
