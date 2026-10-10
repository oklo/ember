"""First observed candidate mass distributions using private Ember tracks.

These are single-object-equivalent masses for unresolved catalogue sources,
not a completeness-corrected birth mass function. The Monte Carlo envelope is
a specified age/initial-state/cloud/photometry sensitivity exercise, not a
calibrated membership posterior. Every associated candidate remains eligible.
"""
from pathlib import Path
import json,sqlite3
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.special import ndtri
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

W=Path(__file__).resolve().parent
inventory=json.loads((W/'track_inventory.json').read_text())
objects=pd.read_csv(W/'catalogues/objects_for_mass_inference.csv',dtype={'gaia_dr3_id':str},low_memory=False)
tracks={}
for r in inventory:
    if r['status']!='complete' or not r.get('completed_target') or not r['cluster']:continue
    if r['mass_Msun']>.25:continue
    cloud='f2' if r['cloud']=='f2' else 'nc'
    name=r['case']
    # Prefer early initial states and complete cold boundaries. Short solver
    # controls and old failed cases are not isochrone inputs.
    score=sum(weight for tag,weight in [('grid_v1',10),('shallow',20),('v2',5),
        ('earlier_ramp',30),('earlier_bounded',40),('cool500',50),('hybrid',50),('primordial',60)] if tag in name)
    key=r['cluster'],cloud,r['mass_Msun']
    if key in tracks and tracks[key]['score']>=score:continue
    history=[json.loads(line) for line in (W/'tracks'/name/'evolution/history.jsonl').read_text().splitlines()]
    tracks[key]=dict(score=score,record=r,age=np.array([s['years'] for s in history]),
                    logL=np.log10([s['luminosity_Lsun'] for s in history]),
                    Teff=np.array([s['Teff_K'] for s in history]))

def relation(cluster,cloud,age,initial_scale):
    chosen={m:t for (c,k,m),t in tracks.items() if c==cluster and k==cloud}
    # Above 0.075 Msun use the warm, clear tables shared by both scenarios.
    for (c,k,m),t in tracks.items():
        if c==cluster and k=='nc' and (cloud=='nc' or m>.075):chosen.setdefault(m,t)
    mass=[];lum=[];used=[]
    for m,t in sorted(chosen.items()):
        elapsed=age-initial_scale*t['record']['initial_contraction_scale_years']
        if not t['age'][0]<=elapsed<=t['age'][-1]:continue
        mass.append(m);lum.append(float(np.interp(elapsed,t['age'],t['logL'])))
        used.append(t['record']['case'])
    mass=np.array(mass);lum=np.array(lum)
    assert len(mass)>=6,(cluster,cloud,len(mass))
    return mass,lum,used

def invert_all_branches(mass,lum,observed,rng):
    """Retain multiple masses at a deuterium-burning fold, never sort L alone.

    Conditional branch weights follow a flat prior in log mass via the inverse
    slope. The overall calculation remains a sensitivity exercise, not a
    jointly fitted cluster mass-function posterior.
    """
    x=np.log10(mass);sign=np.sign(np.diff(lum));assert np.all(sign!=0)
    changes=np.flatnonzero(sign[1:]!=sign[:-1])+1
    limits=np.concatenate([[0],changes,[len(mass)-1]])
    roots=[];weights=[]
    for a,b in zip(limits[:-1],limits[1:]):
        yy=lum[a:b+1];xx=x[a:b+1]
        if yy[-1]<yy[0]:yy=yy[::-1];xx=xx[::-1]
        inverse=PchipInterpolator(yy,xx,extrapolate=False)
        val=inverse(observed);weight=abs(inverse.derivative()(observed))
        weight=np.where(np.isfinite(val),np.maximum(weight,1e-8),0.)
        roots.append(val);weights.append(weight)
    roots=np.array(roots);weights=np.array(weights)
    total=weights.sum(axis=0);valid=total>0
    random=rng.random(len(observed))*total
    which=np.argmax(np.cumsum(weights,axis=0)>random,axis=0)
    result=np.full(len(observed),np.nan)
    result[valid]=10**roots[which[valid],np.flatnonzero(valid)]
    return result,(weights>0).sum(axis=0)>1

ages={'pleiades':dict(nominal=125e6,low=110e6,high=140e6),
      'hyades':dict(nominal=650e6,low=600e6,high=800e6)}
receipt=dict(scope=__doc__,age_scenarios_years=ages,
    initial_age_model='n=1.5 constant-Teff contraction estimate GM^2/(7 R L), nominal scale 1, sensitivity scale uniform 0--2; deuterium timing is not independently calibrated.',
    uncertainty_sampling='256 stratified scenarios: uniform age and initial-age scale, equal clear/cloud-clearing choices, 0.04 dex common model luminosity shift, individual measured and conversion errors.',
    membership='Union of published selections; no combined probability or contamination correction.',
    completeness='No correction for missing faint tail members or unobserved sky. Deep and wide source tags retained.',
    binaries='Unresolved sources treated as single-object-equivalent luminosities. No binary population correction.',
    extinction='Published dereddened Gaia M_G used where available; residual colour/extinction effects are not fully calibrated.',clusters={})
receipt['deuterium_ambiguity']='All branches of a nonmonotonic mass-luminosity relation are retained, with inverse-slope weights conditional on a flat log-mass prior. The remaining joins between primordial-D and post-D initial states are not validated; the faintest bins remain exploratory.'
N=256
bins=np.geomspace(.008,.25,13);width=np.diff(np.log10(bins))
outputs=[];plots={};isochrone_rows=[]
for cl in ('hyades','pleiades'):
    rng=np.random.default_rng(20261009+(cl=='pleiades'))
    strata=(np.arange(N)+.5)/N
    sampled_age=ages[cl]['low']+rng.permutation(strata)*(ages[cl]['high']-ages[cl]['low'])
    sampled_scale=2*rng.permutation(strata)
    sampled_model_shift=.04*ndtri(rng.permutation(strata))
    d=objects[objects.selected_by_any&(objects.cluster==cl)&objects.logL_observed.notna()].copy()
    d=d.reset_index(drop=True)
    luminosity=d.logL_observed.to_numpy()
    measurement=d.e_logL_stat.fillna(.10).clip(lower=.01).to_numpy()
    conversion=d.logL_conversion_systematic_dex.fillna(.04).to_numpy()
    error=np.hypot(measurement,conversion)
    draws=np.full((N,len(d)),np.nan);hist=[];tail_hist=[];ambiguous=np.zeros(len(d))
    inside=d.associations.str.contains('inside_tidal_radius',na=False)
    tail=(d.associations.str.contains('leading_tail|trailing_tail|outside_tidal_radius|greater_pleiades',na=False)&~inside).to_numpy()
    deep=d.source_names.str.contains('schneider|zapatero',na=False).to_numpy()
    used=set()
    for cloud in ('nc','f2'):
        mass,lum,names=relation(cl,cloud,ages[cl]['nominal'],1.)
        for m,l,name in zip(mass,lum,names):
            isochrone_rows.append(dict(cluster=cl,cloud=cloud,assumed_age_years=ages[cl]['nominal'],
                initial_age_scale=1,mass_Msun=m,logL_Lsun=l,case=name))
    for k in range(N):
        age=sampled_age[k];scale=sampled_scale[k]
        cloud='nc' if k%2==0 else 'f2';mass,lum,names=relation(cl,cloud,age,scale);used.update(names)
        observed=luminosity+rng.normal(size=len(d))*error-sampled_model_shift[k]
        # Interpolate only within the actual track mass/luminosity domain.
        draws[k],multiple=invert_all_branches(mass,lum,observed,rng)
        ambiguous+=multiple
        hist.append(np.histogram(draws[k],bins)[0]);tail_hist.append(np.histogram(draws[k,tail],bins)[0])
    valid=np.isfinite(draws);covered=valid.any(axis=0)
    q=np.full((3,len(d)),np.nan)
    q[:,covered]=np.nanquantile(draws[:,covered],[.16,.50,.84],axis=0)
    d['mass_sensitivity_low_Msun']=q[0];d['mass_median_Msun']=q[1];d['mass_sensitivity_high_Msun']=q[2]
    d['fraction_scenarios_in_mass_range']=valid.mean(axis=0)
    d['fraction_scenarios_with_multiple_mass_solutions']=ambiguous/N
    d['tail_or_associated_outside_core']=tail;d['deep_survey_candidate']=deep
    d['mass_interpretation']='single-object-equivalent; observed candidate; age/cloud/initial-state sensitivity'
    outputs.append(d)
    counts=np.quantile(hist,[.16,.50,.84],axis=0)
    tails=np.quantile(tail_hist,[.16,.50,.84],axis=0)
    plots[cl]=(counts,tails)
    receipt['clusters'][cl]=dict(selected_with_luminosity=len(d),objects_with_any_mass_support=int(covered.sum()),
        objects_with_majority_mass_support=int((valid.mean(axis=0)>.5).sum()),
        expected_candidates_in_range=float(np.mean(valid.sum(axis=1))),
        expected_tail_candidates_in_range=float(np.mean(valid[:,tail].sum(axis=1))),
        minimum_mass_model=float(min(m for c,k,m in tracks if c==cl)),
        bins_Msun=bins.tolist(),count_quantiles=counts.tolist(),tail_count_quantiles=tails.tolist(),
        included_runs=sorted(used))
    np.savez_compressed(W/f'{cl}_mass_sensitivity_draws.npz',mass_Msun=draws,object_id=d.object_id.to_numpy(dtype=str))

result=pd.concat(outputs,ignore_index=True)
result.to_csv(W/'catalogues/ember_candidate_masses.csv',index=False)
pd.DataFrame(isochrone_rows).to_csv(W/'cluster_isochrones_provisional.csv',index=False)
with sqlite3.connect(W/'cluster_catalogues.sqlite') as con:
    result.to_sql('ember_candidate_masses',con,if_exists='replace',index=False)
(W/'observed_mass_function_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
fig,axes=plt.subplots(1,2,figsize=(11.5,5),sharey=True)
for ax,cl,color in zip(axes,('hyades','pleiades'),('#267783','#b96535')):
    counts,tails=plots[cl];middle=np.sqrt(bins[:-1]*bins[1:])
    ax.stairs(np.where(counts[1]>0,counts[1]/width,np.nan),bins,baseline=None,color=color,lw=2,label='All identified candidates')
    ax.fill_between(middle,counts[0]/width,counts[2]/width,color=color,alpha=.20,label='Age, clouds and measurement sensitivity')
    ax.stairs(np.where(tails[1]>0,tails[1]/width,np.nan),bins,baseline=None,color=color,lw=1.4,linestyle='--',label='Tails / associations outside core')
    ax.set_xscale('log');ax.set_yscale('log');ax.set_xlim(.008,.25);ax.set_ylim(1,1e4)
    ax.set_xticks([.01,.02,.04,.08,.15,.25]);ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_xlabel(r'Ember mass ($M_\odot$)');ax.grid(axis='y',alpha=.13)
    ax.text(.04,.96,cl.capitalize(),transform=ax.transAxes,fontsize=16,va='top')
    text='[Fe/H] = +0.15; 600–800 Myr' if cl=='hyades' else '[Fe/H] = 0; 110–140 Myr'
    ax.text(.04,.87,text,transform=ax.transAxes,fontsize=10,va='top')
    if cl=='hyades':ax.set_ylabel(r'Observed candidates per dex in mass')
handles,labels=axes[1].get_legend_handles_labels()
fig.legend(handles,labels,loc='upper center',ncol=3,fontsize=9,frameon=False)
fig.text(.08,.037,'Provisional masses. Includes associated members and tails; unresolved sources are treated as single objects.',fontsize=9)
fig.text(.08,.010,'Survey completeness and contamination are not corrected. The lowest-mass bins still need initial-state checks.',fontsize=9)
fig.tight_layout(rect=(0,.065,1,.93))
fig.savefig(W/'observed_mass_functions.png');fig.savefig(W/'observed_mass_functions.pdf')
print(json.dumps({cl:{k:v for k,v in a.items() if k not in ('included_runs','bins_Msun','count_quantiles','tail_count_quantiles')}
    for cl,a in receipt['clusters'].items()},indent=2))
