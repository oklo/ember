"""Mass-distribution check within common observed regions, retaining wide data."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullFormatter
W=Path(__file__).resolve().parent;D=W.parents[1]/'docs/research/independent'
obs=pd.read_csv(W/'observed_mass_function_draws.csv');models=pd.read_csv(W/'model_mass_function_draws.csv')
e=pd.read_csv(W/'structural_endpoints.csv').set_index('run')
rows=[];rng=np.random.default_rng(35190)
plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':190})
fig,axes=plt.subplots(1,2,figsize=(12,5.2),sharey=True)
for ax,cl,color,scenario in zip(axes,['hyades','pleiades'],['#b55225','#246d98'],['central_matched','ir_matched']):
    o=obs[(obs.cluster==cl)&(obs.scenario==scenario)]
    bins=o.groupby('bin')[['mass_lo_Msun','mass_hi_Msun']].first().sort_index()
    edges=np.r_[bins.mass_lo_Msun.to_numpy(),bins.mass_hi_Msun.iloc[-1]]
    mid=np.sqrt(edges[1:]*edges[:-1]);width=np.diff(np.log10(edges))
    ref=(bins.mass_lo_Msun>=.1)&(bins.mass_hi_Msun<=.2)
    a=o.pivot(index='rep',columns='bin',values='N').to_numpy()
    # Binwise counting fluctuations and the reference normalization are drawn
    # together; the normalization bins are not independent of themselves.
    ix=rng.integers(len(a),size=8192);draw=rng.gamma(a[ix]+.5)
    density=draw/draw[:,ref].sum(1)[:,None]/width
    qq=np.quantile(density,[.16,.5,.84],axis=0)
    ax.fill_between(mid,qq[0],qq[2],color='black',alpha=.12)
    ax.plot(mid,qq[1],marker='o',ls='-',lw=1.7,ms=4,color='black',label='Observed candidates',zorder=5)
    for birth,mc,label in [(.63,color,'Common input ratio 0.63'),(.55,'#785b8f','Lower field ratio 0.55')]:
        eligible=e[(e.cluster==cl)&(e.expected_birth_ratio==birth)&e.structural_screen_pass]
        eligible=eligible[~((eligible.binary_mode=='random')&(eligible.binary_fraction>0))]
        if len(eligible)==0:continue
        curves=[]
        for run in eligible.index:
            for cloud in ['nc','f2']:
                m=models[(models.run==run)&(models.cloud==cloud)&(models.scenario==scenario)&(models.faint_extension!='continued_slope')]
                values=m.pivot(index='rep',columns='bin',values='N').to_numpy()
                curve=values.mean(0)/values[:,ref].sum(1).mean()/width
                curves.append(curve)
                for k,n in enumerate(curve):rows.append(dict(run=run,cluster=cl,cloud=cloud,scenario=scenario,expected_birth_ratio=birth,mass_lo_Msun=edges[k],mass_hi_Msun=edges[k+1],normalized_per_dex=n))
        curves=np.array(curves)
        ax.fill_between(mid,curves.min(0),curves.max(0),color=mc,alpha=.16)
        ax.plot(mid,np.median(curves,axis=0),color=mc,lw=2,label=label)
    # The synthetic-light grid ends at 0.25 Msun. Do not interpret the last
    # bin, which lacks scatter into the sample from above that grid boundary.
    ax.set_xscale('log');ax.set_yscale('log');ax.set_xlim(.029,.2);ax.set_ylim(.3,12)
    ax.set_xticks([.03,.05,.08,.1,.2]);ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r'Ember source-equivalent mass ($M_\odot$)')
    ax.text(.03,.97,cl.capitalize()+(' · 10 pc' if cl=='hyades' else ' · DANCe central region'),transform=ax.transAxes,va='top',weight='bold')
    ax.legend(loc='lower left',fontsize=9.5,frameon=False)
axes[0].set_ylabel(r'Number per dex / $N(0.10$–$0.20\,M_\odot)$')
fig.subplots_adjust(left=.085,right=.98,bottom=.20,top=.96,wspace=.12)
fig.text(.085,.08,'Models shown pass the broad mass and size checks. Shading spans realizations and atmosphere choices.',fontsize=10)
fig.text(.085,.035,'Normalization tests shape only; absolute stellar counts and spatial profiles are checked separately.',fontsize=10)
stem=D/'CLUSTER_MASS_SHAPES_20261009';fig.savefig(stem.with_suffix('.png'));fig.savefig(stem.with_suffix('.pdf'))
pd.DataFrame(rows).to_csv(W/'normalized_mass_shapes.csv',index=False)
print(stem.with_suffix('.png'))
