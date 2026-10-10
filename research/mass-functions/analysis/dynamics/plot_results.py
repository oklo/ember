"""Compact figure separating catalogue coverage from cluster evolution."""
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

W=Path(__file__).resolve().parent;OUT=W.parents[1]/'docs/research/independent'
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.labelsize':12,'font.family':'DejaVu Sans','savefig.facecolor':'white'})
colors={'hyades':'#b45424','pleiades':'#266d9f'}
a=pd.read_csv(W/'observed_aperture_ratios.csv');obs=pd.read_csv(W/'observed_ratio_counting_uncertainty.csv')
ends=pd.read_csv(W/'dynamical_endpoints.csv');ends=ends[~ends.run.str.contains('_b25_')];mock=pd.read_csv(W/'mock_survey_ratios.csv')
hist=pd.read_csv(W/'dynamical_history.csv')
fig=plt.figure(figsize=(12.2,8.6));gs=fig.add_gridspec(2,2,height_ratios=[1,1.05],hspace=.36,wspace=.3)
ax=fig.add_subplot(gs[0,0]);bx=fig.add_subplot(gs[0,1]);cx=fig.add_subplot(gs[1,0]);dx=fig.add_subplot(gs[1,1])
for cl in colors:
    q=a[(a.cluster==cl)&(a.aperture=='radius_pc')&(a.limit>=10)&(a.limit<=200)].sort_values('limit')
    ax.plot(q.limit,q.ratio,'o-',color=colors[cl],lw=1.8,ms=4,label=cl.capitalize())
    ax.fill_between(q.limit,q.ratio_lo,q.ratio_hi,color=colors[cl],alpha=.12,lw=0)
ax.axhspan(.5481,.7148,color='#676767',alpha=.13,label='Field reference');ax.axhline(.6291,color='#555555',ls=':',lw=1)
ax.set_xscale('log');ax.set_xticks([10,20,30,60,200]);ax.xaxis.set_major_formatter(ScalarFormatter());ax.set_ylim(0,.8)
ax.set_xlabel('Catalogue aperture around cluster (pc)');ax.set_ylabel(r'$N(0.03–0.08\ M_\odot)\,/\,N(0.10–0.20\ M_\odot)$')
ax.legend(loc='upper right',frameon=False,fontsize=9);ax.text(.02,.97,'a  Observed catalogue',transform=ax.transAxes,va='top',fontweight='bold')
ax.text(.03,.04,'Larger apertures add mostly stellar detections.',transform=ax.transAxes,fontsize=9,color='#333333')

bx.axhspan(.5481,.7148,color='#676767',alpha=.13);bx.axhline(.6291,color='#555555',ls=':',lw=1)
positions={};i=0
for cl in ['pleiades','hyades']:
    q=ends[ends.cluster==cl].sort_values(['initial_rh_pc','binary_fraction','seed'])
    for _,row in q.iterrows():
        x=i;positions[row.run]=x
        m=mock[(mock.run==row.run)&mock.scenario.str.startswith('wide_union')]
        bx.plot(x,row.all_ratio,'o',ms=5,mfc='white',mec=colors[cl],label='All descendants' if i==0 else None)
        bx.plot(x,row.jacobi_ratio,'s',ms=5,mfc=colors[cl],mec=colors[cl],alpha=.55,label='Inside tidal radius' if i==0 else None)
        if len(m):
            lo=m.ratio.min();hi=m.ratio.max();bx.plot([x,x],[lo,hi],color=colors[cl],lw=6,solid_capstyle='round',label='Mock catalogue range' if i==0 else None)
        real=obs[(obs.cluster==cl)&(obs.aperture=='all')].iloc[0]
        bx.errorbar(x+.18,real.ratio,yerr=[[real.ratio-real.q16],[real.q84-real.ratio]],fmt='D',ms=4,color='black',capsize=2,label='Observed candidates' if i==0 else None)
        i+=1
    i+=.6
labels=[]
for name,x in positions.items():
    row=ends[ends.run==name].iloc[0]
    labels.append(f"{row.cluster.capitalize()}\n{row.initial_rh_pc:g} pc"+('\nbinaries' if row.binary_fraction else '\nsingle'))
bx.set_xticks(list(positions.values()),labels,fontsize=8);bx.set_ylim(0,1.1);bx.set_xlim(-.5,max(positions.values(),default=0)+.65);bx.set_ylabel('Brown-dwarf / reference-star ratio')
bx.text(.02,.97,'b  N-body predictions',transform=bx.transAxes,va='top',fontweight='bold');bx.legend(loc='upper left',bbox_to_anchor=(.01,.92),ncol=2,fontsize=8,frameon=False)

for cl in ['hyades','pleiades']:
    for run,q in hist[hist.cluster==cl].groupby('run'):
        if run not in positions:continue
        row=ends[ends.run==run].iloc[0];ls='--' if row.initial_rh_pc==2 else '-';lw=2 if row.binary_fraction else 1.3
        cx.plot(q.time_Myr,q.mass_within_jacobi_Msun,color=colors[cl],ls=ls,lw=lw,alpha=.85)
cx.errorbar(650,290,yerr=[[290-240],[330-290]],fmt='D',color=colors['hyades'],capsize=4,ms=5)
cx.errorbar(125,750,yerr=[[750-580],[850-750]],fmt='D',color=colors['pleiades'],capsize=4,ms=5)
cx.text(.02,.97,'c  Cluster mass',transform=cx.transAxes,va='top',fontweight='bold');cx.set_xlabel('Age (Myr)');cx.set_ylabel(r'Mass inside tidal radius ($M_\odot$)');cx.set_xlim(0,690);cx.set_ylim(0,1350)
cx.text(.03,.04,'Diamonds: approximate literature ranges\nDashed: compact start; solid: extended start',transform=cx.transAxes,fontsize=9)

q=ends[(ends.cluster=='hyades') & (ends.binary_fraction==0)].sort_values('initial_rh_pc')
if not len(q):q=ends.sort_values('binary_fraction')
if len(q):
    row=q.iloc[-1];d=pd.read_csv(W/'runs'/row.run/'final_descendants.csv')
    # Sky-independent projection; preserve every tail particle in data ledger.
    stars=(d.mass_Msun>=.1)&(d.mass_Msun<.2);bd=(d.mass_Msun>=.03)&(d.mass_Msun<.08)
    dx.scatter(d.loc[stars,'y_pc'],d.loc[stars,'x_pc'],s=7,color='#8f8f8f',alpha=.35,lw=0,label='0.10–0.20 solar masses')
    dx.scatter(d.loc[bd,'y_pc'],d.loc[bd,'x_pc'],s=10,color=colors[row.cluster],alpha=.6,lw=0,label='0.03–0.08 solar masses')
    theta=np.linspace(0,2*np.pi,300);dx.plot(18*np.cos(theta),18*np.sin(theta),color='black',lw=.9,ls=':')
    dx.set_xlim(-1100,1100);dx.set_ylim(-250,250)
    dx.text(.02,.97,f'd  {row.cluster.capitalize()} descendants at {row.time_Myr:g} Myr',transform=dx.transAxes,va='top',fontweight='bold',fontsize=10)
    dx.legend(loc='upper right',frameon=False,fontsize=8)
dx.set_xlabel('Galactic y relative to cluster (pc)');dx.set_ylabel('Galactic x relative to cluster (pc)')
fig.text(.07,.013,'Exploratory models: one birth mass distribution; selection brackets, not a fitted completeness correction.',fontsize=10,color='#333333')
fig.subplots_adjust(left=.08,right=.98,top=.97,bottom=.085)
for ext in ['png','pdf']:fig.savefig(OUT/f'CLUSTER_DISSOLUTION_20261009.{ext}',dpi=190)
plt.close(fig)
print(OUT/'CLUSTER_DISSOLUTION_20261009.png')
