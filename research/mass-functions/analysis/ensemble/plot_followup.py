from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter

W=Path(__file__).resolve().parent;D=W.parents[1]/'docs/research/independent'
e=pd.read_csv(W/'structural_endpoints.csv')
e=e[~((e.binary_mode=='random')&(e.binary_fraction>0))]
j=pd.read_csv(W/'joint_checks.csv');obs=pd.read_csv(W/'observed_matched_counts.csv')
mp=pd.read_csv(W/'model_angular_profile_draws.csv');op=pd.read_csv(W/'observed_angular_profile_draws.csv')
plan=json.loads((W/'analysis_plan.json').read_text())
plt.rcParams.update({'font.size':11,'axes.labelsize':12,'axes.spines.top':False,
    'axes.spines.right':False,'savefig.dpi':190,'legend.fontsize':9.5})
fig,axs=plt.subplots(3,2,figsize=(12,12),gridspec_kw={'height_ratios':[1,1,1]})
colors={'hyades':'#b55225','pleiades':'#246d98'}
for col,cl in enumerate(['hyades','pleiades']):
    color=colors[cl];g=e[e.cluster==cl];crit=plan['structural_screen'][cl]
    ax=axs[0,col];mb=crit['mass_excluding_BH_NS_Msun'];rb=crit['half_mass_radius_pc']
    ax.add_patch(Rectangle((rb[0],mb[0]),rb[1]-rb[0],mb[1]-mb[0],color=color,alpha=.10,zorder=0))
    for r in g.itertuples():
        new=r.workspace.endswith('v2');passed=r.structural_screen_pass;low=r.expected_birth_ratio<.6
        ax.scatter(r.aperture_halfmass_pc,r.aperture_mass_excluding_BH_NS_Msun,
            marker='^' if low else ('s' if new else 'o'),s=48,facecolor=('#785b8f' if low else color) if passed else 'white',
            edgecolor=color if new else '#777777',linewidth=1.2,zorder=3)
        if new:
            ax.annotate(str(r.seed),(r.aperture_halfmass_pc,r.aperture_mass_excluding_BH_NS_Msun),
                xytext=(5,5),textcoords='offset points',fontsize=9)
    ax.set_xlim(2.3,6.25)
    ax.set_ylim((50,450) if cl=='hyades' else (500,1020))
    ax.set_xlabel('Half-mass radius in the aperture (pc)')
    ax.set_ylabel(r'Mass in the aperture ($M_\odot$)')
    ax.text(.02,.97,f'{"ab"[col]}  {cl.capitalize()}',transform=ax.transAxes,va='top',weight='bold',fontsize=13)
    ax.text(.98,.05,f'Aperture: {crit["aperture_pc"]:g} pc',transform=ax.transAxes,ha='right',fontsize=10)
    if col==0:ax.legend(handles=[Line2D([],[],marker='s',color='none',markeredgecolor='k',markerfacecolor='white',label='New realizations (number = seed)'),Line2D([],[],marker='o',color='none',markeredgecolor='grey',label='Initial grid')],loc='upper right',frameon=False)

    ax=axs[1,col]
    o=op[(op.cluster==cl)&(op.scenario=='r200_union')&(op.population=='stars')]
    mat=o.pivot(index='rep',columns='bin',values='N').to_numpy()
    cdf=np.cumsum(mat,axis=1)/mat.sum(axis=1)[:,None]
    x=np.array([5,10,20,40,80,180] if cl=='hyades' else [1,2,3,5,10,30,90,180]);q=np.quantile(cdf,[.16,.5,.84],axis=0)
    ax.fill_between(x,q[0],q[2],color=color,alpha=.18)
    ax.plot(x,q[1],color=color,lw=2.5,marker='o',ms=4,label='Observed reference stars',zorder=5)
    for r in g.itertuples():
        p=mp[(mp.run==r.run)&(mp.cloud=='nc')&(mp.scenario=='wide_union_shallow')&(mp.population=='stars')&(mp.faint_extension!='continued_slope')]
        mm=p.groupby('bin').N.mean().to_numpy()
        ax.plot(x,np.cumsum(mm)/mm.sum(),color=('#785b8f' if r.expected_birth_ratio<.6 else 'black') if r.structural_screen_pass else '#777777',
            alpha=.65 if r.structural_screen_pass else .22,lw=1.3 if r.structural_screen_pass else .8,
            ls='-' if r.workspace.endswith('v2') else '--')
    ax.set_xscale('log');ax.set_xticks(x);ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_xlim((4,190) if cl=='hyades' else (.8,190));ax.set_ylim(0,1.04)
    ax.set_xlabel('Angular separation from cluster center (degrees)');ax.set_ylabel('Cumulative fraction of reference stars')
    ax.text(.02,.97,f'{"cd"[col]}  Stellar spatial check',transform=ax.transAxes,va='top',weight='bold')
    if col==0:ax.legend(loc='lower right',frameon=False)
    ax.text(.97,.18,r'$0.10$–$0.20\,M_\odot$',transform=ax.transAxes,ha='right',fontsize=11)

    ax=axs[2,col];ax.axhspan(.55,.71,color='#a0a0a0',alpha=.17);ax.axhline(.63,color='#777777',ls=':',lw=1)
    scenarios=['central_matched','ir_matched','wide'];labels=['10 pc','IR region','Wide union'] if cl=='hyades' else ['Central 3°','DANCe, $K_s<17$','Wide union']
    for idx,scenario in enumerate(scenarios):
        os='r200_union' if scenario=='wide' else scenario
        oo=obs[(obs.cluster==cl)&(obs.scenario==os)].iloc[0]
        for k,r in enumerate(g.itertuples()):
            rr=j[(j.run==r.run)&(j.faint_extension!='continued_slope')]
            rr=rr[rr.scenario.str.startswith('wide')] if scenario=='wide' else rr[rr.scenario==scenario]
            xx=idx+(k-(len(g)-1)/2)*.035
            lo,hi=rr.model_ratio.min(),rr.model_ratio.max()
            ax.plot([xx,xx],[lo,hi],color=('#785b8f' if r.expected_birth_ratio<.6 else color) if r.structural_screen_pass else '#999999',
                alpha=.85 if r.structural_screen_pass else .35,lw=2 if r.structural_screen_pass else 1.1,solid_capstyle='round')
        ax.errorbar(idx+.2,oo.ratio,yerr=[[oo.ratio-oo.ratio_lo],[oo.ratio_hi-oo.ratio]],fmt='D',ms=5,
            color='black',capsize=3,lw=1.3,zorder=6)
    ax.set_xticks(range(3),labels);ax.set_ylim(0,1.06);ax.set_xlim(-.3,2.4)
    ax.set_ylabel('Brown-dwarf / reference-star ratio')
    ax.text(.02,.97,f'{"ef"[col]}  Brown-dwarf prediction',transform=ax.transAxes,va='top',weight='bold')
    if col==0:ax.legend(handles=[Line2D([],[],color=color,lw=3,label='Inside structural range'),Line2D([],[],color='#999999',lw=2,label='Outside range'),Line2D([],[],marker='D',color='k',lw=0,label='Observed candidates')],loc='upper right',frameon=False)
    if col==1:ax.legend(handles=[Line2D([],[],color='#785b8f',lw=3,label='Lower field normalization (0.55)')],loc='upper right',bbox_to_anchor=(1,.90),frameon=False)
fig.subplots_adjust(left=.08,right=.98,bottom=.075,top=.975,hspace=.32,wspace=.24)
fig.text(.08,.032,'All escapers retained. Grey band: field component ratio. Cluster predictions include unresolved binaries.',fontsize=10,color='#444444')
fig.text(.08,.014,'Structural boxes are broad checks; model spans show survey and atmosphere choices, not full uncertainties.',fontsize=10,color='#444444')
stem=D/'CLUSTER_DISSOLUTION_FOLLOWUP_20261009'
fig.savefig(stem.with_suffix('.png'));fig.savefig(stem.with_suffix('.pdf'));print(stem.with_suffix('.png'))
