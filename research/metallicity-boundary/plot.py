"""Render the saved stellar states; no interpolation or extrapolated tails."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
plt.rcParams.update({'font.size':12,'axes.labelsize':14,'legend.fontsize':10.5,
                     'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':200})

def read(name,key):
    groups={}
    with (ROOT/name).open() as f:
        for row in csv.DictReader(f):
            k=float(row[key]) if key=='mass_Msun' else row[key]
            groups.setdefault(k,[]).append({x:float(v) for x,v in row.items() if x!=key})
    return {k:sorted(v,key=lambda r:r['age_Gyr']) for k,v in groups.items()}

def line(ax,rows,y,**style):
    ax.plot([r['age_Gyr'] for r in rows],[r[y] for r in rows],**style)

tracks=read('tracks.csv','mass_Msun');early=read('early.csv','mass_Msun')
fig=plt.figure(figsize=(11.5,10));ax=fig.add_axes([.085,.43,.695,.52]);bottom=fig.add_axes([.085,.095,.695,.24])
colors=dict(zip(sorted(tracks),plt.cm.viridis(np.linspace(.08,.85,len(tracks)))))
for mass,h in sorted(tracks.items()):
    color=colors[mass];line(ax,h,'nuclear_fraction',c=color,lw=2 if mass==.058 else 1.55,label=f'{mass:.4g}')
    end=h[-1];complete=end['age_Gyr']>=1000*(1-1e-10)
    ax.plot(end['age_Gyr'],end['nuclear_fraction'],'o',ms=5.5,mec=color,mfc=color if complete else 'white')
ax.set(xscale='log',xlim=(.05,1150),ylim=(-.025,1.07),xlabel='Age (Gyr)',ylabel=r'$L_{\rm nuc}/L$')
ax.set_xticks([.1,1,10,100,1000],labels=['0.1','1','10','100','1000'])
ax.axhline(1,c='.35',ls=':',lw=.8);ax.grid(alpha=.13)
ax.text(.03,.895,r'$3Z_{\rm ref}$ · cloudy, $f_{\rm sed}=2$',transform=ax.transAxes,va='top',fontsize=14)
ax.legend(title=r'Mass ($M_\odot$)',loc='upper left',bbox_to_anchor=(1.015,1),frameon=False)
ax.text(1.03,.25,'Filled: 1 Tyr\nOpen: incomplete',transform=ax.transAxes,fontsize=10.5,linespacing=1.5)
for mass,h in sorted(early.items()):
    line(bottom,tracks[mass],'nuclear_fraction',c=colors[mass],alpha=.45,lw=1.2,ls='--')
    line(bottom,h,'nuclear_fraction',c=colors[mass],lw=2,label=f'{mass:.4g}')
bottom.set(xlim=(.5,5),ylim=(.025,.91),xlabel='Age (Gyr)',ylabel=r'$L_{\rm nuc}/L$');bottom.grid(alpha=.13)
bottom.legend(title=r'Mass ($M_\odot$)',loc='upper left',bbox_to_anchor=(1.015,1),frameon=False)
bottom.text(1.03,.36,'Solid: 50 Myr steps\nDashed: coarser steps',transform=bottom.transAxes,fontsize=10.5,linespacing=1.5)
fig.text(.085,.365,'Early evolution — independent timestep checks',fontsize=12)
fig.text(.085,.025,'Open endpoints are incomplete. The early checks are not joined to the older tracks.',fontsize=10,color='.25')
fig.savefig(ROOT/'boundary.png');plt.close(fig)

data=read('timestep_comparison.csv','case')
fig,axes=plt.subplots(1,2,figsize=(12,5.6));fig.subplots_adjust(left=.075,right=.975,bottom=.24,top=.92,wspace=.27)
for case,color,style in [('2 Gyr','#b35c1e','--'),('0.5 Gyr','#b35c1e','-'),('0.25 Gyr','#5b327e','-'),('No microscopic transport','#287998','-')]:
    label='Microscopic transport disabled' if case.startswith('No') else f'Settling; {case} maximum step'
    for a,key in zip(axes,['nuclear_fraction','central_T_MK']):
        line(a,data[case],key,c=color,ls=style,lw=2,label=label)
for a in axes:a.set(xlim=(5,100),xlabel='Age (Gyr)');a.grid(alpha=.13)
axes[0].set(ylabel=r'$L_{\rm nuc}/L$',ylim=(-.03,1.5));axes[1].set(ylabel=r'Central temperature ($10^6$ K)')
axes[0].text(.035,.96,r'$0.058\,M_\odot$ · $3Z_{\rm ref}$ · $f_{\rm sed}=2$',transform=axes[0].transAxes,va='top')
handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=2,frameon=False,fontsize=10.5)
fig.savefig(ROOT/'timestep_comparison.png');plt.close(fig)
