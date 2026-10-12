"""Render the saved stellar states; no interpolation or extrapolated tails."""
from pathlib import Path
import argparse,csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--omit-masses',nargs='*',type=float,default=[])
parser.add_argument('--output-stem',type=Path,default=ROOT/'boundary')
args=parser.parse_args()
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

tracks=read('tracks.csv','mass_Msun')
palette=dict(zip(sorted(tracks),plt.cm.viridis(np.linspace(.08,.85,max(2,len(tracks))))))
tracks={m:h for m,h in tracks.items()
        if not any(abs(m-omit)<1e-9 for omit in args.omit_masses)}
if not tracks:parser.error('No tracks remain after filtering.')
metadata={e['mass_Msun']:e for e in json.loads((ROOT/'provenance.json').read_text())['endpoints']}
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, LinearSegmentedColormap
fig,(ax,lum)=plt.subplots(2,1,figsize=(11.5,10.4),sharex=True,
                         gridspec_kw={'height_ratios':[1,1.05]})
fig.subplots_adjust(left=.10,right=.78,bottom=.09,top=.965,hspace=.065)
rows=[];points=[]
tnorm=Normalize(650,2700)
tcmap=LinearSegmentedColormap.from_list('stellar_temperature',plt.cm.plasma(np.linspace(.04,.90,256)))
endpoints=[]
for m,h in sorted(tracks.items()):
 color=palette[m]
 x=np.array([a['age_Gyr'] for a in h]);y=np.array([a['nuclear_fraction']for a in h]);mask=x>=.05
 ax.plot(x[mask],y[mask],color=color,lw=2 if m==.05825 else 1.6,label=f'{m:.4g}')
 last=h[-1];full=last['age_Gyr']>=1000*(1-1e-10) or metadata[m].get('stop','').startswith('nuclear fraction below');marker='s' if metadata[m].get('running') else 'o'
 ax.plot(x[-1],y[-1],marker=marker,ms=5,mec=color,mfc=color if full else 'white',zorder=5)
 ll=np.log10([a['luminosity_Lsun'] for a in h]);tt=np.array([a['Teff_K']for a in h])
 xy=np.column_stack((x[mask],ll[mask]));segments=np.stack((xy[:-1],xy[1:]),axis=1)
 lc=LineCollection(segments,cmap=tcmap,norm=tnorm,linewidths=1.8 if m!=.05825 else 2.5,capstyle='round')
 lc.set_array(.5*(tt[mask][:-1]+tt[mask][1:]));lum.add_collection(lc)
 color_end=tcmap(tnorm(tt[-1]))
 lum.plot(x[-1],ll[-1],marker=marker,ms=5,mec=color_end,mfc=color_end if full else 'white',zorder=5)
 endpoints.append((m,x[-1],ll[-1],full))
ax.set(xscale='log',xlim=(.05,1150),ylim=(-.025,max(1.12 if args.omit_masses else 1.07,1.035*max(a['nuclear_fraction'] for h in tracks.values() for a in h if a['age_Gyr']>=.05))),ylabel=r'$L_{\rm nuc}/L$')
ax.axhline(1,color='.35',ls=':',lw=.8)
ax.text(.025,.955,r'$3Z_{\rm ref}$ · cloudy, $f_{\rm sed}=2$',transform=ax.transAxes,va='top',fontsize=14)
ax.legend(title=r'Mass ($M_\odot$)',loc='upper left',bbox_to_anchor=(1.015,1),frameon=False,labelspacing=.38)
ax.text(1.03,-.10,'Filled: 1 Tyr or\nfaded below 0.1%\nOpen: incomplete',transform=ax.transAxes,fontsize=10.5,linespacing=1.5,va='top')
if .05815 in tracks:
 last=tracks[.05815][-1]
 ax.annotate('0.05815',xy=(last['age_Gyr'],last['nuclear_fraction']),
             xytext=(13,24),textcoords='offset points',fontsize=10,va='center',
             arrowprops=dict(arrowstyle='-',lw=.65,color='.35',shrinkA=2,shrinkB=4))
lum.set(ylim=(-6.2,-2.55),xlabel='Age (Gyr)',ylabel=r'$\log_{10}(L/L_\odot)$')
lum.set_xticks([.1,1,10,100,1000],labels=['0.1','1','10','100','1000'])
for a in [ax,lum]:a.grid(alpha=.13)
# Identify the warming borderline track and the mass range without a second large legend.
for mass,offset in [(.054,(-48,0)),(.058,(-36,17)),(.0582,(-70,25)),(.05825,(4,20)),(.0583,(28,-2)),(.05835,(-48,-10)),(.0585,(-5,-25)),(max(tracks),(-42,20))]:
 end=next((e for e in endpoints if e[0]==mass),None)
 if end is None:continue
 m,x,y,_=end
 lum.annotate(f'{m:.4g}',xy=(x,y),xytext=offset,textcoords='offset points',fontsize=10,
              ha='center',va='center',color='black',bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=.8),
              arrowprops=dict(arrowstyle='-',lw=.65,color='.35',shrinkA=2,shrinkB=4))
lum.text(.50,.95,r'Labels: mass in $M_\odot$',transform=lum.transAxes,va='top',fontsize=10.5)
cax=fig.add_axes([.815,.18,.017,.245])
cb=fig.colorbar(plt.cm.ScalarMappable(norm=tnorm,cmap=tcmap),cax=cax)
cb.set_ticks([700,1000,1500,2000,2500]);cb.set_label(r'$T_{\rm eff}$ (K)',labelpad=8)
note=('Instantaneous convection; no semiconvection. Mixing and cloud clearing remain uncertain.'
      if args.omit_masses else
      'Instantaneous convection; no semiconvection. Sharp merger features remain unresolved.')
fig.text(.10,.025,note,fontsize=10.5,color='.25')
args.output_stem.parent.mkdir(parents=True,exist_ok=True)
fig.savefig(args.output_stem.with_suffix('.png'))
fig.savefig(args.output_stem.with_suffix('.pdf'))
plt.close(fig)
if args.omit_masses:raise SystemExit(0)

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
