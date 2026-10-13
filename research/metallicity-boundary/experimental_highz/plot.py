"""Plot the retained flux-checked five-times-reference-metallicity controls."""
from pathlib import Path
import csv, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

W=Path(__file__).resolve().parent
rows=list(csv.DictReader((W/'tracks.csv').open()))
meta=json.loads((W/'provenance.json').read_text())
plt.rcParams.update({'font.size':12.5,'axes.labelsize':15,'legend.fontsize':11.5,
 'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.facecolor':'white'})
fig,(ax,lum)=plt.subplots(2,1,figsize=(11.4,8.8),sharex=True)
fig.subplots_adjust(left=.10,right=.79,bottom=.145,top=.975,hspace=.065)
cmap=plt.get_cmap('plasma');norm=Normalize(900,2700)
color={k:plt.get_cmap('viridis')(v) for k,v in
       [('006',.15),('0061',.28),('006125',.41),('00614',.48),('00615',.55),('0062',.69),('0064',.80)]}
handles=[]
for r in meta['source_cases']:
 case=r['case'];main=r['tau']==3
 rr=[v for v in rows if v['case']==case and float(v['age_Gyr'])>=.1]
 x=np.array([float(v['age_Gyr'])for v in rr]);y=np.array([float(v['nuclear_fraction'])for v in rr])
 ll=np.log10([float(v['luminosity_Lsun'])for v in rr]);tt=np.array([float(v['Teff_K'])for v in rr])
 if not main:
  for a,b in [(ax,y),(lum,ll)]:a.plot(x,b,color='.4',lw=1.5,ls='--',alpha=.55,zorder=1)
  handles.append(Line2D([],[],c='.4',lw=1.5,ls='--',alpha=.65,label=r'0.06; $\tau_R=10$'))
  continue
 ax.plot(x,y,color=color[case],lw=2.2,zorder=3)
 xy=np.column_stack([x,ll]);lc=LineCollection(np.stack([xy[:-1],xy[1:]],axis=1),cmap=cmap,
                 norm=norm,linewidths=2.4,capstyle='round')
 lc.set_array(.5*(tt[:-1]+tt[1:]));lum.add_collection(lc)
 stopped=not bool(r['report']['requested_age_reached'])
 for a,yy,c in [(ax,y,color[case]),(lum,ll,cmap(norm(tt[-1])))]:
  a.plot(x[-1],yy[-1],marker='o',ms=5.5,mfc='white'if stopped else c,mec=c,mew=1.2,zorder=5)
 handles.append(Line2D([],[],c=color[case],lw=2.2,label=f"{r['mass_Msun']:.4g}"))
ax.axhline(1,color='.35',lw=.8,ls=':',zorder=0)
ax.set(xscale='log',xlim=(.1,1300),ylim=(-.025,1.13),ylabel=r'$L_{\rm nuc}/L$')
ax.text(.025,.955,r'$5Z_{\rm ref}$ · clear',transform=ax.transAxes,fontsize=14,va='top')
ax.text(.025,.84,r'$Z=0.08241,\quad Y=0.27$',transform=ax.transAxes,fontsize=12)
ax.legend(handles=handles,title=r'Mass ($M_\odot$)',loc='upper left',
          bbox_to_anchor=(1.015,1),frameon=False,handlelength=2.1)
ax.text(1.045,.31,'Experimental tracks\n\nOpen circles:\ntable coverage limits',
        transform=ax.transAxes,fontsize=11.5,linespacing=1.5,va='top')
lum.set(ylim=(-5.60,-2.75),xlabel='Age (Gyr)',ylabel=r'$\log_{10}(L/L_\odot)$')
lum.set_xticks([.1,1,10,100,1000],labels=['0.1','1','10','100','1000'])
for axis in(ax,lum):axis.grid(alpha=.13,lw=.6)
cax=fig.add_axes([.835,.18,.018,.265]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax)
cb.set_label(r'$T_{\rm eff}$ (K)');cb.set_ticks([900,1300,1700,2100,2500])
fig.text(.10,.045,'Post-deuterium starts; solid curves match the atmosphere at Rosseland optical depth 3.',fontsize=11)
fig.text(.10,.021,'Local H/He atmosphere response is provisional.',fontsize=11)
for ext in('png','pdf'):fig.savefig(W/f'Ember_boundary_5Z_experimental.{ext}',dpi=220)
plt.close(fig)
