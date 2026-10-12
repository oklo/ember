"""Plot retained current-physics cloudy models from the accompanying CSV."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

W=Path(__file__).resolve().parent
rows=list(csv.DictReader((W/'tracks.csv').open()))
plt.rcParams.update({'font.size':12,'axes.labelsize':14,'legend.fontsize':11,
                     'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':200})
fig,(ax,lum)=plt.subplots(2,1,figsize=(11.4,9.1),sharex=True)
fig.subplots_adjust(left=.10,right=.81,bottom=.12,top=.975,hspace=.07)
colors=['#5364ad','#319889','#b17926']
cmap=plt.get_cmap('plasma'); norm=Normalize(700,2800)
ends=[]
for mass,color in zip([.060,.0605,.061],colors):
    rr=[r for r in rows if float(r['mass_Msun'])==mass and float(r['age_Gyr'])>=.1]
    x=np.array([float(r['age_Gyr']) for r in rr]); y=np.array([float(r['nuclear_fraction']) for r in rr])
    ll=np.log10([float(r['luminosity_Lsun']) for r in rr]); tt=np.array([float(r['Teff_K']) for r in rr])
    provisional=np.array([int(r['provisional']) for r in rr],dtype=bool)
    for special in [False,True]:
        keep=np.flatnonzero(provisional==special)
        if not len(keep):continue
        if special and keep[0]>0:keep=np.r_[keep[0]-1,keep]
        ax.plot(x[keep],y[keep],color=color,lw=2.2,ls='--' if special else '-')
        xy=np.column_stack([x[keep],ll[keep]])
        if special:
            lum.plot(x[keep],ll[keep],color=cmap(norm(np.mean(tt[keep]))),lw=2.3,ls='--')
        else:
            lc=LineCollection(np.stack([xy[:-1],xy[1:]],axis=1),cmap=cmap,norm=norm,
                              linewidths=2.3,capstyle='round')
            lc.set_array(.5*(tt[keep][:-1]+tt[keep][1:]));lum.add_collection(lc)
    marker='s' if mass==.061 else 'o'
    for axis,yy,cc in [(ax,y,color),(lum,ll,cmap(norm(tt[-1])))]:
        axis.plot(x[-1],yy[-1],marker=marker,ms=6,mec=cc,
                  mfc='white' if mass==.0605 else cc,zorder=5)
    label_age={.060:x[-1],.0605:463.,.061:180.}[mass]
    j=int(np.argmin(abs(x-label_age)))
    ends.append((mass,x[j],ll[j]))
ax.axhline(1,color='.4',lw=.8,ls=':')
ax.set(xscale='log',xlim=(.1,1600),ylim=(-.025,1.055),ylabel=r'$L_{\rm nuc}/L$')
ax.text(.025,.88,r'$3Z_{\rm ref}$ · cloudy, $f_{\rm sed}=2$',transform=ax.transAxes,fontsize=14)
handles=[Line2D([],[],c=c,lw=2.2,label=f'{m:.4g}'+(' *' if m==.0605 else ''))
         for m,c in zip([.060,.0605,.061],colors)]
ax.legend(handles=handles,title=r'Mass ($M_\odot$)',loc='upper left',bbox_to_anchor=(1.015,1),frameon=False)
ax.text(1.035,.57,'* Extended EOS',transform=ax.transAxes,fontsize=10.5)
ax.text(1.035,.36,'Dashed:\natmosphere check\npending',transform=ax.transAxes,fontsize=10.5,linespacing=1.45)
lum.set(ylim=(-6.05,-2.8),xlabel='Age (Gyr)',ylabel=r'$\log_{10}(L/L_\odot)$')
lum.set_xticks([.1,1,10,100,1000],labels=['0.1','1','10','100','1000'])
for m,x,y in ends:
    offset={.060:(-20,18),.0605:(-5,23),.061:(0,-24)}[m]
    lum.annotate(f'{m:.4g}',(x,y),xytext=offset,textcoords='offset points',fontsize=11,
                 ha='center',va='center',bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=1),
                 arrowprops=dict(arrowstyle='-',color='.4',lw=.7))
for axis in [ax,lum]:axis.grid(alpha=.14)
cax=fig.add_axes([.845,.17,.018,.25]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax)
cb.set_label(r'$T_{\rm eff}$ (K)');cb.set_ticks([800,1200,1600,2000,2400,2800])
fig.text(.10,.035,'Metal ionization included. Curves stop at computed states.',fontsize=11)
fig.savefig(W/'Ember_boundary_3Z_cloudy_current.png')
fig.savefig(W/'Ember_boundary_3Z_cloudy_current.pdf')
