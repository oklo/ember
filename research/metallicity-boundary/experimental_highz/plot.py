"""Plot the retained five-times-reference-metallicity controls."""
from pathlib import Path
import csv, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

W = Path(__file__).resolve().parent
rows = list(csv.DictReader((W/'tracks.csv').open()))
meta = json.loads((W/'provenance.json').read_text())
plt.rcParams.update({'font.size':12.5, 'axes.labelsize':15,
                     'legend.fontsize':11.5, 'axes.spines.top':False,
                     'axes.spines.right':False, 'pdf.fonttype':42,
                     'savefig.facecolor':'white'})
fig, (ax, lum) = plt.subplots(2, 1, figsize=(11.4,8.8), sharex=True)
fig.subplots_adjust(left=.10, right=.79, bottom=.115, top=.975, hspace=.065)
cmap = plt.get_cmap('plasma'); norm = Normalize(1000,2700)
colors = {'tau1':'#27867a', 'tau3':'#355f8b', 'tau3_original':'.55'}
for case in ('tau3_original','tau1','tau3'):
    rr = [r for r in rows if r['case']==case and float(r['age_Gyr'])>=.1]
    x = np.array([float(r['age_Gyr']) for r in rr])
    y = np.array([float(r['nuclear_fraction']) for r in rr])
    ll = np.log10([float(r['luminosity_Lsun']) for r in rr])
    tt = np.array([float(r['Teff_K']) for r in rr])
    if case=='tau3_original':
        ax.plot(x,y,color=colors[case],lw=3.5,alpha=.38,ls='--',zorder=1)
        lum.plot(x,ll,color=colors[case],lw=3.5,alpha=.3,ls='--',zorder=1)
        continue
    ax.plot(x,y,color=colors[case],lw=2.3,zorder=3)
    xy = np.column_stack([x,ll])
    lc = LineCollection(np.stack([xy[:-1],xy[1:]],axis=1),cmap=cmap,
                        norm=norm,linewidths=2.5 if case=='tau3' else 1.8,
                        linestyles='solid' if case=='tau3' else 'dashed',
                        capstyle='round',alpha=1 if case=='tau3' else .8)
    lc.set_array(.5*(tt[:-1]+tt[1:]));lum.add_collection(lc)
    for axis, yy, color in ((ax,y,colors[case]),(lum,ll,cmap(norm(tt[-1])))):
        stopped = meta['endpoints'][case]['stop_reason'] != 'requested age'
        axis.plot(x[-1],yy[-1],marker='o',ms=6,
                  mfc='white' if stopped else color,mec=color,mew=1.4,zorder=5)
    if case=='tau1':
        ax.annotate(f'Composition limit · {x[-1]:.4g} Gyr',
                    (x[-1],y[-1]),xytext=(35,.77),textcoords='data',fontsize=11.5,
                    arrowprops=dict(arrowstyle='-|>',lw=.85,color='.2',mutation_scale=10),
                    ha='left',va='center')
    else:
        lum.annotate(f'{x[-1]:.4g} Gyr · {tt[-1]:.4g} K',
                     (x[-1],ll[-1]),xytext=(-20,28),textcoords='offset points',
                     fontsize=11.5,ha='right',va='center',
                     arrowprops=dict(arrowstyle='-|>',lw=.85,color='.2',mutation_scale=10))
ax.axhline(1,color='.35',lw=.8,ls=':',zorder=0)
ax.set(xscale='log',xlim=(.1,1300),ylim=(-.025,1.13),ylabel=r'$L_{\rm nuc}/L$')
ax.text(.025,.955,r'$5Z_{\rm ref}$ · clear · $0.06\,M_\odot$',
        transform=ax.transAxes,fontsize=14,va='top')
ax.text(.025,.82,r'$Z=0.08241,\quad Y=0.27$',transform=ax.transAxes,fontsize=12)
handles=[Line2D([],[],c=colors['tau3'],lw=2.3,label=r'$\tau_R=3$'),
         Line2D([],[],c=colors['tau1'],lw=2.3,label=r'$\tau_R=1$'),
         Line2D([],[],c='.55',lw=2,alpha=.6,ls='--',label='Original EOS')]
ax.legend(handles=handles,title='Atmosphere match',loc='upper left',
          bbox_to_anchor=(1.02,1),frameon=False,handlelength=2.2)
ax.text(1.045,.47,'Experimental tracks\n\nOpen circles:\ntable limits',transform=ax.transAxes,fontsize=11.5,
        linespacing=1.45)
lum.set(ylim=(-5.20,-2.75),xlabel='Age (Gyr)',ylabel=r'$\log_{10}(L/L_\odot)$')
lum.set_xticks([.1,1,10,100,1000],labels=['0.1','1','10','100','1000'])
for axis in (ax,lum):axis.grid(alpha=.13,lw=.6)
cax = fig.add_axes([.835,.18,.018,.265])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax)
cb.set_label(r'$T_{\rm eff}$ (K)');cb.set_ticks([1000,1400,1800,2200,2600])
fig.text(.10,.035,'Post-deuterium starts; computed states only. Atmosphere and EOS source checks continue.',fontsize=11)
for ext in ('png','pdf'):fig.savefig(W/f'Ember_boundary_5Z_experimental.{ext}',dpi=220)
plt.close(fig)
