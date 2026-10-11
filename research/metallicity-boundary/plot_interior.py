"""Render the compact accepted-model profiles for the borderline object."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
DEST=Path(__file__).resolve().parent
data=np.load(DEST/'interior.npz')
ages=data['age_Gyr'];q=data['q'];conv=data['convective']
raw={'T_K':data['T_K'],'X':data['X']}
kind=np.where(conv,0,np.where(data['conductive_share']>=.5,1,2))
mass_edges=np.r_[0.,.5*(q[:-1]+q[1:]),1.]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,
                     'axes.labelsize':14,'axes.spines.top':False,
                     'axes.spines.right':False,'savefig.dpi':200})
fig, axes = plt.subplots(3, 1, figsize=(11.2, 10.6), sharex=True,
                          gridspec_kw={'height_ratios':[.85,1.4,1.3]})
fig.subplots_adjust(left=.105, right=.855, top=.90, bottom=.085, hspace=.09)
a, b, c = axes
f = data['nuclear_fraction']
teff = data['Teff_K']
a.plot(ages, f, color='#292929', lw=2, label=r'$L_{\rm nuc}/L$')
a.axhline(1, color='.45', lw=.8, ls=':')
a.set(ylabel=r'$L_{\rm nuc}/L$', ylim=(0,1.10))
a.set_yticks([0,.5,1])
a.grid(axis='y', alpha=.15)
# The earliest Hayashi temperature is not relevant to the late warming here.
twin = a.twinx()
twin.spines['right'].set_visible(True)
warm = ages >= .5
twin.plot(ages[warm], teff[warm], color='#a35b26', lw=1.4, ls='--')
twin.set(ylabel=r'$T_{\rm eff}$ (K)', ylim=(950,2300))
twin.set_yticks([1000,1500,2000])
a.text(.52,.48, r'Black: $L_{\rm nuc}/L$   ·   dashed: $T_{\rm eff}$',
       transform=a.transAxes, fontsize=10.5,
       bbox={'facecolor':'white','edgecolor':'none','alpha':.8,'pad':1})

palette = ListedColormap(['#c0d5e7','#edc99d','#ffffff'])
norm = BoundaryNorm([-.5,.5,1.5,2.5],3)
time_edges = np.r_[ages[0], .5*(ages[:-1]+ages[1:]), ages[-1]]
b.pcolormesh(time_edges, mass_edges, kind.T, cmap=palette, norm=norm,
             shading='flat', rasterized=True)
if np.any(~conv):
    b.contour(ages, q, conv.T.astype(float), levels=[.5], colors='#344859', linewidths=.8)
temperature = b.contour(ages, q, (raw['T_K']/1e6).T, levels=[.5,1,2,3],
                        colors='.3', linewidths=.7, linestyles='--')
positions=[]
for level in [.5,1,2]:
    eligible=np.flatnonzero((np.min(raw['T_K'],axis=1)/1e6<level)&(np.max(raw['T_K'],axis=1)/1e6>level)&(ages>.7*ages[-1]))
    if len(eligible):
        row=eligible[np.argmin(np.abs(ages[eligible]-.85*ages[-1]))]
        temp=raw['T_K'][row]/1e6
        positions.append((ages[row],np.interp(level,temp[::-1],q[::-1])))
b.clabel(temperature,fmt=lambda v:f'{v:g} MK',fontsize=10,inline=True,manual=positions)
if np.max(raw['T_K'][-1])/1e6>3:
    row=len(ages)-1
    pos=np.interp(3,(raw['T_K'][row]/1e6)[::-1],q[::-1])
    b.annotate('3 MK',xy=(ages[row],pos),xytext=(-42,16),textcoords='offset points',
               ha='center',fontsize=10,color='.3',arrowprops=dict(arrowstyle='-',lw=.65,color='.3'))
b.set(ylabel=r'Enclosed mass $m/M$', ylim=(0,1))
for v in [.25,.45]:b.axhline(v,color='.5',lw=.55,ls=':')

im = c.pcolormesh(time_edges, mass_edges, raw['X'].T, cmap='viridis',
                  vmin=.56, vmax=.75, shading='flat', rasterized=True)
if np.any(~conv):
    c.contour(ages, q, conv.T.astype(float), levels=[.5], colors='white', linewidths=1.05)
c.set(ylabel=r'Enclosed mass $m/M$', ylim=(.25,.45), xlabel='Age (Gyr)')
c.set_yticks([.25,.30,.35,.40,.45])
c.text(.025,.95,'Hydrogen across the separating layer',transform=c.transAxes,
       va='top',fontsize=12,color='black',
       bbox={'facecolor':'white','edgecolor':'none','alpha':.82,'pad':3})
cax=fig.add_axes([.89,.095,.017,.265])
bar=fig.colorbar(im,cax=cax);bar.set_label('Hydrogen mass fraction X',fontsize=12)
bar.set_ticks([.56,.60,.65,.70,.75])

last=ages[-1]
a.set_xlim(0,last)
if last>100:c.set_xticks(np.arange(0,last+1,20))
cross=np.flatnonzero((f>=1)&(ages>1))
if len(cross):
    j=cross[0]
    for panel in axes:panel.axvline(ages[j],color='.25',lw=.7,ls=':',alpha=.65)
    a.annotate('Nuclear balance',xy=(ages[j],f[j]),xytext=(ages[j]-15,.70),fontsize=10.5,
               arrowprops=dict(arrowstyle='-|>',lw=.8,mutation_scale=9,color='.2'))
if last>80:
    j=np.argmax(np.where(ages>1,f,-1))
    a.annotate(f'{f[j]:.4g}',xy=(ages[j],f[j]),xytext=(ages[j]-12,.83),fontsize=10.5,
               arrowprops=dict(arrowstyle='-|>',lw=.8,mutation_scale=9,color='.2'))
fig.text(.105,.966,r'$0.05825\,M_\odot$ · $3Z_{\rm ref}$ · cloudy, $f_{\rm sed}=2$',fontsize=16)
fig.legend(handles=[Patch(facecolor=palette(0),label='Convective'),
                    Patch(facecolor=palette(1),label='Nonconvective: conduction'),
                    Patch(facecolor='white',edgecolor='.5',label='Nonconvective: radiation')],
           loc='upper left',bbox_to_anchor=(.098,.951),ncol=3,frameon=False,fontsize=10.5)
fig.text(.105,.025,'Lower panel enlarges the dotted mass range; white edges mark the convection boundary.',fontsize=10.5)
for ext in ['png','pdf']:fig.savefig(DEST/f'Ember_boundary_interior.{ext}')
plt.close(fig)
