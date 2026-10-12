"""Independent meshes; lines connect accepted states without smoothing."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

W=Path(__file__).resolve().parent
plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,2,figsize=(12,7),sharex='col',
                     gridspec_kw={'width_ratios':[1.15,1]})
for prefix,label,color in [('mesh','Local convection','#b86b2c'),
                           ('mesh_coupled','Heat + composition mixing','#287f8e')]:
    for n,style,width,alpha in [(536,'--',1.4,.65),(1071,'-',1.8,1.)]:
        a=np.genfromtxt(W/f'{prefix}_{n}.csv',delimiter=',',names=True)
        for row,key in enumerate(['nuclear_fraction','luminosity_Lsun']):
            y=a[key] if key=='nuclear_fraction' else np.log10(a[key])
            for col in ([0,1] if prefix=='mesh_coupled' else [0]):
                ax=axes[row,col]
                ax.plot(a['age_Gyr'],y,color=color,lw=width,ls=style,alpha=alpha,
                        label=f'{label}, {n} points')
                ax.plot(a['age_Gyr'][-1],y[-1],'o',color=color,alpha=alpha,ms=4)
for ax in axes[0]:ax.axhline(1,color='.5',ls=':',lw=1)
axes[0,0].set_ylabel(r'$L_{\rm nuc}/L$')
axes[0,0].legend(frameon=False,loc='lower left',fontsize=9)
axes[1,0].set(ylabel=r'$\log_{10}(L/L_\odot)$',xlabel='Age (Gyr)',xlim=(20,122))
axes[1,1].set(xlabel='Age (Gyr)',xlim=(120,1010))
axes[0,0].set_ylim(.45,1.1)
axes[1,0].set_ylim(-5.05,-4.30)
axes[0,1].set_ylim(.965,1.012)
axes[1,1].set_ylim(-4.95,-4.65)
axes[1,1].text(.97,.05,r'At 1 Tyr: $\Delta L/L = 4.945\%$',
               transform=axes[1,1].transAxes,ha='right',fontsize=10)
for ax in axes.flat:ax.grid(alpha=.15,lw=.6)
fig.tight_layout(rect=(0,.045,1,1))
fig.text(.075,.015,'0.05820 solar masses. Layered mixing assumes H/d = 25; its transition law remains uncertain.',fontsize=10)
for suffix in ['png','pdf']:fig.savefig(W/f'mesh_comparison.{suffix}',dpi=170)
plt.close(fig)
