"""Independent meshes; lines connect accepted states without smoothing."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

W=Path(__file__).resolve().parent
plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,1,figsize=(9,7),sharex=True)
for n,color in [(536,'#c47733'),(1071,'#347f91')]:
    a=np.genfromtxt(W/f'mesh_{n}.csv',delimiter=',',names=True)
    for ax,key in zip(axes,['nuclear_fraction','luminosity_Lsun']):
        y=a[key] if key=='nuclear_fraction' else np.log10(a[key])
        ax.plot(a['age_Gyr'],y,color=color,lw=1.5,label=f'{n} points')
        ax.plot(a['age_Gyr'][-1],y[-1],'o',color=color,ms=5)
axes[0].axhline(1,color='.5',ls=':',lw=1)
axes[0].set_ylabel(r'$L_{\rm nuc}/L$');axes[0].legend(frameon=False,loc='lower left')
axes[1].set(ylabel=r'$\log_{10}(L/L_\odot)$',xlabel='Age (Gyr)',xlim=(20,122))
for ax in axes:ax.grid(alpha=.15,lw=.6)
fig.tight_layout()
for suffix in ['png','pdf']:fig.savefig(W/f'mesh_comparison.{suffix}',dpi=170)
plt.close(fig)
