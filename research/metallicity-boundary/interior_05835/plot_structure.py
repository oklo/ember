"""Plot the retained 0.05835-solar-mass track and its verified interior replay."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

W = Path(__file__).resolve().parent


def main():
    z = np.load(W/'interior_05835.npz')
    meta = json.loads((W/'provenance.json').read_text())
    history = [json.loads(s) for s in (W/'original_history.jsonl').read_text().splitlines()]
    whole_t = np.array([s['years']/1e9 for s in history])
    t, q = z['age_Gyr'], z['q']
    edges = np.r_[0., .5*(q[:-1]+q[1:]), 1.]
    conv = z['convective']
    kind = np.where(conv, 0, np.where(z['conductive_share']>=.5, 1, 2))
    palette = ListedColormap(['#b9d2e6', '#eed0a5', '#fffdfa'])
    norm = BoundaryNorm([-.5,.5,1.5,2.5],3)
    # Never paint across unavailable time intervals, including partial previews.
    breaks = (np.r_[0, np.flatnonzero(np.diff(t)>.4)+1, len(t)] if meta['partial']
              else np.array([0,len(t)]))
    chunks = [np.arange(a,b) for a,b in zip(breaks[:-1],breaks[1:]) if b-a>1]
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':12.5,
        'axes.labelsize':14, 'axes.spines.top':False, 'axes.spines.right':False,
        'savefig.dpi':190, 'hatch.linewidth':.6})
    fig, axes = plt.subplots(3,2,figsize=(14.5,10.8),sharex='col',sharey='row',
        gridspec_kw={'width_ratios':[1.75,1], 'height_ratios':[.92,1.25,1.25]})
    fig.subplots_adjust(left=.085,right=.90,bottom=.09,top=.88,wspace=.14,hspace=.105)
    limits = [(0,897),(850,897)]
    for col,(lo,hi) in enumerate(limits):
        top, transport, composition = axes[:,col]
        top.plot(whole_t,[s['nuclear_fraction']for s in history],color='#242424',lw=2)
        top.axhline(1,color='.35',ls=':',lw=.85)
        top.set(ylim=(0,1.15),yticks=[0,.5,1])
        top.grid(axis='y',alpha=.13)
        thermal=top.twinx()
        thermal.plot(whole_t,np.array([s['central_T_K']for s in history])/1e6,
                     color='#ab5128',ls='--',lw=1.65)
        thermal.set(ylim=(0,3.3),yticks=[0,1,2,3])
        thermal.spines['right'].set_visible(col==1)
        if col==0:thermal.set_yticklabels([]);thermal.tick_params(right=False)
        else:thermal.set_ylabel('Central temperature (MK)',color='black',labelpad=9)
        for ix in chunks:
            if t[ix[-1]]<lo or t[ix[0]]>hi:continue
            x=t[ix]
            te=np.r_[x[0],.5*(x[:-1]+x[1:]),x[-1]]
            transport.pcolormesh(te,edges,kind[ix].T,cmap=palette,norm=norm,
                                shading='flat',rasterized=True)
            transport.contour(x,q,conv[ix].T.astype(float),levels=[.5],
                              colors='#3b4c59',linewidths=.7)
            contours=transport.contour(x,q,z['T_K'][ix].T/1e6,
                levels=[.5,1.,2.],colors='#76533a',linestyles='--',linewidths=.85)
            # Label temperature contours in the broad view near 600 Gyr, and
            # near the cooling endpoint in the enlargement, using actual profiles.
            target=600 if col==0 else 887
            if x[0]<=target<=x[-1]:
                row=ix[np.argmin(abs(x-target))];positions=[]
                for level in [.5,1,2]:
                    temp=z['T_K'][row]/1e6
                    if temp.min()<level<temp.max():
                        positions.append((t[row],np.interp(level,temp[::-1],q[::-1])))
                transport.clabel(contours,fmt=lambda v:f'{v:g} MK',manual=positions,
                                  fontsize=10.5,inline_spacing=4)
            im=composition.pcolormesh(te,edges,z['X'][ix].T,cmap='viridis',
                vmin=.42,vmax=.80,shading='flat',rasterized=True)
            composition.contour(x,q,conv[ix].T.astype(float),levels=[.5],
                                colors='white',linewidths=.7,alpha=.8)
            metals=composition.contour(x,q,z['Z'][ix].T,levels=[.075,.10,.15],
                                       colors='black',linestyles=':',linewidths=1.05)
            positions=[]
            if col==0:
                for level,target in [(.075,220),(.10,490),(.15,870)]:
                    if not x[0]<=target<=x[-1]:continue
                    row=ix[np.argmin(abs(x-target))]
                    profile=z['Z'][row]
                    crossings=np.flatnonzero((profile[:-1]-level)*(profile[1:]-level)<0)
                    if len(crossings):
                        j=crossings[0]
                        pos=q[j]+(q[j+1]-q[j])*(level-profile[j])/(profile[j+1]-profile[j])
                        positions.append((t[row],pos))
                composition.clabel(metals,fmt=lambda v:f'Z={v:g}',manual=positions,
                                    fontsize=10.5,inline_spacing=4)
        for panel in [transport,composition]:
            panel.set(ylim=(0,1),yticks=[0,.25,.5,.75,1])
            panel.tick_params(axis='both',labelsize=12)
        composition.set_xlabel('Age (Gyr)')
        top.set_xlim(lo,hi)
        if col==0:composition.set_xticks([0,200,400,600,800])
        else:composition.set_xticks([850,860,870,880,890,897])
    axes[0,0].set_ylabel(r'$L_{\rm nuc}/L$')
    axes[1,0].set_ylabel(r'Enclosed mass $m/M$')
    axes[2,0].set_ylabel(r'Enclosed mass $m/M$')
    arrow=dict(arrowstyle='-|>',lw=.85,color='.2',mutation_scale=11,shrinkA=4,shrinkB=3)
    axes[0,0].annotate('99.72%',xy=(363.6516,.997188),xytext=(230,.67),
        fontsize=12,arrowprops=arrow)
    axes[0,0].text(.50,.34,'Core warms, then cools',transform=axes[0,0].transAxes,fontsize=11.5)
    axes[0,1].annotate('Burning fades',xy=(880,.1498),xytext=(860,.45),
        fontsize=12,arrowprops=arrow)
    axes[1,0].text(.055,.78,'Convective envelope',transform=axes[1,0].transAxes,fontsize=12)
    axes[1,0].text(.055,.12,'Convective core',transform=axes[1,0].transAxes,fontsize=12)
    axes[1,0].annotate('Conductive separating layer',xy=(450,.331),xytext=(100,.44),
        fontsize=11.5,arrowprops=arrow)
    axes[1,1].text(.62,.18,'Conductive\ninterior',transform=axes[1,1].transAxes,
                   ha='center',fontsize=12)
    axes[2,0].text(.03,.92,'Hydrogen',transform=axes[2,0].transAxes,fontsize=12,
        bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=3))
    # The explicit cloud rule is marked on the zoom, not presented as a
    # demonstrated cause of the preceding centuries of evolution.
    for ax in axes[:,1]:
        for age in [875.647,891.863]:ax.axvline(age,color='.4',lw=.7,ls=':',alpha=.7)
    axes[0,1].text(.97,.97,'Cloud clearing: 1200 → 900 K',ha='right',va='top',
        transform=axes[0,1].transAxes,fontsize=10)
    axes[2,1].text(.07,.12,'Center at 897 Gyr\nX = 0.4333 · Z = 0.1838',
        transform=axes[2,1].transAxes,fontsize=11,
        bbox=dict(facecolor='white',edgecolor='none',alpha=.88,pad=4))
    cax=fig.add_axes([.93,.10,.014,.275])
    bar=fig.colorbar(im,cax=cax);bar.set_label('Hydrogen mass fraction X',fontsize=12)
    bar.set_ticks([.45,.55,.65,.75]);bar.ax.tick_params(labelsize=11)
    fig.text(.085,.96,r'$0.05835\,M_\odot$ · $Z=0.04945$ · cloudy, $f_{\rm sed}=2$',fontsize=17)
    fig.legend(handles=[Line2D([],[],color='#242424',lw=2,label=r'$L_{\rm nuc}/L$'),
        Line2D([],[],color='#ab5128',ls='--',lw=1.65,label='Central temperature'),
        Patch(facecolor=palette(0),edgecolor='#3b4c59',lw=.6,label='Convection'),
        Patch(facecolor=palette(1),label='Nonconvective: conduction'),
        Patch(facecolor=palette(2),edgecolor='.6',lw=.6,label='Nonconvective: radiation')],
        loc='upper left',bbox_to_anchor=(.08,.943),ncol=5,frameon=False,fontsize=10.9,
        columnspacing=1.3,handlelength=1.7)
    fig.text(.085,.025,'Neutral-metal core EOS; instantaneous convection without semiconvection; prescribed cloud clearing.',fontsize=11)
    fig.text(.607,.893,'Late evolution, enlarged',fontsize=12)
    if meta['partial']:
        fig.text(.5,.48,'PREVIEW: profile reconstruction in progress',ha='center',fontsize=12,
                 bbox=dict(facecolor='white',edgecolor='.7',alpha=.9,pad=5))
    for extension in ['png','pdf']:
        fig.savefig(W/f'Ember_boundary_05835_structure.{extension}')
    plt.close(fig)


if __name__ == '__main__':
    main()
