"""Compare the two Ember flash histories using saved surface states."""
from pathlib import Path
import argparse
import csv
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import numpy as np


def plot(data, output):
    meta = json.loads((data/'inputs.json').read_text())
    # Endpoints may have been retained from CSV or appended from native JSON.
    for endpoint in meta['endpoints'].values():
        for key in ['pulse_years', 'Teff_K']:
            endpoint[key] = float(endpoint[key])
    for name,digest in meta['csv_sha256'].items():
        assert hashlib.sha256((data/name).read_bytes()).hexdigest() == digest

    def load(name):
        with (data/name).open() as f:
            return list(csv.DictReader(f))

    main, onset, older, pre = [load(name+'.csv') for name in ('main_sequence','onset','older','pre_flash')]
    def arrays(rows):
        return np.array([float(r['Teff_K']) for r in rows]), np.array([float(r['log_L_Lsun']) for r in rows])
    mt, ml = arrays(main)
    pt, pl = arrays(pre)
    ot, ol = arrays(onset)
    lt, ll = arrays(older)
    orange, purple, grey = '#b65125', '#6a5792', '#797979'
    plt.rcParams.update({'font.size': 11.5, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42})
    fig, (full, zoom) = plt.subplots(1,2,figsize=(12.4,6.3),
                                   gridspec_kw={'width_ratios':[1.23,1]})
    fig.subplots_adjust(left=.075,right=.98,bottom=.19,top=.87,wspace=.25)
    full.plot(mt,ml,color=grey,lw=1.5,zorder=2)
    full.plot(pt,pl,color=grey,lw=1.5,zorder=2)
    full.plot([mt[-1],pt[0]],[ml[-1],pl[0]],color=grey,ls=':',lw=1.35,zorder=2)
    full.plot([pt[-1],ot[0]],[pl[-1],ol[0]],color=grey,ls=':',lw=1.35,zorder=2)
    for ax in (full,zoom):
        ax.plot(lt,ll,color=purple,lw=1.8,alpha=.88,zorder=3)
        ax.plot(ot,ol,color=orange,lw=2.,zorder=4)
        ax.plot(ot[0],ol[0],'o',ms=5,mfc='white',mec='black',mew=.8,zorder=6)
        ax.plot(ot[-1],ol[-1],'o',ms=5,color=orange,zorder=7)
        ax.set_xlabel('Effective surface temperature (K)')
        ax.set_ylabel(r'$\log_{10}(L/L_\odot)$')
        ax.grid(color='0.85',alpha=.5,lw=.55)
        for axis in (ax.xaxis,ax.yaxis):
            axis.set_major_formatter(FuncFormatter(lambda x,_: f'{x:.4g}'))
    full.plot(lt[-1],ll[-1],'o',ms=5,color=purple,zorder=6)
    full.set(xlim=(6700,2500),ylim=(-3.3,-.64))
    zoom.set(xlim=(max(4950,float(ot.max())+35),4500),
             ylim=(-2.87,max(-2.19,float(ol.max())+.12)))
    full.set_xticks([6500,5500,4500,3500,2500])
    tick_spacing = 100 if abs(np.diff(zoom.get_xlim())[0]) < 700 else 200
    zoom.set_xticks(np.arange(4500, np.floor(zoom.get_xlim()[0]/tick_spacing)*tick_spacing+1, tick_spacing))
    def callout(ax,text,point,position,ha='left',relpos=(.5,.5)):
        ax.annotate(text,xy=point,xytext=position,ha=ha,va='center',fontsize=10.2,
                    arrowprops={'arrowstyle':'-|>','color':'0.22','lw':.85,
                                'mutation_scale':11,'shrinkA':4,'shrinkB':4,
                                'relpos':relpos},zorder=9)
    e=meta['endpoints']['older']
    callout(full,f"{e['pulse_years']/1e6:.4g} Myr\n{e['Teff_K']:.4g} K",
            (lt[-1],ll[-1]),(min(6400,lt[-1]+500),min(-2.48,ll[-1]-.28)),
            relpos=(1.,1.))
    peak=int(np.argmax(lt))
    callout(full,f"Temperature maximum\n{lt[peak]:.4g} K",(lt[peak],ll[peak]),(5580,-1.15))
    e=meta['endpoints']['onset']
    unit,scale=('Myr',1e6) if e['pulse_years']>=1e6 else ('kyr',1e3)
    callout(zoom,f"{e['pulse_years']/scale:.4g} {unit}\n{e['Teff_K']:.4g} K",
            (ot[-1],ol[-1]),(max(4630,min(5950,ot[-1]-450)),ol[-1]-.12))
    callout(zoom,'Common starting model',(ot[0],ol[0]),(4520,-2.81),ha='right')
    full.text(2820,-3.18,r'$0.1\,M_\odot$',ha='center',fontsize=11)
    full.annotate('Atmosphere adjustment',xy=((mt[-1]+pt[0])/2,(ml[-1]+pl[0])/2),
                  xytext=(3650,-2.94),fontsize=9.5,ha='center',va='center',
                  arrowprops={'arrowstyle':'-|>','color':'0.35','lw':.7,'mutation_scale':9})

    def direction(ax,rows,age,color,span=.055):
        tx,ly=arrays(rows)
        ages=np.array([float(r['pulse_years']) for r in rows])
        i=int(np.argmin(abs(ages-age)))
        sx=abs(np.diff(ax.get_xlim())[0]); sy=abs(np.diff(ax.get_ylim())[0])
        for j in range(i+1,len(rows)):
            if np.hypot((tx[j]-tx[i])/sx,(ly[j]-ly[i])/sy)>span:
                ax.annotate('',xy=(tx[j],ly[j]),xytext=(tx[i],ly[i]),
                            arrowprops={'arrowstyle':'-|>','color':color,'lw':1.3,
                                        'mutation_scale':12},zorder=8)
                return
    for age in (9e4,7e5,1.8e6):
        direction(full,older,age,purple,.045)
    direction(full,onset,7e5,orange,.025)
    direction(zoom,onset,5000,orange,.09)
    direction(zoom,onset,13000,orange,.06)
    direction(zoom,older,8000,purple,.08)
    direction(zoom,older,18000,purple,.055)
    # A quiet rectangle identifies the enlarged part of the same HR plane.
    from matplotlib.patches import Rectangle
    x0, x1 = sorted(zoom.get_xlim())
    y0, y1 = zoom.get_ylim()
    full.add_patch(Rectangle((x0,y0),x1-x0,y1-y0,fill=False,edgecolor='0.65',lw=.7,ls='--',zorder=1))
    handles=[Line2D([],[],color=grey,lw=1.5,label='Before the flash'),
             Line2D([],[],color=orange,lw=2.,label='Finite mixing from onset'),
             Line2D([],[],color=purple,lw=1.8,label='Alternative mixing history')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.52,.985),
               ncol=3,frameon=False,fontsize=11,handlelength=2.7,columnspacing=2.)
    fig.text(.075,.09,'Right: enlargement of the flash region. Times are measured from the saved model before shell convection.',fontsize=9.5)
    fig.text(.075,.054,'These are alternative histories. Lines join saved models; short adjustments are represented by their endpoints.',fontsize=9.5)
    fig.text(.075,.018,'Dotted connections mark atmosphere or starting-model changes, not resolved evolution. The complete flash remains unverified.',fontsize=9.5)
    output.mkdir(parents=True,exist_ok=True)
    fig.savefig(output/'ember_two_sequences_hr.png',dpi=200)
    fig.savefig(output/'ember_two_sequences_hr.pdf',metadata={'CreationDate':None})
    plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    plot(args.data,args.output)
