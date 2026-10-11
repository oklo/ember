"""Metallicity boundary with measured mass intervals and unresolved trials."""
from pathlib import Path
from datetime import datetime,timezone
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator,FuncFormatter,NullFormatter

OUT=Path(__file__).resolve().parent

def main():
    bounds=json.loads((OUT/'transition_mass.json').read_text())['brackets']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,
        'axes.labelsize':14,'axes.spines.top':False,'axes.spines.right':False,
        'pdf.fonttype':42,'savefig.facecolor':'white'})
    fig,ax=plt.subplots(figsize=(10,7.1))
    fig.subplots_adjust(left=.12,right=.97,bottom=.19,top=.91)
    colors={'nc':'#355f8b','f2':'#b76b32'}
    for cloud in ['nc','f2']:
        bs=[b for b in bounds if b['cloud']==cloud];color=colors[cloud]
        x=[];y=[]
        for b in bs:
            z=b['Z'];lo=b['lower_cooling_mass'];hi=b['upper_burning_mass']
            if lo is not None and hi is not None:
                mid=(lo+hi)/2
                _,_,bars=ax.errorbar(z,mid,yerr=[[mid-lo],[hi-mid]],fmt='none',
                            color=color,elinewidth=2.5,capsize=6,capthick=1.5,zorder=3)
                if b['unmatched_eos_or_envelope']:
                    for bar in bars:bar.set_linestyle('--')
                x.append(z);y.append(mid)
                label=f'{lo:.4g}–{hi:.4g}'+('*' if b['unmatched_eos_or_envelope'] else '')
            elif hi is not None:
                ax.errorbar(z,hi,yerr=.0012,uplims=True,fmt='none',
                            color=color,elinewidth=1.8,capsize=5,zorder=3)
                label=f'≤ {hi:.4g}'
            else:
                continue
            ax.annotate(label,(z,hi),xytext=(0,11),textcoords='offset points',
                        ha='center',fontsize=10.5,color='black')
            for pending in b['unresolved']:
                ax.plot(z,pending['mass_Msun'],marker='D',ms=5,mfc='white',
                        mec=color,mew=1.1,zorder=5)
        ax.plot(x,y,color=color,lw=1,ls=':',alpha=.55,zorder=1)
    ax.set(xscale='log',xlim=(.0045,.06),ylim=(.0558,.0785),
           xlabel='Initial metal mass fraction, Z',ylabel=r'Transition mass ($M_\odot$)')
    zvals=sorted({b['Z'] for b in bounds})
    ax.xaxis.set_major_locator(FixedLocator(zvals))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda z,p:f'{z:.4g}'))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.grid(axis='y',alpha=.14,lw=.65)
    ax.text(.03,.035,r'$Y=0.27$ • solar metal ratios',transform=ax.transAxes,fontsize=11)
    legend=[Line2D([],[],color=colors['nc'],lw=2.5,label='Clear'),
            Line2D([],[],color=colors['f2'],lw=2.5,label=r'Cloudy, $f_{\rm sed}=2$'),
            Line2D([],[],color='.3',marker='D',mfc='white',ls='none',ms=5,label='Unfinished mass trial')]
    ax.legend(handles=legend,loc='lower left',bbox_to_anchor=(.015,.105),frameon=False,fontsize=11)
    fig.text(.12,.958,'Ember  |  Hydrogen-burning boundary',fontsize=17)
    fig.text(.12,.10,'Bars span demonstrated cooling and sustained burning at 1 Tyr; they are not statistical error bars.',fontsize=10.5)
    fig.text(.12,.067,'* EOS/envelope settings differ at interval ends; matched repeats pending. Physical uncertainties are separate.',fontsize=10)
    for ext in ['png','pdf']:
        fig.savefig(OUT/f'transition_mass.{ext}',dpi=220)
    plt.close(fig)

if __name__=='__main__':main()
