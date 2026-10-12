"""Metallicity boundary with measured mass intervals and unresolved trials."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator,FuncFormatter,NullFormatter

OUT=Path(__file__).resolve().parent

def main():
    data=json.loads((OUT/'transition_mass.json').read_text())
    bounds=data['brackets']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,
        'axes.labelsize':14,'axes.spines.top':False,'axes.spines.right':False,
        'pdf.fonttype':42,'savefig.facecolor':'white'})
    fig,ax=plt.subplots(figsize=(10.6,7.5))
    fig.subplots_adjust(left=.11,right=.965,bottom=.225,top=.90)
    colors={'nc':'#355f8b','f2':'#b76b32'}
    for cloud in ['nc','f2']:
        bs=[b for b in bounds if b['cloud']==cloud];color=colors[cloud]
        x=[];y=[]
        for b in bs:
            z=b['Z'];lo=b['lower_cooling_mass'];hi=b['upper_burning_mass']
            if lo is not None and hi is not None:
                mid=(lo+hi)/2
                _,_,bars=ax.errorbar(z,mid,yerr=[[mid-lo],[hi-mid]],fmt='none',
                            color=color,alpha=.45,elinewidth=1.3,capsize=6,capthick=1.1,zorder=3)
                if b['unmatched_eos_or_envelope']:
                    for bar in bars:bar.set_linestyle('--')
                x.append(z);y.append(mid)
                label=f'{lo:.4g}–{hi:.4g}'+('*' if b.get('extended_tables') or b['unmatched_eos_or_envelope'] else '')
            elif hi is not None:
                ax.errorbar(z,hi,yerr=.0012,uplims=True,fmt='none',
                            color=color,alpha=.45,elinewidth=1.3,capsize=5,zorder=3)
                label=f'≤ {hi:.4g}'
            else:
                continue
            offset=(0,11) if lo is not None else (10,-2)
            if b['metallicity']=='z100':offset=(0,-20)
            if b['metallicity']=='z033':offset=(0,-23)
            if b['metallicity']=='z058' and cloud=='nc':offset=(0,-23)
            label_mass=hi
            if b['metallicity']=='z058' and cloud=='f2':
                label_mass=lo;offset=(0,-20)
            if b['metallicity']=='z173':
                label_mass=lo;offset=(0,-20)
            if b.get('extended_tables'):offset=(0,-22)
            ax.annotate(label,(z,label_mass),xytext=offset,textcoords='offset points',
                        ha='center' if lo is not None else 'left',fontsize=10.5,color='.4')
            for pending in b['unresolved']+b.get('unresolved_below_cooling',[]):
                ax.plot(z,pending['mass_Msun'],marker='D',ms=5,mfc='white',
                        mec=color,mew=1.1,alpha=.5,zorder=5)
        ax.plot(x,y,color=color,lw=1,ls=':',alpha=.35,zorder=1)
    paired={(b['Z'],b['cloud']):b for b in data.get('corrected_brackets',[])}
    for b in paired.values():
        ax.plot([b['Z']]*2,[b['lower_cooling_mass'],b['upper_burning_mass']],
                color=colors[b['cloud']],lw=3,zorder=6)
        for pending in b.get('unresolved',[]):
            ax.plot(b['Z'],pending['mass_Msun'],marker='D',ms=5.5,
                    mfc='white',mec=colors[b['cloud']],mew=1.3,zorder=8)
    for model in data['corrected_burning_models']:
        assert model['status']=='supported_1Tyr' and not model['included_in_cutoff']
        color=colors[model['cloud']]
        ax.plot(model['Z'],model['mass_Msun'],marker='s',ms=7.5,
                mfc=color,mec='white',mew=.8,zorder=7)
        solar=model['Z']<.02
        pair=paired.get((model['Z'],model['cloud']))
        label=(f"{pair['lower_cooling_mass']:.4g}–{pair['upper_burning_mass']:.4g}"
               if pair else f"{model['mass_Msun']:.4g}")
        low_z=model['Z']<.007
        ax.annotate(label,
                    (model['Z'],model['mass_Msun']),xytext=(10 if low_z else -10,7 if solar else 0),
                    textcoords='offset points',ha='left' if low_z else 'right',
                    va='bottom' if solar else 'center',fontsize=10.5)
    for model in data.get('corrected_cooling_models',[]):
        assert model['status']=='cooling' and not model['included_in_cutoff']
        color=colors[model['cloud']]
        ax.plot(model['Z'],model['mass_Msun'],marker='o',ms=7.5,
                mfc='white',mec=color,mew=1.7,zorder=7)
        if (model['Z'],model['cloud']) not in paired:
            low_z=model['Z']<.007
            ax.annotate(f"{model['mass_Msun']:.4g}",
                        (model['Z'],model['mass_Msun']),xytext=(10 if low_z else -10,0),
                        textcoords='offset points',ha='left' if low_z else 'right',va='center',fontsize=10.5)
    for model in data.get('corrected_unresolved_models',[]):
        assert not model['included_in_cutoff']
        ax.plot(model['Z'],model['mass_Msun'],marker='D',ms=6,
                mfc='white',mec=colors[model['cloud']],mew=1.4,zorder=8)
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
            Line2D([],[],color='.3',marker='D',mfc='white',ls='none',ms=5,label='Unfinished mass trial'),
            Line2D([],[],color='.3',marker='s',ls='none',ms=6,
                   label='Metal ionization: burning at 1 Tyr'),
            Line2D([],[],color='.3',marker='o',mfc='white',ls='none',ms=6,
                   label='Metal ionization: cooling, nuclear fraction < 0.1%')]
    ax.legend(handles=legend,loc='lower left',bbox_to_anchor=(.015,.105),frameon=False,fontsize=11)
    fig.text(.11,.95,'Ember  |  Hydrogen-burning boundary',fontsize=17)
    fig.text(.11,.135,'Thin bars: neutral metals. Thick bars: metal ionization. Tested mass intervals, not statistical errors.',fontsize=10.3)
    fig.text(.11,.098,'* Provisional high-Z pair: 0.05835 cools by 897 Gyr; lower-mass warm tracks remain unfinished.',fontsize=10.3)
    fig.text(.11,.061,'Squares burn at 1 Tyr; circles cool below 0.1% nuclear support. Nearby masses remain under test.',fontsize=10.3)
    for ext in ['png','pdf']:
        fig.savefig(OUT/f'transition_mass.{ext}',dpi=220)
    plt.close(fig)

if __name__=='__main__':main()
