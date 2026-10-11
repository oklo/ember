"""Plot the saved mixing and metallicity comparisons; no stellar tables required."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEST = Path(__file__).resolve().parent


def series():
    groups = [{}, {}]
    with (DEST / 'mixing_metallicity.csv').open() as f:
        for record in csv.DictReader(f):
            panel = int(record['panel'])
            label = record['treatment']
            rows = groups[panel].setdefault(label, [])
            rows.append(dict(years=float(record['age_Gyr'])*1e9,
                             nuclear_fraction=float(record['nuclear_fraction']),
                             Teff_K=float(record['Teff_K'])))
    styles = [[('#777777', '-'), ('#cb722c', '--'), ('#347f91', '-')],
              [('#346b9c', '--'), ('#777777', '-'), ('#a84e49', ':')]]
    result = []
    for group, choices in zip(groups, styles):
        result.append([dict(label=label, rows=sorted(rows, key=lambda r:r['years']),
                            color=color, style=style)
                       for (label, rows), (color, style) in zip(group.items(), choices)])
    return result


if __name__ == '__main__':
    groups = series()
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':12,
                         'axes.labelsize':13, 'legend.fontsize':10.5,
                         'axes.spines.top':False,'axes.spines.right':False})
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.1), sharex=True,
                             gridspec_kw={'height_ratios':[1,1]})
    fig.subplots_adjust(left=.095, right=.98, top=.95, bottom=.14, wspace=.22, hspace=.10)
    latest = max(row['years']/1e9 for group in groups for curve in group for row in curve['rows'])
    upper_age = max(150, latest*1.12)
    visible = [row for group in groups for curve in group for row in curve['rows'] if row['years'] >= 1e10]
    upper_surplus = max(16, 5*np.ceil(max(100*(row['nuclear_fraction']-1) for row in visible)/5)+5)
    lower_surplus = min(-15, 5*np.floor(min(100*(row['nuclear_fraction']-1) for row in visible)/5)-2)
    upper_teff = max(1530, 50*np.ceil(max(row['Teff_K'] for row in visible)*1.02/50))
    for column, curves in enumerate(groups):
        for curve in curves:
            h = curve['rows']
            x = np.array([row['years']/1e9 for row in h])
            mask = x >= 10
            frac = np.array([row['nuclear_fraction'] for row in h])
            teff = np.array([row['Teff_K'] for row in h])
            for row, y in enumerate([100*(frac-1),teff]):
                ax = axes[row,column]
                ax.plot(x[mask],y[mask],color=curve['color'],ls=curve['style'],lw=2.1,
                        label=curve['label'])
                ax.plot(x[-1],y[-1],marker='o',ms=5,mec=curve['color'],
                        mfc=curve['color'] if x[-1]>=1000*(1-1e-10) else 'white',zorder=4)
        axes[0,column].axhline(0,color='.4',lw=.8,ls=':')
        axes[0,column].set_ylim(lower_surplus,upper_surplus)
        axes[0,column].legend(loc='lower left' if column == 0 else 'lower right',frameon=False)
        axes[1,column].set_xlabel('Age (Gyr)')
        axes[1,column].set_ylim(1080,upper_teff)
        for ax in axes[:,column]:
            ax.set_xscale('log')
            ax.set_xlim(10,upper_age)
            ticks = [t for t in [10,20,50,100,200,500,1000] if t <= upper_age]
            ax.set_xticks(ticks, labels=[str(t) for t in ticks])
            ax.grid(alpha=.15)
    axes[0,0].set_ylabel(r'$L_{\rm nuc}/L-1$ (%)')
    axes[1,0].set_ylabel(r'$T_{\rm eff}$ (K)')
    axes[0,0].text(.025,.955,r'$0.05825\,M_\odot$ · $3Z_{\rm ref}$',
                   transform=axes[0,0].transAxes,va='top',fontsize=13)
    axes[0,1].text(.025,.955,r'$0.05815\,M_\odot$ · fixed helium',
                   transform=axes[0,1].transAxes,va='top',fontsize=13)
    fig.text(.095,.061,'Cloudy atmospheres; each sequence starts on the Hayashi track. Filled endpoints: 1 Tyr. Open: incomplete.',fontsize=10.5,color='.3')
    fig.text(.095,.032,'Weak mixing: composition only. Heat + composition: assumed layer height H/d = 25 and marginal mixing between regimes.',fontsize=10.5,color='.3')
    for ext in ['png','pdf']:
        fig.savefig(DEST / f'mixing_metallicity.{ext}',dpi=200)
