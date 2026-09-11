#!/usr/bin/env python3
"""Compare the LBA97 and Ember photospheres over maps on one opacity scale."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter
import numpy as np

HERE = Path(__file__).resolve().parent
EMBER, LBA = '#ad431f', '#357ba3'


def sampled(x, y, count):
    distance = np.r_[0, np.cumsum(np.hypot(np.diff(x)/np.ptp(x), np.diff(y)/np.ptp(y)))]
    return np.unique(np.searchsorted(distance, np.linspace(0, distance[-1], count)))


def main():
    record = json.loads((HERE/'photospheric_evolution_provenance.json').read_text())
    for name, key in [('photospheric_evolution.csv', 'photosphere_csv_sha256'),
                      ('photospheric_opacity_map.json', 'opacity_map_sha256'),
                      ('evolution_latest.csv', 'history_sha256')]:
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == record[key]
    rows = list(csv.DictReader((HERE/'photospheric_evolution.csv').open()))
    star = {key: np.array([float(r[key]) for r in rows]) for key in rows[0]}
    source = json.loads((HERE/'photospheric_opacity_map.json').read_text())
    historical = json.loads((HERE/'lba97_figure6_extraction.json').read_text())
    x, y = np.log10(star['photosphere_rho_g_cm3']), np.log10(star['photosphere_T_K'])
    radius, age = star['R_Rsun'], star['age_yr']/1e12
    chosen = sampled(x, y, 22)
    track = historical['track']
    hx, hy, hr = [np.array(track[k]) for k in ['log10_density_g_cm3', 'log10_temperature_K', 'R_Rsun']]
    hchosen = sampled(hx, hy, 48)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                        'axes.spines.right': False, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=.075, right=.875, bottom=.24, top=.96, wspace=.15)
    norm, cmap = Normalize(-3, 3), 'Greys'
    patches = []
    for cell in historical['cells']:
        a, b = cell['logrho']; c, d = cell['logT']
        patches.append(Rectangle((a, c), b-a, d-c))
    background = PatchCollection(patches, cmap=cmap, norm=norm, linewidth=0, rasterized=True)
    background.set_array(np.array([c['estimated_log10_kappa_cm2_g'] for c in historical['cells']]))
    axes[0].add_collection(background)
    shade = axes[1].pcolormesh(source['log10_density_g_cm3'], source['log10_temperature_K'],
                              source['log10_kappa_cm2_g'], shading='gouraud',
                              cmap=cmap, norm=norm, rasterized=True)
    # Keep the full published axes; missing modern source coverage is explicit.
    xmax = max(source['log10_density_g_cm3'])
    ymin, ymax = min(source['log10_temperature_K']), max(source['log10_temperature_K'])
    axes[1].set_facecolor('#f4f4f4')
    for left, bottom, width, height in [(xmax, 3, -2-xmax, 1),
                                        (-8, 3, xmax+8, ymin-3),
                                        (-8, ymax, xmax+8, 4-ymax)]:
        axes[1].add_patch(Rectangle((left, bottom), width, height, facecolor='none',
                                   edgecolor='#b8b8b8', hatch='///', lw=0, zorder=2))
    for i, ax in enumerate(axes):
        ax.plot(hx, hy, color=LBA, lw=.7, alpha=.48, zorder=3)
        ax.scatter(hx[hchosen], hy[hchosen], s=(11*hr[hchosen]/.14)**2,
                   facecolors='none', edgecolors=LBA, linewidths=.6, alpha=.5, zorder=4)
        ax.plot(x, y, color='white', lw=2.7, zorder=5)
        ax.plot(x, y, color=EMBER, lw=1.1, zorder=6)
        ax.scatter(x[chosen], y[chosen], s=(11*radius[chosen]/.14)**2,
                   facecolors='none', edgecolors=EMBER, linewidths=.9, zorder=7)
        ax.set(xlim=(-8, -2), ylim=(3, 4),
               xlabel=r'$\log_{10}\,[\rho_{\rm ph}/(\mathrm{g\,cm^{-3}})]$')
        ax.text(.035, .955, '('+'ab'[i]+')', transform=ax.transAxes,
                va='top', bbox={'facecolor':'white', 'edgecolor':'none', 'alpha':.75, 'pad':2})
        ax.text(-7.8, 3.065, 'Grains included' if i == 0 else 'Grains not yet included',
                fontsize=9, bbox={'facecolor':'white', 'edgecolor':'none', 'alpha':.8, 'pad':3})
        ax.annotate(f'{age[-1]:.4g} trillion yr', (x[-1], y[-1]), xytext=(-4.15, 3.69),
                    fontsize=8, color=EMBER, arrowprops={'arrowstyle':'-', 'lw':.7, 'color':EMBER})
        for axis in [ax.xaxis, ax.yaxis]:
            axis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:.4g}'))
    axes[0].set_ylabel(r'$\log_{10}\,(T_{\rm ph}/\mathrm{K})$')
    cax = fig.add_axes([.904, .24, .018, .72])
    colorbar = fig.colorbar(shade, cax=cax, extend='both')
    colorbar.set_label(r'$\log_{10}\,[\kappa_{\rm R}/(\mathrm{cm^2\,g^{-1}})]$')
    handles = [Line2D([], [], color=LBA, lw=.8, alpha=.5, label=r'LBA97, $0.1\,M_\odot$'),
               Line2D([], [], color=EMBER, lw=1.1, label=r'Ember, $0.1\,M_\odot$')]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.45, .105), ncol=2, frameon=False)
    sizes = [Line2D([], [], color='#555555', marker='o', markerfacecolor='none',
                    markersize=11*r/.14, linestyle='none', label=f'{r:.2g} '+r'$R_\odot$')
             for r in [.06, .14, .28]]
    fig.legend(handles=sizes, loc='lower center', bbox_to_anchor=(.45, .015), ncol=3, frameon=False)
    for suffix in ['pdf', 'png']:
        metadata = {'Creator':'Ember', 'CreationDate':None, 'ModDate':None} if suffix == 'pdf' else {'Software':'Ember'}
        fig.savefig(HERE/f'photospheric_evolution.{suffix}', dpi=200, metadata=metadata)
    plt.close(fig)
    print(f'Plotted {len(rows)} Ember states and {len(hx)} recovered LBA97 markers with one opacity normalization.')


if __name__ == '__main__':
    main()
