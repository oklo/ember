#!/usr/bin/env python3
"""Plot the recovered photosphere over a fixed hydrogen-rich opacity map."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter
import numpy as np

HERE = Path(__file__).resolve().parent
EMBER = '#ad431f'


def main():
    record = json.loads((HERE/'photospheric_evolution_provenance.json').read_text())
    for name, key in [('photospheric_evolution.csv', 'photosphere_csv_sha256'),
                      ('photospheric_opacity_map.json', 'opacity_map_sha256'),
                      ('evolution_latest.csv', 'history_sha256')]:
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == record[key]
    rows = list(csv.DictReader((HERE/'photospheric_evolution.csv').open()))
    star = {key: np.array([float(r[key]) for r in rows]) for key in rows[0]}
    source = json.loads((HERE/'photospheric_opacity_map.json').read_text())
    x = np.log10(star['photosphere_rho_g_cm3'])
    y = np.log10(star['photosphere_T_K'])
    radius, age = star['R_Rsun'], star['age_yr']/1e12
    # Sparse actual accepted models, separated along the displayed curve.
    distance = np.r_[0, np.cumsum(np.hypot(np.diff(x)/np.ptp(x), np.diff(y)/np.ptp(y)))]
    chosen = np.unique(np.searchsorted(distance, np.linspace(0, distance[-1], 22)))
    xlim = (x.min()-.065, x.max()+.065)
    ylim = (y.min()-.026, y.max()+.026)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.5))
    fig.subplots_adjust(left=.075, right=.875, bottom=.22, top=.96, wspace=.29)
    for ax in axes:
        shade = ax.pcolormesh(source['log10_density_g_cm3'], source['log10_temperature_K'],
                             source['log10_kappa_cm2_g'], shading='gouraud',
                             cmap='Greys', vmin=-4, vmax=2, rasterized=True)
        ax.plot(x, y, color='white', lw=3.0, zorder=3)
        ax.plot(x, y, color=EMBER, lw=1.1, zorder=4)
        # scatter area is in squared points: its circle diameter scales as R.
        ax.scatter(x[chosen], y[chosen], s=(12*radius[chosen]/.14)**2,
                   facecolors='none', edgecolors=EMBER, linewidths=.9, zorder=5)
        ax.set(xlabel=r'$\log_{10}\,[\rho_{\rm ph}/(\mathrm{g\,cm^{-3}})]$',
               ylabel=r'$\log_{10}\,(T_{\rm ph}/\mathrm{K})$')
        for axis in [ax.xaxis, ax.yaxis]:
            axis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:.4g}'))
    axes[0].set(xlim=(-8, -2), ylim=(3, 4))
    # Show the LBA97 axis extent without inventing opacity beyond the source table.
    axes[0].set_facecolor('#f4f4f4')
    xmax = max(source['log10_density_g_cm3'])
    ymin, ymax = min(source['log10_temperature_K']), max(source['log10_temperature_K'])
    for left, bottom, width, height in [(xmax, 3, -2-xmax, 1),
                                        (-8, 3, xmax+8, ymin-3),
                                        (-8, ymax, xmax+8, 4-ymax)]:
        axes[0].add_patch(Rectangle((left, bottom), width, height, facecolor='none',
                                   edgecolor='#b8b8b8', hatch='///', lw=0, zorder=2))
    axes[0].add_patch(Rectangle((xlim[0], ylim[0]), xlim[1]-xlim[0], ylim[1]-ylim[0],
                               fill=False, edgecolor='#a66048', linestyle='--', lw=.8))
    axes[0].annotate('Enlarged at right', (xlim[1], ylim[1]), xytext=(-5.9, 3.79),
                     fontsize=9, arrowprops={'arrowstyle': '-', 'lw': .7, 'color': '#805040'})
    axes[1].set(xlim=xlim, ylim=ylim)
    marks = [0, int(np.argmin(abs(age-.71))), int(np.argmin(abs(age-2))),
             int(np.argmin(abs(age-3.5))), len(age)-1]
    offsets = [(-8, -20), (14, 0), (12, -2), (12, -3), (11, 9)]
    for i, offset in zip(marks, offsets):
        label = 'Start' if i == 0 else f'{age[i]:.4g} trillion yr'
        axes[1].annotate(label, (x[i], y[i]), xytext=offset, textcoords='offset points',
                         fontsize=8, ha='right' if i == 0 else 'left',
                         arrowprops={'arrowstyle': '-', 'lw': .6, 'color': EMBER})
    cax = fig.add_axes([.904, .22, .018, .74])
    colorbar = fig.colorbar(shade, cax=cax)
    colorbar.set_label(r'$\log_{10}\,[\kappa_{\rm R,abs}/(\mathrm{cm^2\,g^{-1}})]$')
    handles = [Line2D([], [], color=EMBER, marker='o', markerfacecolor='none',
                      markersize=12*r/.14, linestyle='none', label=f'{r:.2g} '+r'$R_\odot$')
               for r in [.12, .14, .16]]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.43, .025), ncol=3,
               frameon=False)
    for suffix in ['pdf', 'png']:
        metadata = {'Creator': 'Ember', 'CreationDate': None, 'ModDate': None} if suffix == 'pdf' else {'Software': 'Ember'}
        fig.savefig(HERE/f'photospheric_evolution.{suffix}', dpi=200, metadata=metadata)
    plt.close(fig)
    print(f'Plotted {len(rows)} photospheric states, {len(chosen)} radius-scaled circles.')


if __name__ == '__main__':
    main()
