"""Plot the frozen MESA comparison histories beside the retained Ember track."""
from pathlib import Path
import argparse
import csv
import gzip
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import numpy as np


def plot(directory):
    ember = [json.loads(line) for line in gzip.decompress(
        (directory/'metal_history.jsonl.gz').read_bytes()).decode().splitlines()]
    def read_csv(name):
        with (directory/name).open() as stream:
            return [{k: float(v) for k, v in r.items() if v != ''}
                    for r in csv.DictReader(stream)]
    mesa = read_csv('mesa_grey_no_diffusion.csv')
    diffusion = read_csv('mesa_grey_diffusion.csv')
    pp = read_csv('pp_comparison_history.csv')
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42})
    fig, (hr, fuel) = plt.subplots(1, 2, figsize=(7.5, 4.2))
    fig.subplots_adjust(left=.105, right=.985, bottom=.16, top=.78, wspace=.35)
    orange, blue, gray = '#ad4d28', '#326b91', '#777777'
    hr.plot([r['Teff_K'] for r in pp], np.log10([r['L_Lsun'] for r in pp]),
            color=gray, alpha=.38, lw=1, ls='--')
    hr.plot([r['Teff'] for r in ember], np.log10([r['luminosity']/3.828e33 for r in ember]),
            color=orange, lw=1.8, zorder=4)
    for rows, ls, alpha in [(mesa, '-', .9), (diffusion, ':', .7)]:
        hr.plot([10**r['log_Teff'] for r in rows], [r['log_L'] for r in rows],
                color=blue, lw=1.35, ls=ls, alpha=alpha)
        fuel.plot([r['star_age']/1e12 for r in rows],
                  [r['total_mass_h1'] for r in rows], color=blue, lw=1.35, ls=ls, alpha=alpha)
    def hydrogen_mass(row):
        if 'global_species' in row:
            return .1*row['global_species'][0]
        # The common prefix is homogeneous; its older rows did not store
        # a separate whole-star inventory. Never apply this to stratified rows.
        assert abs(row['mixed_mass_fraction']-1) < 1e-12
        assert row['central_H'] == row['surface_H']
        return .1*row['central_H']
    fuel.plot([r['age_years']/1e12 for r in ember],
              [hydrogen_mass(r) for r in ember], color=orange, lw=1.8)
    peak = max(mesa, key=lambda r: r['log_Teff'])
    hr.plot(10**peak['log_Teff'], peak['log_L'], 'o', ms=3.5, color=blue)
    hr.annotate(f"{10**peak['log_Teff']:.4g} K", xy=(10**peak['log_Teff'], peak['log_L']),
                xytext=(5700, -2.85), fontsize=8, color=blue,
                arrowprops=dict(arrowstyle='-', color=blue, lw=.7))
    last = mesa[-1]
    hr.plot(10**last['log_Teff'], last['log_L'], 'o', ms=3.5, color=blue)
    hr.annotate(f"{10**last['log_Teff']:.4g} K", xy=(10**last['log_Teff'], last['log_L']),
                xytext=(2650, -7.05), fontsize=8, color=blue,
                arrowprops=dict(arrowstyle='-', color=blue, lw=.7))
    hr.set(xlim=(7200, 450), ylim=(-7.45, -1.85), xlabel='Effective surface temperature (K)',
           ylabel=r'$\log_{10}(L/L_\odot)$')
    fuel.set(xlim=(0, 4.05), ylim=(0, .072), xlabel='Age (Tyr)',
             ylabel=r'Remaining hydrogen mass ($M_\odot$)')
    for ax in (hr, fuel):
        ax.grid(alpha=.13)
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.4g}'))
    handles = [Line2D([], [], color=orange, lw=1.8, label='Ember: settling metals'),
               Line2D([], [], color=blue, lw=1.35, label='MESA: grey, no diffusion'),
               Line2D([], [], color=gray, alpha=.38, lw=1, ls='--', label='Ember pp comparison (HR only)'),
               Line2D([], [], color=blue, alpha=.7, lw=1.35, ls=':', label='MESA: grey, with diffusion')]
    fig.legend(handles=handles, loc='upper center', ncol=2, frameon=False, fontsize=8.5)
    fig.savefig(directory/'mesa_comparison.pdf', metadata={'CreationDate': None})
    fig.savefig(directory/'mesa_comparison.png', dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    plot(parser.parse_args().directory)
