#!/usr/bin/env python3
"""Plot the completed convection transition using the published CSV files."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    with (HERE/name).open() as stream:
        rows = list(csv.DictReader(stream))
    return {key: np.array([float(r[key]) if r[key] else np.nan for r in rows]) for key in rows[0]}


def main():
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 2, figsize=(7, 3.15))
    fig.subplots_adjust(left=.095, right=.98, bottom=.18, top=.83, wspace=.35)
    for points, history_name, color, style in [
            (512, 'evolution_history.csv', '#ad4d28', '-'),
            (1024, 'evolution_1024_continuation.csv', '#706078', '--')]:
        history = read(history_name)
        profile = read(f'evolution_profile_{points}.csv')
        assert history['age_yr'][-1] == 3.56e12
        axes[0].plot(history['age_yr']/1e12, history['convective_mass_fraction'],
                     color=color, ls=style, lw=1.3, label=f'{points} mass points')
        axes[1].plot(profile['mass_g']/profile['mass_g'][-1], profile['X'],
                     color=color, ls=style, lw=1.3)
    axes[0].set(xlim=(3.54, 3.56), ylim=(.59, 1.015),
                xlabel='Age (trillion yr)', ylabel='Convective mass fraction')
    axes[1].set(xlim=(0, 1), ylim=(.145, .18),
                xlabel='Enclosed fraction of stellar mass', ylabel='Hydrogen mass fraction',
                title='Interior at 3.560 trillion yr')
    for ax in axes:
        ax.grid(alpha=.15)
        for axis in [ax.xaxis, ax.yaxis]:
            axis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.4g}'))
    fig.legend(*axes[0].get_legend_handles_labels(), loc='upper center', ncol=2, frameon=False)
    fig.savefig(HERE/'convection_transition.pdf', metadata={'CreationDate': None})
    fig.savefig(HERE/'convection_transition.png', dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    main()
