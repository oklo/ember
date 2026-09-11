#!/usr/bin/env python3
"""Convection and the final hydrogen profile from the current plotted track."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
from current_evolution import current_history

HERE = Path(__file__).resolve().parent


def main():
    columns, rows = current_history()
    star = {k: np.array([r[i] for r in rows]) for i, k in enumerate(columns)}
    record = json.loads((HERE/'evolution_latest_provenance.json').read_text())
    path = HERE/record['profile_csv']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == record['profile_csv_sha256']
    profile = list(csv.DictReader(path.open()))
    mass = np.array([float(r['mass_g']) for r in profile]); mass /= mass[-1]
    hydrogen = np.array([float(r['X']) for r in profile])
    age = star['age_yr']/1e12
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 2, figsize=(7, 3.25))
    fig.subplots_adjust(left=.10, right=.98, bottom=.20, top=.87, wspace=.37)
    axes[0].plot(age, star['convective_mass_fraction'], color='#ad4d28', lw=1.4)
    axes[0].set(xlim=(3.4, age[-1]+.02), ylim=(0, 1.04), xlabel='Age (trillion yr)',
                ylabel='Convective mass fraction')
    axes[1].plot(mass, hydrogen, color='#ad4d28', lw=1.5)
    boundary = next(r['mass_fraction'] for r in record['mixing_regions'] if not r['convective'])
    axes[1].axvspan(0, boundary, color='#718dab', alpha=.12)
    axes[1].axvline(boundary, color='#718dab', lw=.7, alpha=.7)
    axes[1].set(xlim=(0, 1), ylim=(0, .19), xlabel='Enclosed mass / stellar mass',
                ylabel='Hydrogen mass fraction')
    for ax in axes:
        ax.grid(alpha=.15)
        for axis in [ax.xaxis, ax.yaxis]:
            axis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:.4g}'))
    for suffix in ['pdf', 'png']:
        meta = {'CreationDate': None, 'ModDate': None} if suffix == 'pdf' else {'Software': 'Ember'}
        fig.savefig(HERE/f'convection_current.{suffix}', dpi=180, metadata=meta)
    plt.close(fig)


if __name__ == '__main__':
    main()
