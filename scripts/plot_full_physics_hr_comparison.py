"""Plot four frozen H-R histories without the source-computation map.

The MESA history is the tabulated-atmosphere/diffusion configuration, not
the separate grey track that reaches remnant cooling. Its limited temporal
coverage is stated on the figure and in the accompanying caption.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np


def plot(snapshot, output):
    meta = json.loads((snapshot / 'inputs.json').read_text())
    for name, digest in meta['snapshot_sha256'].items():
        assert hashlib.sha256((snapshot / name).read_bytes()).hexdigest() == digest
    def rows(name):
        with (snapshot / name).open() as stream:
            return [{k: float(v) for k, v in row.items()}
                    for row in csv.DictReader(stream)]
    ember, mesa = rows('ember.csv'), rows('mesa.csv')
    eg, mg = rows('ember_grey.csv'), rows('mesa_grey.csv')
    et = np.array([r['Teff_K'] for r in ember])
    el = np.array([r['log_L_Lsun'] for r in ember])
    mt = np.array([r['Teff_K'] for r in mesa])
    ml = np.array([r['log_L_Lsun'] for r in mesa])
    assert all(np.isfinite(v).all() for v in (et, el, mt, ml))
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42})
    pulse = rows('pulse.csv') if (snapshot / 'pulse.csv').exists() else []
    pulse_prefix = rows('pulse_prefix.csv') if (snapshot / 'pulse_prefix.csv').exists() else []
    fig, ax = plt.subplots(figsize=(8.2, 6.6 if pulse else 6.2))
    fig.subplots_adjust(left=.115, right=.97, bottom=.12, top=.81 if pulse else .84)
    orange, blue = '#ad4d28', '#326b91'
    ax.plot(et, el, color=orange, lw=1.9, zorder=4,
            label='Ember: non-grey, element settling')
    ax.plot(mt, ml, color=blue, lw=1.2, alpha=.6, zorder=3,
            label='MESA: atmosphere table, diffusion')
    for seq, color, name in [(eg, orange, 'Ember'), (mg, blue, 'MESA')]:
        tx = np.array([r['Teff_K'] for r in seq])
        ly = np.array([r['log_L_Lsun'] for r in seq])
        assert np.isfinite(tx).all() and np.isfinite(ly).all()
        ax.plot(tx, ly, color=color, ls='--', lw=1.5, alpha=.8,
                label=f'{name}: grey, no diffusion', zorder=2)
    if not pulse:
        ax.plot(et[-1], el[-1], 'o', ms=4.8, color=orange, zorder=6)
    ax.plot(mt[-1], ml[-1], 'o', ms=4.8, color=blue, zorder=6)
    if not pulse:
        ax.annotate(f"Ember: {ember[-1]['age_years']/1e12:.4g} Tyr\n"
                    f"{et[-1]:.4g} K", xy=(et[-1], el[-1]), xytext=(4390, -2.99),
                    ha='left', va='center', color=orange, fontsize=8,
                    arrowprops={'arrowstyle': '-|>', 'color': orange, 'lw': .8,
                                'mutation_scale': 10, 'shrinkA': 4, 'shrinkB': 5})
    ax.annotate(f"MESA: {mesa[-1]['age_years']/1e6:.4g} Myr\n"
                f"{mt[-1]:.4g} K", xy=(mt[-1], ml[-1]), xytext=(3100, -2.23),
                ha='center', va='center', color=blue, fontsize=8,
                arrowprops={'arrowstyle': '-|>', 'color': blue, 'lw': .8,
                            'mutation_scale': 10, 'shrinkA': 4, 'shrinkB': 5})
    # Arrows follow actual saved points; no smoothing, resampling or curve fit.
    for age in (3.60e12, 3.90e12):
        ages = np.array([r['age_years'] for r in ember])
        i = int(np.argmin(abs(ages - age)))
        # Choose a visible segment in the same chronological direction.
        scale = np.hypot((et-et[i])/3600, (el-el[i])/1.3)
        choices = np.where((np.arange(len(ember)) > i) & (scale > .024))[0]
        if len(choices):
            j = int(choices[0])
            ax.annotate('', xy=(et[j], el[j]), xytext=(et[i], el[i]),
                        arrowprops={'arrowstyle': '-|>', 'color': orange,
                                    'lw': 1.3, 'mutation_scale': 11}, zorder=5)
    if eg[-1]['Teff_K'] >= 2500:
        ax.plot(eg[-1]['Teff_K'], eg[-1]['log_L_Lsun'], 'o', ms=4,
                markerfacecolor='white', markeredgecolor=orange, zorder=6)
        ax.annotate(f"Ember grey: {eg[-1]['age_years']/1e12:.4g} Tyr\n"
                    f"{eg[-1]['Teff_K']:.4g} K",
                    xy=(eg[-1]['Teff_K'], eg[-1]['log_L_Lsun']),
                    xytext=(4550, -4.10), ha='left', va='center', color=orange,
                    fontsize=8, arrowprops={'arrowstyle':'-', 'color':orange, 'lw':.7})
    if pulse:
        # The dotted connection indicates ordering across the changed boundary;
        # its interior is not a sequence of computed stellar models.
        first = (pulse_prefix or pulse)[0]
        ax.plot([et[-1], first['Teff_K']], [el[-1], first['log_L_Lsun']],
                color=orange, lw=1.2, ls=':', alpha=.85, zorder=5)
        def connector_point(f):
            return (et[-1] + f * (first['Teff_K'] - et[-1]),
                    el[-1] + f * (first['log_L_Lsun'] - el[-1]))
        ax.annotate('', xy=connector_point(.78), xytext=connector_point(.28),
                    arrowprops={'arrowstyle': '-|>', 'color': orange,
                                'lw': 1.1, 'mutation_scale': 10}, zorder=6)
        for index, seq in enumerate([pulse_prefix, pulse]):
            if not seq:
                continue
            ax.plot([r['Teff_K'] for r in seq], [r['log_L_Lsun'] for r in seq],
                    color=orange, lw=1.7, ls='-', alpha=.85, zorder=5,
                    label='_nolegend_')
        end = pulse[-1]
        ax.plot(end['Teff_K'], end['log_L_Lsun'], 'o', ms=5,
                color=orange, zorder=7)
        ax.annotate(f"Ember: {end['age_years']/1e12:.4g} Tyr\n"
                    f"{end['Teff_K']:.4g} K",
                    xy=(end['Teff_K'], end['log_L_Lsun']), xytext=(4390, -2.71),
                    ha='left', va='center', fontsize=8, color=orange,
                    arrowprops={'arrowstyle': '-|>', 'color': orange, 'lw': .8,
                                'mutation_scale': 10, 'shrinkA': 4, 'shrinkB': 5})
    for seq, color, age in [(eg, orange, 3.61e12), (mg, blue, 2.40e12)]:
        i = min(range(len(seq)), key=lambda j: abs(seq[j]['age_years']-age))
        for j in range(i+1, len(seq)):
            if np.hypot((seq[j]['Teff_K']-seq[i]['Teff_K'])/4700,
                        (seq[j]['log_L_Lsun']-seq[i]['log_L_Lsun'])/2.8) > .025:
                ax.annotate('', xy=(seq[j]['Teff_K'], seq[j]['log_L_Lsun']),
                            xytext=(seq[i]['Teff_K'], seq[i]['log_L_Lsun']),
                            arrowprops={'arrowstyle':'-|>', 'color':color,
                                        'lw':1, 'mutation_scale':10}, zorder=4)
                break
    # Direction on each cooling branch, using chronological saved points
    # after the temperature maximum. The two arrows are offset along the
    # tracks so their arrowheads remain distinct.
    for seq, color, level in [(eg, orange, -2.60), (mg, blue, -2.88)]:
        peak = max(range(len(seq)), key=lambda i: seq[i]['Teff_K'])
        candidates = [i for i in range(peak + 1, len(seq) - 1)
                      if 3200 < seq[i]['Teff_K'] < 6500]
        if not candidates:
            continue
        i = min(candidates, key=lambda j: abs(seq[j]['log_L_Lsun'] - level))
        for j in range(i + 1, len(seq)):
            if np.hypot((seq[j]['Teff_K'] - seq[i]['Teff_K']) / 4700,
                        (seq[j]['log_L_Lsun'] - seq[i]['log_L_Lsun']) / 2.8) > .035:
                ax.annotate('', xy=(seq[j]['Teff_K'], seq[j]['log_L_Lsun']),
                            xytext=(seq[i]['Teff_K'], seq[i]['log_L_Lsun']),
                            arrowprops={'arrowstyle': '-|>', 'color': color,
                                        'lw': 1.2, 'mutation_scale': 11}, zorder=5)
                break
    # Show continued cooling at the frame edge. Arrow endpoints follow the
    # saved cooling sequence; interpolate only its crossing of the plot limit.
    for seq, color in [(eg, orange), (mg, blue)]:
        for i in range(1, len(seq)):
            prev, end = seq[i - 1], seq[i]
            if (prev['Teff_K'] >= 2500 > end['Teff_K']
                    and end['log_L_Lsun'] < -4):
                f = (2500 - prev['Teff_K']) / (end['Teff_K'] - prev['Teff_K'])
                edge_y = prev['log_L_Lsun'] + f * (end['log_L_Lsun'] - prev['log_L_Lsun'])
                candidates = [r for r in seq[:i] if r['log_L_Lsun'] < -4]
                start = min(candidates, key=lambda r: abs(r['Teff_K'] - 2730))
                ax.annotate('', xy=(2500, edge_y),
                            xytext=(start['Teff_K'], start['log_L_Lsun']),
                            arrowprops={'arrowstyle': '-|>', 'color': color,
                                        'lw': 1.3, 'mutation_scale': 12,
                                        'shrinkA': 0, 'shrinkB': 0}, zorder=5)
                break
    ax.set(xlim=(7200, 2500), ylim=(-4.65, -1.85),
           xlabel='Effective surface temperature (K)',
           ylabel=r'$\log_{10}(L/L_\odot)$')
    ax.grid(alpha=.13)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(FuncFormatter(lambda value, _: f'{value:.4g}'))
    handles, labels = ax.get_legend_handles_labels()
    order = [0, 2, 1, 3]
    fig.legend([handles[i] for i in order], [labels[i] for i in order],
               loc='upper left', bbox_to_anchor=(.105, .985), ncol=2,
               columnspacing=1.7, frameon=False, fontsize=8.5)
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output/'full_physics_hr_comparison.pdf', metadata={'CreationDate': None})
    fig.savefig(output/'full_physics_hr_comparison.png', dpi=200)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('snapshot', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    plot(args.snapshot, args.output)
