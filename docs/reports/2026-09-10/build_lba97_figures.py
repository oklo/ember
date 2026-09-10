#!/usr/bin/env python3
"""Rebuild the LBA97 comparison from recorded figure coordinates and Ember data."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import numpy as np

HERE = Path(__file__).resolve().parent
EMBER = '#ad4d28'
LBA = '#3f7887'


def save(fig, name):
    fig.savefig(HERE/(name+'.pdf'), metadata={'CreationDate': None})
    fig.savefig(HERE/(name+'.png'), dpi=180)
    plt.close(fig)


def main():
    source = json.loads((HERE/'lba97_figure1_digitization.json').read_text())
    with (HERE/'evolution_history.csv').open() as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
    star = {k: np.array([float(r[k]) for r in rows]) for k in
            ['age_yr', 'Teff_K', 'L_Lsun', 'central_X', 'central_Y3']}
    assert np.all(np.diff(star['age_yr']) > 0)
    with (HERE/'evolution_1024_continuation.csv').open() as stream:
        fine_rows = list(csv.DictReader(stream))
    fine = {k: np.array([float(r[k]) for r in fine_rows]) for k in star}
    assert np.all(np.diff(fine['age_yr']) > 0)
    hr = np.asarray(source['pixels']['hr'])
    c = source['calibration']['hr']
    temperature = c['left_Teff_K']+(hr[:, 0]-c['left_px'])*(c['right_Teff_K']-c['left_Teff_K'])/(c['right_px']-c['left_px'])
    loglum = c['top_log10_L_Lsun']+(hr[:, 1]-c['top_px'])*(c['bottom_log10_L_Lsun']-c['top_log10_L_Lsun'])/(c['bottom_px']-c['top_px'])
    with (HERE/'lba97_digitized_hr.csv').open('w', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(['sample', 'pixel_x', 'pixel_y', 'Teff_K', 'log10_L_Lsun'])
        writer.writerows((i, *p, t, l) for i, (p, t, l) in enumerate(zip(hr, temperature, loglum)))
    c = source['calibration']['composition']
    composition = {}
    with (HERE/'lba97_digitized_composition.csv').open('w', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(['species', 'sample', 'pixel_x', 'pixel_y', 'age_yr', 'mass_fraction'])
        for species in ['hydrogen', 'helium3']:
            pixels = np.asarray(source['pixels'][species])
            age = (pixels[:, 0]-c['zero_age_px'])*6e12/(c['age_6T_px']-c['zero_age_px'])
            fraction = (pixels[:, 1]-c['zero_fraction_px'])/(c['unit_fraction_px']-c['zero_fraction_px'])
            assert np.all(np.diff(age) > 0)
            composition[species] = age, fraction
            writer.writerows((species, i, *p, a, f) for i, (p, a, f) in enumerate(zip(pixels, age, fraction)))
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42})
    handles = [Line2D([], [], color=EMBER, lw=1.8),
               Line2D([], [], color=LBA, lw=1.1, alpha=.5),
               Line2D([], [], color=LBA, marker='D', linestyle='none', ms=4),
               Line2D([], [], color='#706078', lw=1.1, linestyle='--', marker='o', markerfacecolor='none', ms=4)]
    labels = ['Ember: 512 mass points', 'LBA97: read from Figure 1', 'LBA97: stated value',
              'Ember: 1024 mass points']
    fig, axes = plt.subplots(1, 2, figsize=(7, 3.65))
    fig.subplots_adjust(left=.105, right=.985, bottom=.17, top=.79, wspace=.34)
    for ax in axes:
        ax.plot(temperature, loglum, color=LBA, lw=1.1, alpha=.5)
        ax.plot(star['Teff_K'], np.log10(star['L_Lsun']), color=EMBER, lw=1.8)
        ax.plot(star['Teff_K'][-1], np.log10(star['L_Lsun'][-1]), 'o', color=EMBER, ms=4)
        ax.plot(fine['Teff_K'], np.log10(fine['L_Lsun']), '--', color='#706078', lw=1.1)
        ax.plot(fine['Teff_K'][-1], np.log10(fine['L_Lsun'][-1]), 'o', color='#706078', mfc='none', ms=5)
        for name in ['main_sequence_start', 'central_radiative_core', 'cooling_endpoint']:
            point = source['published_points'][name]
            ax.plot(point['Teff_K'], point['log10_L_Lsun'], 'D', color=LBA, ms=3.5)
        ax.set(xlabel='Effective surface temperature (K)', ylabel=r'$\log_{10}(L/L_\odot)$')
    axes[0].set(xlim=(6100, 1500), ylim=(-5.4, -2.05), title='Evolution through cooling')
    axes[1].set(xlim=(3550, 2100), ylim=(-3.48, -2.43), title='Hydrogen-burning detail')
    axes[1].annotate(f"Ember at {star['age_yr'][-1]/1e12:.4g} trillion yr",
                     (star['Teff_K'][-1], np.log10(star['L_Lsun'][-1])),
                     xytext=(12, -22), textcoords='offset points', fontsize=8,
                     arrowprops={'arrowstyle': '-', 'color': EMBER, 'lw': .7})
    axes[1].annotate('LBA97: central radiative core', (3450, -2.54),
                     xytext=(18, 12), textcoords='offset points', fontsize=7.5)
    for ax in axes:
        ax.grid(alpha=.15)
        for axis in [ax.xaxis, ax.yaxis]:
            axis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.4g}'))
    fig.legend(handles, labels, loc='upper center', ncol=2, frameon=False, fontsize=8)
    save(fig, 'lba97_hr')

    fig, axes = plt.subplots(1, 2, figsize=(7, 3.2))
    fig.subplots_adjust(left=.10, right=.98, bottom=.18, top=.80, wspace=.34)
    reading = source['reading_tolerance_px']
    abundance_scale = reading/(c['zero_fraction_px']-c['unit_fraction_px'])
    age_scale = reading*6/(c['age_6T_px']-c['zero_age_px'])
    for ax, species, key, ylabel in zip(axes, ['hydrogen', 'helium3'],
                                       ['central_X', 'central_Y3'],
                                       ['Central hydrogen mass fraction', r'Central $^3$He mass fraction']):
        age, fraction = composition[species]
        ax.plot(age/1e12, fraction, color=LBA, lw=1.1, alpha=.5)
        ax.plot(star['age_yr']/1e12, star[key], color=EMBER, lw=1.8)
        ax.plot(star['age_yr'][-1]/1e12, star[key][-1], 'o', color=EMBER, ms=4)
        ax.plot(fine['age_yr']/1e12, fine[key], '--', color='#706078', lw=1.1)
        ax.plot(fine['age_yr'][-1]/1e12, fine[key][-1], 'o', color='#706078', mfc='none', ms=5)
        # A representative reading scale is clearer than a confidence band:
        # errors in tracing a printed curve are correlated and not statistical.
        example = 4 if species == 'hydrogen' else 6
        ax.errorbar(age[example]/1e12, fraction[example], xerr=age_scale,
                    yerr=abundance_scale, fmt='none', color=LBA, alpha=.6, lw=.8, capsize=2)
        ax.set(xlim=(0, 6.35), xlabel='Reported age (trillion yr)', ylabel=ylabel)
        ax.grid(alpha=.15)
        for axis in [ax.xaxis, ax.yaxis]:
            axis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.4g}'))
    core = source['published_points']['central_radiative_core']
    axes[0].plot(core['age_yr']/1e12, core['central_X'], 'D', color=LBA, ms=4)
    peak = source['published_points']['helium3_peak']
    axes[1].plot(peak['age_yr']/1e12, peak['central_Y3'], 'D', color=LBA, ms=4)
    axes[0].set_ylim(-.025, .75)
    axes[1].set_ylim(-.004, .115)
    fig.legend(handles, labels, loc='upper center', ncol=2, frameon=False, fontsize=8)
    save(fig, 'lba97_composition')
    print(f"Ember: {len(rows)} states through {star['age_yr'][-1]/1e12:.4g} trillion yr; "
          f"LBA97: {len(hr)} HR samples, {sum(len(v[0]) for v in composition.values())} abundance samples")


if __name__ == '__main__':
    main()
