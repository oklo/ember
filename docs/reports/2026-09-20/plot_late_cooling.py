"""Plot the checked temperature turn and luminosity decline."""
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


def plot(data, output):
    meta = json.loads((data/'late_cooling_inputs.json').read_text())
    assert hashlib.sha256((data/'late_cooling.csv').read_bytes()).hexdigest() == meta['csv_sha256']
    with (data/'late_cooling.csv').open() as stream:
        rows = [{k: float(v) for k,v in row.items()} for row in csv.DictReader(stream)]
    t = np.array([r['pulse_years']/1e6 for r in rows])
    temperature = np.array([r['Teff_K'] for r in rows])
    light = np.array([r['photon_luminosity_erg_s']/3.828e33 for r in rows])
    nuclear = np.array([r['nuclear_power_erg_s']/3.828e33 for r in rows])
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42})
    fig, (hr, time) = plt.subplots(1, 2, figsize=(9.5, 4.25), layout='constrained')
    orange = '#ad4d28'
    hr.plot(temperature, np.log10(light), color=orange, lw=1.8)
    peak = int(np.argmax(temperature))
    hr.plot(temperature[peak], np.log10(light[peak]), 'o', ms=4, color=orange)
    hr.annotate(f'{t[peak]:.4g} Myr\n{temperature[peak]:.4g} K',
                xy=(temperature[peak], np.log10(light[peak])), xytext=(6376, -1.36),
                fontsize=9, ha='left', va='center',
                arrowprops=dict(arrowstyle='-|>', lw=.8, color='black', mutation_scale=10))
    hr.plot(temperature[-1], np.log10(light[-1]), 'o', ms=4, color=orange)
    hr.annotate(f'{t[-1]:.4g} Myr\n{temperature[-1]:.4g} K',
                xy=(temperature[-1], np.log10(light[-1])), xytext=(6382, -1.58),
                fontsize=9, ha='left', va='center',
                arrowprops=dict(arrowstyle='-|>', lw=.8, color='black', mutation_scale=10))
    i, j = len(rows)*2//3, len(rows)*2//3+4
    hr.annotate('', xy=(temperature[j], np.log10(light[j])),
                xytext=(temperature[i], np.log10(light[i])),
                arrowprops=dict(arrowstyle='-|>', color=orange, lw=1.2, mutation_scale=12))
    hr.set(xlim=(6410, min(6305, temperature[-1]-16)), ylim=(-1.69, -1.28),
           xlabel='Effective surface temperature (K)', ylabel=r'$\log_{10}(L/L_\odot)$')
    time.plot(t, light, color=orange, lw=1.8, label='Surface luminosity')
    time.plot(t, nuclear, color='0.3', lw=1.5, ls='--', label='Nuclear power')
    time.set(xlabel='Time from the model before shell convection\n(Myr)',
             ylabel=r'Luminosity ($L_\odot$)', xlim=(1.02, t[-1]+.03), ylim=(0, .052))
    time.legend(frameon=False, fontsize=9, loc='upper right')
    for ax in (hr, time):
        ax.grid(alpha=.12)
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.4g}'))
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output/'late_cooling.pdf', metadata={'CreationDate': None})
    fig.savefig(output/'late_cooling.png', dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    plot(args.data, args.output)
