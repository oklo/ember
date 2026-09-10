#!/usr/bin/env python3
"""Save a compact reference, convergence comparisons and optional track figure.

Pass completed 512/1024/2048-point runs followed by the tighter 512-point run.
Only --plot requires matplotlib. Full profiles remain in the input files.
"""
import argparse
import hashlib
import json
from pathlib import Path


def record(path):
    raw = path.read_bytes()
    data = json.loads(raw)
    if not data['converged'] or data['history'][-1][0] != 1e12:
        raise ValueError(f'{path}: a completed trillion-year calculation is required')
    history = data['history']
    return data, {
        'input': str(path), 'sha256': hashlib.sha256(raw).hexdigest(),
        'points': data['points'], 'step_error_tolerances': data['step_error_tolerances'],
        'accepted_macrosteps': len(history) - 1, 'rejected_macrosteps': data['rejected_steps'],
        'final': dict(zip(data['columns'], history[-1], strict=True)),
        'peak_He3': {'age_yr': max(history, key=lambda r: r[6])[0],
                     'mass_fraction': max(r[6] for r in history)},
        'max_last_halfstep_luminosity_imbalance': max(abs(r[8]) for r in history[1:]),
        'max_last_halfstep_nuclear_mass_imbalance': max(abs(r[9]) for r in history[1:]),
        'minimum_accepted_convective_mass_fraction': min(r[11] for r in history[1:]),
    }


def difference(first, second):
    a, b = first['final'], second['final']
    return {
        'reference': first['input'], 'comparison': second['input'],
        'relative': {k: b[k] / a[k] - 1 for k in
                     ['R_Rsun', 'L_Lsun', 'Teff_K', 'central_T_K', 'central_rho']},
        'absolute': {k: b[k] - a[k] for k in ['central_X', 'central_Y3']},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs', nargs=4, type=Path)
    parser.add_argument('--output', type=Path, default=Path('docs/results'))
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    loaded = [record(p) for p in args.runs]
    data, records = zip(*loaded, strict=True)
    if [r['points'] for r in records] != [512, 1024, 2048, 512]:
        raise ValueError('expected 512, 1024, 2048, then tighter 512-point runs')
    for run in data:
        if (run['nuclear_model'], run['transport_model'], run['atmosphere_model']) != (
                'sfii-svh', 'wd', 'cond-corrected'):
            raise ValueError('reference summary requires SFII/SVH, wd conduction and corrected COND')
    for run, factor in zip(data, [10, 10, 10, 2.5], strict=True):
        expected = dict(zip(['log_structure', 'absolute_abundance', 'relative_surface_luminosity'],
                            [1e-5 * factor, 1e-8 * factor, 1e-4 * factor], strict=True))
        if run['step_error_tolerances'] != expected:
            raise ValueError('unexpected timestep tolerance scales')
    args.output.mkdir(parents=True, exist_ok=True)
    convergence = {
        'description': 'Same physical prescriptions; numerical comparisons, not physical error bars',
        'runs': records,
        'mesh_comparisons': [difference(records[0], records[1]), difference(records[1], records[2])],
        'timestep_comparison': difference(records[0], records[3]),
    }
    (args.output / 'evolution_1tyr_convergence.json').write_text(
        json.dumps(convergence, indent=2, allow_nan=False) + '\n')
    reference = data[2]
    compact = {k: v for k, v in reference.items() if k not in ['history', 'profile', 'profile_columns', 'columns']}
    compact.update(records[2])
    # Retain every accepted macrostep, but omit step diagnostics from the plot data.
    # Earlier run binaries used an initial convective-fraction placeholder; it is
    # not represented as a measurement here. Accepted fractions are audited above.
    columns = [0, 2, 3, 4, 5, 6, 12, 13]
    compact['history_columns'] = [reference['columns'][i] for i in columns]
    compact['history'] = [[row[i] for i in columns] for row in reference['history']]
    compact['atmosphere_integration_used'] = {
        'log_state_tolerance': 2e-8, 'sensitivity_tolerance': 2e-8,
        'tau_top': .001, 'alpha': 1.9, 'henyey_y': 1 / 3,
    }
    compact['reference_run_provenance_note'] = (
        'Long references began with the Z=.01/.02 TOPS source pair; the later Z=.03 plane '
        'does not enter their mapped-Z interpolation. Subsequent driver changes add '
        'atmosphere settings and a measured initial convective fraction to JSON. '
        'The original initial-fraction placeholder is omitted from this compact history.')
    (args.output / 'evolution_m010_1tyr.json').write_text(
        json.dumps(compact, indent=2, allow_nan=False) + '\n')
    if args.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
        fig, axes = plt.subplots(2, 2, figsize=(9, 6), sharex=True, constrained_layout=True)
        history = reference['history']
        age = [r[0] / 1e12 for r in history]
        for i, label, color in [(5, 'Hydrogen-1', '#b34b24'), (6, 'Helium-3', '#246e91')]:
            axes[0, 0].plot(age, [r[i] for r in history], label=label, color=color)
        axes[0, 0].plot(age, [.98 - r[5] - r[6] for r in history],
                        label='Helium-4', color='#657344')
        axes[0, 0].set_ylabel('Baryonic mass fraction')
        axes[0, 0].legend(frameon=False, loc='center left', bbox_to_anchor=(.02, .58))
        for ax, i, label, scale in [
            (axes[0, 1], 3, r'Luminosity ($10^{-3}\,L_\odot$)', 1000),
            (axes[1, 0], 2, r'Radius ($R_\odot$)', 1),
            (axes[1, 1], 4, 'Effective temperature (K)', 1),
        ]:
            ax.plot(age, [scale * r[i] for r in history], color='#246e91')
            ax.set_ylabel(label)
        for ax in axes.flat:
            ax.grid(alpha=.18)
            ax.set_xlim(0, 1)
        for ax in axes[1]:
            ax.set_xlabel('Elapsed time (trillion years)')
        fig.suptitle('0.1 solar mass: first trillion years of the specified main-sequence model\n'
                     '2048 points · SFII/SVH · COND with convective composition correction', fontsize=12)
        fig.savefig(args.output / 'evolution_m010_1tyr.pdf')
        fig.savefig(args.output / 'evolution_m010_1tyr.png', dpi=180)


if __name__ == '__main__':
    main()
