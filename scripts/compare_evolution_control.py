#!/usr/bin/env python3
"""Compare completed stellar models at exactly the same age and mass mesh.

This records a fresh-start/source-family control, not a convergence error
bound. No interpolation in age, relabelling of a restart, or extrapolation
is performed. Different intermediate stopping ages can change step phasing.
"""
import argparse
import json
import math
from pathlib import Path

from prepare_nongrey_sources import digest


def compare(old, new):
    for data in [old, new]:
        if data.get('converged') is not True or not data.get('history') or not data.get('profile'):
            raise ValueError('both stellar calculations must have completed')
        if len(set(data['columns'])) != len(data['columns']):
            raise ValueError('ambiguous history columns')
    required = ['mass_basis', 'mass_Msun', 'age_origin', 'points',
                'metal_inventory', 'convection_criterion', 'transport_model', 'nuclear_model']
    if any(old[k] != new[k] for k in required):
        raise ValueError('control changes mass, clock, mesh count or non-input-family physics')
    if old.get('thermal_neutrino_model', 'none') != new.get('thermal_neutrino_model', 'none'):
        raise ValueError('control changes the thermal-neutrino prescription')
    histories = [dict(zip(d['columns'], d['history'][-1], strict=True)) for d in [old, new]]
    if histories[0]['age_yr'] != histories[1]['age_yr']:
        raise ValueError('control must stop at exactly the same age')
    metrics = {}
    for name in ['R_Rsun', 'L_Lsun', 'Teff_K', 'central_X', 'central_Y3',
                 'central_T_K', 'central_rho', 'surface_X', 'surface_Y3', 'convective_mass_fraction']:
        x, y = (h[name] for h in histories)
        if not all(math.isfinite(v) for v in [x, y]):
            raise ValueError('nonfinite final history')
        metrics[name] = {'reference': x, 'control': y, 'difference': y-x,
                         'relative_difference': (y-x)/x if x else None}
    if old['profile_columns'] != new['profile_columns'] or len(old['profile']) != len(new['profile']):
        raise ValueError('control must retain the same profile definition and mass mesh')
    columns = old['profile_columns']; profiles = {}
    mass_index = columns.index('mass_g')
    for a, b in zip(old['profile'], new['profile'], strict=True):
        if len(a) != len(columns) or len(b) != len(columns) or a[mass_index] != b[mass_index]:
            raise ValueError('control changes the mass mesh')
        if not all(math.isfinite(v) for v in a+b):
            raise ValueError('nonfinite stellar profile')
    for i, name in enumerate(columns):
        differences = [abs(b[i]-a[i]) for a, b in zip(old['profile'], new['profile'], strict=True)]
        point = max(range(len(differences)), key=differences.__getitem__)
        fractions = [(abs((b[i]-a[i])/a[i]), j) for j, (a, b) in enumerate(zip(old['profile'], new['profile'], strict=True)) if a[i]]
        relative, location = max(fractions) if fractions else (None, None)
        profiles[name] = {'max_absolute_difference': differences[point], 'absolute_max_point': point,
                          'max_relative_difference_where_reference_nonzero': relative,
                          'relative_max_point': location}
    return {'age_yr': histories[0]['age_yr'], 'unchanged_assumptions': {k: old[k] for k in required},
            'final_history': metrics, 'profile': profiles,
            'input_descriptions': [{k: d.get(k) for k in ['eos_model', 'atmosphere_model', 'opacity', 'opacity_directory']}
                                   for d in [old, new]],
            'step_error_tolerances': [d['step_error_tolerances'] for d in [old, new]],
            'energy_diagnostics': [{k: v for k, v in h.items() if 'balance' in k or 'luminosity' in k}
                                   for h in histories]}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['reference', 'control', 'output']:
        p.add_argument(name, type=Path)
    a = p.parse_args()
    inputs = [a.reference, a.control]; hashes = {str(f.resolve()): digest(f) for f in inputs}
    result = compare(*(json.loads(f.read_text()) for f in inputs))
    if any(digest(path) != value for path, value in hashes.items()):
        raise ValueError('stellar models changed during comparison')
    report = {'scope': __doc__, 'script_sha256': digest(__file__), 'input_files_sha256': hashes, **result}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(report['final_history'], indent=2))


if __name__ == '__main__':
    main()
