#!/usr/bin/env python3
"""Measure the heat-transport effect of failed radiative-ratio comparisons.

Keep the original radiative failures. Evaluate the existing conduction table
with the same atomic/baryonic mass mapping, both with the selected temperature
activation and, separately, without that activation. The latter is a physical
sensitivity calculation, not an adopted change to the stellar prescription.
"""
import argparse
from decimal import Decimal, localcontext
import json
import math
from pathlib import Path
import subprocess

from audit_tops_electron_dispersion import KEV, KB
from fetch_tops_composition import digest
from reduce_tops_group_factors import add_inputs, verify


def combined_change(delta, radiative, conduction):
    return delta / (1 + (1 + delta) * radiative / conduction)


def controls():
    maximum = 0.
    count = 0
    with localcontext() as context:
        context.prec = 90
        for rad in (1., 1e30):
            for cond in (1e-10, 1., 1e20):
                for delta in (-.07, .1, 1e-8):
                    r, c, d = [Decimal(str(x)) for x in (rad, cond, delta)]
                    old = 1 / (1 / r + 1 / c)
                    new = 1 / (1 / (r * (1 + d)) + 1 / c)
                    exact = float(new / old - 1)
                    value = combined_change(delta, rad, cond)
                    maximum = max(maximum, abs(value / exact - 1))
                    count += 1
    if maximum > 1e-12:
        raise ValueError('harmonic-sum comparison failed')
    return {'cases': count, 'maximum_relative_error': maximum}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('comparison', 'reduction', 'conduction_reference', 'conduction_probe',
                 'mass_probe', 'scratch', 'output'):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    if a.output.exists() or a.scratch.exists():
        raise FileExistsError('preserve completed sensitivity outputs')
    comparison, reduction, conduction_ref = [json.loads(path.read_text()) for path in
                                            (a.comparison, a.reduction, a.conduction_reference)]
    inputs = {}
    for report in (comparison, reduction, conduction_ref):
        add_inputs(inputs, report['input_sha256'])
        add_inputs(inputs, report.get('output_sha256', {}))
    for path in (a.comparison, a.reduction, a.conduction_reference, a.mass_probe,
                 Path('/tmp/ember-opacity-mass-basis-v1.cpp'), Path(__file__),
                 Path('src/conduction_table.cpp')):
        add_inputs(inputs, {str(path.resolve()): digest(path)})
    if inputs.get(str(a.conduction_probe.resolve())) != digest(a.conduction_probe):
        raise ValueError('unidentified conduction executable')
    verify(inputs)
    control = controls()
    criterion = comparison['relative_criterion']
    selected = [r for r in comparison['comparisons'] if abs(r['relative_error']) > criterion]
    if not selected or any((r['X'], r['Z']) != (reduction['X'], reduction['Z']) for r in selected):
        raise ValueError('expected failed comparisons at the reference composition')
    records = {(r['temperature_keV'], r['density_atomic_g_cm3']): r for r in reduction['records']}
    for row in selected:
        direct = records[row['temperature_keV'], row['density_atomic_g_cm3']]
        if (direct['source'] != row['reference_source'] or direct['refractive_ratio'] != row['direct_ratio']):
            raise ValueError('different direct reference source')
    hscale, hescale, zscale = map(float, subprocess.check_output([str(a.mass_probe)], text=True).split())
    x, z = reduction['X'], reduction['Z']
    scale = 1 / (x / hscale + (1 - x - z) / hescale + z / zscale)
    xb, zb = scale * x / hscale, scale * z / zscale
    if max(abs(xb * hscale / scale - x), abs(zb * zscale / scale - z)) > 1e-14:
        raise ValueError('mass-basis round trip failed')
    a.scratch.mkdir(parents=True)
    request = a.scratch / 'conduction.txt'
    request.write_text(''.join(
        f"1 1 {r['density_atomic_g_cm3']/scale:.17g} {r['temperature_keV']*KEV/KB:.17g} 1 {xb:.17g} 0 {1-xb-zb:.17g}\n"
        for r in selected))
    command = [str(a.conduction_probe), 'data/conduction/condtab21wd_metals.dat', str(request), format(zb, '.17g')]
    response = subprocess.run(command, text=True, capture_output=True, timeout=60)
    (a.scratch / 'conduction.stdout').write_text(response.stdout)
    (a.scratch / 'conduction.stderr').write_text(response.stderr)
    response.check_returncode()
    lines = response.stdout.splitlines()
    if len(lines) != len(selected) + 1 or not lines[-1].startswith('checksum '):
        raise ValueError('incomplete conduction response')
    rows = []
    for row, line in zip(selected, lines[:-1], strict=True):
        raw = float(line.split()[0]) / scale
        if not math.isfinite(raw) or raw <= 0:
            raise ValueError('invalid conductive opacity')
        direct = records[row['temperature_keV'], row['density_atomic_g_cm3']]
        rad = direct['native_uncut_times_ratio_atomic_cm2_g']
        temperature = row['temperature_keV'] * KEV / KB
        u = min(1., max(0., math.log(temperature / 3e5) / math.log(1e6 / 3e5)))
        w = u*u*(3-2*u)
        delta = row['relative_error']
        selected_change = delta if w == 0 else combined_change(delta, rad, raw / w)
        full_change = combined_change(delta, rad, raw)
        rows.append({**row, 'temperature_K': temperature,
                     'direct_radiative_opacity_atomic_cm2_g': rad,
                     'conductive_opacity_without_activation_atomic_cm2_g': raw,
                     'selected_conduction_activation': w,
                     'selected_combined_relative_error': selected_change,
                     'full_conduction_combined_relative_error': full_change})
    verify(inputs)
    result = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
              'reference_X_atomic': x, 'reference_Z_atomic': z,
              'source_comparisons': comparison['source_comparisons'],
              'unchanged_radiative_criterion': criterion,
              'original_radiative_comparison_passed': comparison['sampled_comparison_passed'],
              'radiative_failures_evaluated': len(rows),
              'radiative_passes_need_no_conduction_assumption': comparison['source_comparisons'] - len(rows),
              'selected_combined_failure_count': sum(abs(r['selected_combined_relative_error']) > criterion for r in rows),
              'full_conduction_combined_failure_count': sum(abs(r['full_conduction_combined_relative_error']) > criterion for r in rows),
              'maximum_selected_combined_error_among_radiative_failures': max(abs(r['selected_combined_relative_error']) for r in rows),
              'maximum_full_conduction_combined_error_among_radiative_failures': max(abs(r['full_conduction_combined_relative_error']) for r in rows),
              'mass_scale_atomic_over_baryonic': scale, 'harmonic_sum_controls': control,
              'conduction_command': command, 'input_sha256': inputs, 'records': rows,
              'limitations': ['Already-passing radiative errors bound their combined-opacity errors for any nonnegative conduction.',
                              'The selected temperature-only activation suppresses conduction below 300000 K, even in dense material.',
                              'Removing that activation here does not validate complete ionization or all low-temperature conduction physics.',
                              'Neither weighting nor source normalization converts the original radiative failure into a pass.',
                              'This is a material-state sensitivity, not a new stellar trajectory or family acceptance.'],
              'scratch_sha256': {p.name: digest(p) for p in a.scratch.iterdir() if p.is_file()}}
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in
                     ('scope', 'input_sha256', 'records', 'limitations', 'scratch_sha256')}), flush=True)


if __name__ == '__main__':
    main()
