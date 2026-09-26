#!/usr/bin/env python3
"""Locate cold opacity interpolation errors and test conduction sensitivity.

Retain all independent radiative comparisons. The unactivated conduction
table describes fully ionized material; evaluating it here is a sensitivity,
not an accepted prescription for the partly ionized envelope. Source electron
counts quantify that limitation. Saved stellar profiles provide context, not
new evolved models or an opacity acceptance domain.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess

import numpy as np

from audit_tops_electron_dispersion import KEV, KB, ME, C, HBAR
from audit_tops_factor_combined_transport import combined_change, controls
from fetch_tops_composition import digest
from reduce_tops_group_factors import add_inputs, verify

MU = 1.66053906660e-24


def ion_inventory(source):
    text = (Path(source) / 'source.txt').read_text()
    block = text.split('No. Fraction Mass Fraction  At. No.  Chem. Sym.  Mat ID.\n')[1].split('Temperature grid')[0]
    rows = [line.split() for line in block.splitlines() if len(line.split()) == 5]
    number_sum = sum(float(row[0]) for row in rows)
    mass_sum = sum(float(row[1]) for row in rows)
    if max(abs(number_sum - 1), abs(mass_sum - 1)) > 5e-5:
        raise ValueError('source printed composition does not close')
    helium = [row for row in rows if row[3] == 'He']
    if len(helium) != 1:
        raise ValueError('expected one helium component')
    # Use the same source convention as tops_groups.source_groups, without
    # renormalizing its printed number fractions or electron counts.
    mean_mass = 4.002602 * float(helium[0][0]) / float(helium[0][1])
    full_charge = sum(float(row[0]) * int(row[2]) for row in rows)
    return mean_mass, full_charge


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scratch', required=True, type=Path)
    parser.add_argument('--report', required=True, type=Path)
    args = parser.parse_args()
    if args.scratch.exists() or args.report.exists():
        raise FileExistsError('preserve completed transport-regime outputs')
    inputs = {}

    def pin(path):
        path = Path(path)
        add_inputs(inputs, {str(path.resolve()): digest(path)})
        return path

    def report(path):
        result = json.loads(pin(path).read_text())
        add_inputs(inputs, result['input_sha256'])
        add_inputs(inputs, result.get('output_sha256', {}))
        return result

    comparison = report('docs/results/tops_cool_interpolation_components_v2.json')
    report('docs/results/tops_ratio_family_hydrogen_x010_transport_v1.json')
    for name in ('audit_tops_cool_transport_regime.py', 'audit_tops_factor_combined_transport.py',
                 'audit_tops_electron_dispersion.py', 'reduce_tops_group_factors.py'):
        pin(Path(__file__).with_name(name))
    probe = Path('/tmp/ember-conduction-composition-probe-v1')
    mass_probe = Path('/tmp/ember-opacity-mass-basis-v1')
    for path in (probe, mass_probe, Path('data/conduction/condtab21wd_metals.dat')):
        if inputs.get(str(path.resolve())) != digest(path):
            raise ValueError('unidentified conduction or mass-basis input')
    refs = {}
    for label in ('010', '030', '070'):
        reduction = report(f'docs/results/tops_cool_ratio_x{label}-z020-control_reduction_v1.json')
        refs.update({(reduction['X'], reduction['Z'], r['temperature_keV'], r['density_atomic_g_cm3']): r
                     for r in reduction['records'] if r['temperature_keV'] < .025})
    if len(refs) != 420 or len(comparison['comparisons']) != 420:
        raise ValueError('expected all 420 cold reference comparisons')
    early = json.loads(Path('/tmp/ember-hermite-value-benchmark-v1/run-1.json').read_text())
    profiles = {
        'early_10_Gyr': [dict(zip(early['profile_columns'], row, strict=True)) for row in early['profile']],
        'evolved_3890_Gyr': [{k: float(v) for k, v in row.items()} for row in
                           csv.DictReader(Path('docs/reports/2026-09-11/evolution_latest_profile.csv').open())]}
    verify(inputs)
    harmonic_controls = controls()
    hscale, hescale, zscale = map(float, subprocess.check_output([str(mass_probe)], text=True).split())
    args.scratch.mkdir()
    rows, commands = [], []
    inventory_cache = {}
    for x in (.1, .3, .7):
        z = .02
        scale = 1 / (x / hscale + (1 - x - z) / hescale + z / zscale)
        xb, zb = scale * x / hscale, scale * z / zscale
        selected = [r for r in comparison['comparisons'] if r['X'] == x]
        # The only omitted grid density is 1e-10: it is outside the conductive
        # source domain. Do not extrapolate it. Its radiative errors bound the
        # combined-opacity errors for any nonnegative conduction.
        supported = [r for r in selected if r['density_atomic_g_cm3'] >= 1e-5]
        omitted = [r for r in selected if r not in supported]
        if any(r['density_atomic_g_cm3'] != 1e-10 for r in omitted):
            raise ValueError('unexpected omitted conductive source coordinate')
        request = args.scratch / f'x{x:g}-conduction.txt'
        request.write_text(''.join(
            f"1 1 {r['density_atomic_g_cm3']/scale:.17g} {r['temperature_keV']*KEV/KB:.17g} 1 {xb:.17g} 0 {1-xb-zb:.17g}\n"
            for r in supported))
        command = [str(probe), 'data/conduction/condtab21wd_metals.dat', str(request), format(zb, '.17g')]
        commands.append(command)
        response = subprocess.run(command, text=True, capture_output=True, timeout=60)
        (args.scratch / f'x{x:g}-conduction.stdout').write_text(response.stdout)
        (args.scratch / f'x{x:g}-conduction.stderr').write_text(response.stderr)
        response.check_returncode()
        lines = response.stdout.splitlines()
        if len(lines) != len(supported) + 1 or not lines[-1].startswith('checksum '):
            raise ValueError('incomplete conduction query')
        conduction = {(r['temperature_keV'], r['density_atomic_g_cm3']): float(line.split()[0]) / scale
                      for r, line in zip(supported, lines[:-1], strict=True)}
        for row in selected:
            t, rho = row['temperature_keV'], row['density_atomic_g_cm3']
            ref = refs[x, z, t, rho]
            if ref['source'] != row['reference_source']:
                raise ValueError('source identity differs')
            source = ref['source']
            if source not in inventory_cache:
                inventory_cache[source] = ion_inventory(source)
            mean_mass, full_charge = inventory_cache[source]
            full_ne = rho / (MU * mean_mass) * full_charge
            ne = ref['free_electron_density_cm3']
            fraction = ne / full_ne
            if not 0 < fraction <= 1.00005:
                raise ValueError('invalid source free-electron fraction')
            temperature = t * KEV / KB
            if temperature >= 3e5:
                raise ValueError('expected no selected conduction below its activation temperature')
            q = HBAR * (3 * math.pi**2 * ne)**(1/3) / (ME * C)
            tf = ME * C * C / KB * q * q / (math.sqrt(1 + q*q) + 1)
            rad = ref['native_uncut_times_ratio_atomic_cm2_g']
            cond = conduction.get((t, rho))
            if cond is not None and not math.isfinite(cond) or cond is not None and cond <= 0:
                raise ValueError('invalid conduction response')
            context = {}
            for label, profile in profiles.items():
                ordered = sorted(profile, key=lambda r: r['temperature_K'])
                tt = np.log([r['temperature_K'] for r in ordered])
                if not tt[0] <= math.log(temperature) <= tt[-1]:
                    raise ValueError('profile interpolation would extrapolate')
                density = float(np.exp(np.interp(math.log(temperature), tt,
                                   np.log([r['density_g_cm3'] for r in ordered]))))
                context[label] = {'density_baryonic_g_cm3': density,
                                  'reference_over_profile_density': rho / scale / density}
            rows.append({**row, 'temperature_K': temperature, 'mass_scale_atomic_over_baryonic': scale,
                         'source_free_electron_density_cm3': ne,
                         'fully_ionized_electron_density_from_printed_source_cm3': full_ne,
                         'source_free_over_fully_ionized_electrons': fraction,
                         'temperature_over_free_electron_fermi_temperature': temperature / tf,
                         'cutoff_energy_over_kT': ref['cutoff_u'],
                         'direct_radiative_opacity_atomic_cm2_g': rad,
                         'bare_conductive_opacity_sensitivity_atomic_cm2_g': cond,
                         'radiative_flux_fraction_if_bare_conduction_applied': None if cond is None else cond / (rad + cond),
                         'selected_combined_relative_error': row['diagnostic_product_relative_error'],
                         'bare_conduction_combined_relative_error_sensitivity': None if cond is None else combined_change(
                             row['diagnostic_product_relative_error'], rad, cond),
                         'bare_conduction_omission_reason': 'outside source density coverage' if cond is None else None,
                         'saved_profile_context': context})
    failures = [r for r in rows if abs(r['diagnostic_product_relative_error']) > .005]
    if any(r['bare_conductive_opacity_sensitivity_atomic_cm2_g'] is None for r in failures):
        raise ValueError('an original failure lacks a conduction sensitivity')
    summary = {'reference_states': len(rows), 'bare_conduction_queries': sum(
        r['bare_conductive_opacity_sensitivity_atomic_cm2_g'] is not None for r in rows),
        'selected_combined_failures': len(failures),
        'bare_conduction_sensitivity_failures': sum(abs(r['bare_conduction_combined_relative_error_sensitivity']) > .005
                                                   for r in failures),
        'maximum_bare_conduction_sensitivity_error_among_original_failures': max(
            abs(r['bare_conduction_combined_relative_error_sensitivity']) for r in failures),
        'free_electron_fraction_range_among_original_failures': [min(r['source_free_over_fully_ionized_electrons'] for r in failures),
                                                                max(r['source_free_over_fully_ionized_electrons'] for r in failures)],
        'reference_over_early_profile_density_range_among_original_failures': [
            min(r['saved_profile_context']['early_10_Gyr']['reference_over_profile_density'] for r in failures),
            max(r['saved_profile_context']['early_10_Gyr']['reference_over_profile_density'] for r in failures)]}
    verify(inputs)
    result = {'scope': __doc__, 'analysis_complete': True, 'accepted_for_stellar_opacity': False,
              'accepted_change_to_conduction': False, 'radiative_relative_criterion': .005,
              'summary': summary, 'harmonic_sum_controls': harmonic_controls, 'commands': commands,
              'records': rows, 'input_sha256': inputs,
              'limitations': [
                  'The fully ionized conductive calculation is not valid automatically when source electrons remain bound.',
                  'At every reference temperature the selected conduction prescription is off; the radiative failures remain failures.',
                  'The sampled source electron fractions have the precision of the printed source composition.',
                  'The source free-electron Fermi temperature diagnoses degeneracy but does not validate conductivity.',
                  'Saved profiles have different compositions and ages; density comparisons are context, not acceptance criteria.',
                  'The early profile is at 10 Gyr, not the exact initial structure.',
                  'The radiative diagnostic product is not a resampled stellar opacity table.'],
              'output_sha256': {str(p.resolve()): digest(p) for p in args.scratch.iterdir() if p.is_file()}}
    args.report.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
