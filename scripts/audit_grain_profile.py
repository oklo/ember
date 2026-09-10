#!/usr/bin/env python3
"""Evaluate suspended-grain opacity on an independently checked gas profile.

This fixed-structure experiment does not solve the atmosphere with grains.
Every nonzero condensate must have an explicit optical-material assignment.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from audit_condensates import closure
from generate_nongrey_grid import composition
from grain_opacity import condensate_mass_fractions, prepared, spectrum
from import_nongrey_grid import read_text
from prepare_nongrey_sources import digest, SOURCES


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path); p.add_argument('material', type=Path)
    p.add_argument('manifest', type=Path); p.add_argument('chemistry_report', type=Path)
    p.add_argument('work', type=Path); p.add_argument('--coordinates', nargs=4, type=float, required=True)
    p.add_argument('--radii', nargs='+', type=float, required=True)
    p.add_argument('--report', type=Path)
    a = p.parse_args(); source = prepared(a.source); material = json.loads(a.material.read_text())
    if material['optical_data_file'] != 'Al2O3.dat' or material['shape'] != 'independent compact spheres':
        raise ValueError('this first profile control explicitly covers alumina only')
    manifest = json.loads(a.manifest.read_text()); chemistry_report = json.loads(a.chemistry_report.read_text())
    matches = [m for m in manifest['models'] if [m[k] for k in ['XH', 'X3', 'teff_K', 'log_g']] == a.coordinates]
    records = [r for r in chemistry_report['records'] if r['coordinates'] == a.coordinates and r['mode'] == 'equilibrium']
    if len(matches) != 1 or len(records) != 1:
        raise ValueError('one atmosphere and one equilibrium chemistry calculation required')
    model, record = matches[0], records[0]
    logfile = a.manifest.parent / model['log']
    chemical = Path(chemistry_report['archive']) / record['source']
    if digest(logfile) != model['log_sha256'] or digest(chemical) != record['source_sha256']:
        raise ValueError('atmosphere or chemistry source checksum differs')
    profile = []
    for line in read_text(logfile).rsplit('FINAL MODEL ATMOSPHERE', 1)[1].splitlines():
        words = line.replace('D', 'E').split()
        if len(words) == 11 and words[0].isdigit():
            profile.append(list(map(float, words)))
    profile.reverse(); chemistry = json.loads(read_text(chemical))
    abundance, weights = composition(*a.coordinates[:2], manifest['metals'])
    elements = json.loads((SOURCES / 'synple-elements.json').read_text())
    numbers, masses = {}, {}
    for symbol, number, mass in zip(elements['symbol'], abundance, weights, strict=True):
        if number > 1e-90:
            numbers[symbol.capitalize()] = number
            masses[symbol.capitalize()] = mass * 1.67333e-24 / 1.66053906660e-24
    audit = closure(chemistry, profile, numbers, masses)
    if any(audit[k] != record[k] for k in audit):
        raise ValueError('independent chemistry audit differs from recorded profile')
    fractions = condensate_mass_fractions(chemistry, masses)
    if set(fractions) != {'Al2O3(s,l)'}:
        raise ValueError('additional condensates need optical data: ' + ', '.join(sorted(set(fractions) - {'Al2O3(s,l)'})))
    fraction = fractions['Al2O3(s,l)']
    if abs(fraction.max() / audit['thresholds']['0']['max_condensed_mass_fraction'] - 1) > 1e-12:
        raise ValueError('grain mass normalization differs from independent chemistry audit')
    a.work.mkdir(parents=True, exist_ok=False)
    profile = np.asarray(profile); order = np.argsort(profile[:, 1]); column = profile[order, 1]
    if np.any(np.diff(column) <= 0):
        raise ValueError('invalid atmospheric mass columns')
    grain_column = float(np.trapezoid(fraction[order], column))
    wave = np.geomspace(*material['wavelength_range_um'], material['wavelength_count'])
    results = []
    for index, radius in enumerate(a.radii):
        optical = spectrum(source, material['optical_data_file'], material['density_g_cm3'], [radius], [1.], wave)
        absorption = optical['absorption_cm2_per_g_dust'] * grain_column
        scattering = optical['scattering_cm2_per_g_dust'] * grain_column
        transport = optical['transport_scattering_cm2_per_g_dust'] * grain_column
        output = a.work / f'optical-depth-{index:02d}.csv'
        np.savetxt(output, np.column_stack([wave, absorption, scattering, transport, optical['asymmetry_parameter']]),
            delimiter=',', header='wavelength_um,grain_absorption_tau,grain_scattering_tau,grain_transport_scattering_tau,asymmetry_parameter', comments='', fmt='%.17g')
        j = int(np.argmax(absorption + scattering))
        results.append({'radius_um': radius, 'maximum_extinction_tau': float((absorption + scattering)[j]),
            'wavelength_at_maximum_um': float(wave[j]), 'maximum_absorption_tau': float(absorption.max()),
            'maximum_scattering_tau': float(scattering.max()), 'output': str(output), 'output_sha256': digest(output)})
    report = {'scope': 'All equilibrium condensates retained as suspended alumina on this fixed gas-only structure. No radiative feedback, gas-opacity update, grain growth, or settling calculation. Not an installed atmosphere.',
        'coordinates': a.coordinates, 'material_assumptions': material,
        'maximum_condensed_mass_fraction': float(fraction.max()), 'grain_column_g_cm2': grain_column,
        'chemistry': audit, 'optical_controls': results,
        'integration_scope': 'Between the top and bottom tabulated mass columns; no unmodeled material is added above the top.',
        'wavelength_coverage': 'Optical control restricted to the compiled alumina data. The source atmosphere extends to 0.09 microns, below the 0.2-micron optical-data limit; full transfer coupling must resolve that gap.',
        'input_sha256': {str(path): digest(path) for path in [a.source, a.material, a.manifest, a.chemistry_report, logfile, chemical]},
        'scripts_sha256': {name: digest(Path(__file__).with_name(name)) for name in ['audit_grain_profile.py', 'grain_opacity.py']}}
    (a.work / 'audit.json').write_text(json.dumps(report, indent=2) + '\n')
    if a.report:
        a.report.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'maximum_condensed_mass_fraction': float(fraction.max()), 'grain_column_g_cm2': grain_column,
                      'maximum_extinction_tau': [r['maximum_extinction_tau'] for r in results]}))


if __name__ == '__main__':
    main()
