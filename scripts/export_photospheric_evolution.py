#!/usr/bin/env python3
"""Recover photospheric states from the accepted atmosphere source structures.

This is a plotting diagnostic, not a change to the evolution boundary at tau=100.
The atmosphere source's density is retained at Rosseland optical depth 2/3.
"""
import argparse
from bisect import bisect_right
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess

import numpy as np

from generate_nongrey_grid import temperatures, sequence
from import_nongrey_grid import read_text, source_state
from nongrey_opacity import read_table


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def weights(axes, values, available):
    """Use the same complete-cell choice as CompositionAtmosphereGrid."""
    if any(not a[0] <= q <= a[-1] for a, q in zip(axes, values)):
        raise ValueError('photospheric query outside atmosphere axes')
    preferred = [min(bisect_right(a, q)-1, len(a)-2) for a, q in zip(axes, values)]
    for alternative in range(16):
        base = preferred.copy()
        for k in range(4):
            if alternative & (1 << k):
                if not base[k] or values[k] != axes[k][base[k]]:
                    break
                base[k] -= 1
        else:
            corners = list(itertools.product(*[(j, j+1) for j in base]))
            if not all(c in available for c in corners):
                continue
            u = [(q-a[j])/(a[j+1]-a[j]) for a, q, j in zip(axes, values, base)]
            return [(c, math.prod(u[k] if c[k] != base[k] else 1-u[k]
                                  for k in range(4))) for c in corners]
    raise ValueError('missing complete source cell for photospheric query')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('paper', type=Path)
    p.add_argument('atmosphere_manifest', type=Path)
    p.add_argument('solar_opacity', type=Path)
    p.add_argument('--probe', type=Path, required=True)
    a = p.parse_args()
    history_record = json.loads((a.paper/'evolution_latest_provenance.json').read_text())
    history_path = a.paper/'evolution_latest.csv'
    assert sha(history_path) == history_record['history_csv_sha256']
    history = [{k: float(v) if v else None for k, v in r.items()}
               for r in csv.DictReader(history_path.open())]
    manifest = json.loads(a.atmosphere_manifest.read_text())
    table_path = Path(history_record['atmosphere'].split(':', 1)[1])
    assert sha(table_path) == manifest['table_sha256']
    source_hashes = {str(a.atmosphere_manifest): sha(a.atmosphere_manifest)}
    base_path = Path(manifest['base_manifest'])
    assert sha(base_path) == manifest['input_files_sha256'][str(base_path)]
    base = json.loads(base_path.read_text())
    records = [(base_path.parent, r, base, r['diagnostics']) for r in base['models']]
    for entry in manifest['extension']:
        work = Path(entry['work'])
        record_path = work/'final/validated.json'
        assert sha(record_path) == entry['validation_sha256']
        record = json.loads(record_path.read_text())
        spec = entry['source_specification']
        spec = {**spec, 'opacity': {'temperature_K': temperatures(spec),
                                   'density_g_cm3': sequence(spec['log_density'])}}
        records.append((work, record, spec, entry['state']))
        source_hashes[str(record_path)] = sha(record_path)
    assert len(records) == manifest['accepted_states']
    axes = []
    for line in table_path.read_text().splitlines():
        if line.split()[0] in ['hydrogen', 'helium3', 'log_teff', 'log_g']:
            axes.append([float(v) for v in line.split()[2:]])
    states, nodes = {}, []
    for root, record, spec, expected in records:
        text = {}
        for kind in ['log', 'convergence']:
            path = root/record[kind]
            assert sha(path) == record[kind+'_sha256']
            source_hashes[str(path)] = sha(path)
            text[kind] = read_text(path)
        common = (text['log'], text['convergence'], record['teff_K'], record['log_g'], spec['opacity'])
        deep = source_state(*common, tau=100)
        photo = source_state(*common, tau=2/3)
        for key in ['T', 'Pgas', 'source_density']:
            assert abs(deep[key]/expected[key]-1) < 1e-12
        q = [record['XH'], record['X3'], math.log10(record['teff_K']), record['log_g']]
        index = tuple(min(range(len(ax)), key=lambda j: abs(ax[j]-v)) for ax, v in zip(axes, q))
        assert all(abs(ax[j]-v) < 1e-13 for ax, j, v in zip(axes, index, q))
        assert index not in states
        states[index] = [math.log(photo[k]) for k in ['T', 'source_density', 'Pgas']]
        states[index] += [math.log(deep[k]) for k in ['T', 'Pgas']]
        nodes.append({'coordinates': q, 'photosphere': photo, 'deep_matching_state': deep})
    output, check_inputs, check_expected = [], [], []
    samples = set(np.linspace(0, len(history)-1, 65).round().astype(int))
    for i, row in enumerate(history):
        logg = math.log10(.1*1.3271244e26/(row['R_Rsun']*6.957e10)**2)
        q = [row['surface_X'], row['surface_Y3'], math.log10(row['Teff_K']), logg]
        stencil = weights(axes, q, states)
        fields = [math.exp(sum(w*states[c][k] for c, w in stencil)) for k in range(5)]
        output.append([row['age_yr'], row['R_Rsun'], row['Teff_K'], logg,
                       row['surface_X'], row['surface_Y3'], *fields[:3]])
        if i in samples:
            check_inputs.append([q[0], q[1], row['Teff_K'], logg])
            check_expected.append(fields[3:])
    result = subprocess.run([str(a.probe), history_record['EOS'].split(':', 1)[1], str(table_path)],
                            input=''.join(' '.join(format(v, '.17g') for v in r)+'\n' for r in check_inputs),
                            capture_output=True, text=True, check=True)
    actual = [json.loads(line) for line in result.stdout.splitlines()]
    assert len(actual) == len(check_expected) and all(r['covered'] for r in actual)
    runtime_error = max(abs(r[k]/expected[j]-1) for r, expected in zip(actual, check_expected)
                        for j, k in enumerate(['T', 'Pgas']))
    assert runtime_error < 1e-11
    columns = ['age_yr', 'R_Rsun', 'Teff_K', 'log_g', 'surface_X', 'surface_Y3',
               'photosphere_T_K', 'photosphere_rho_g_cm3', 'photosphere_Pgas']
    csv_path = a.paper/'photospheric_evolution.csv'
    with csv_path.open('w', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n'); writer.writerow(columns); writer.writerows(output)

    # A fixed hydrogen-rich reference map, as in LBA97 Figure 6. These are
    # absorption-only means over the actual sampled wavelength interval.
    plane = next(r for r in base['opacity_planes'] if r['XH'] == .7 and r['X3'] == 0)
    assert sha(a.solar_opacity) == plane['sha256']
    table = read_table(a.solar_opacity)
    nu = np.asarray(table['frequency'])[::-1]
    logk = np.asarray(table['log_opacity'], dtype=float).reshape(table['shape'])[::-1]
    means, weight_fractions = [], []
    for j, lt in enumerate(table['log_temperature']):
        temperature = math.exp(lt)
        x = 6.6256e-27*nu/(1.38054e-16*temperature)
        # Normalized dB_nu/dT weight; constants independent of nu cancel.
        w = x**4*np.exp(-x)/(-np.expm1(-x))**2
        norm = np.trapezoid(w, x)
        inverse = np.trapezoid(w[:, None]*np.exp(-logk[:, :, j]), x, axis=0)
        means.append(norm/inverse)
        weight_fractions.append(float(norm/(4*math.pi**4/15)))
    means = np.asarray(means)
    assert np.all(np.isfinite(means)) and np.all(means > 0)
    background = {'hydrogen': .7, 'helium3': 0., 'Z': .02, 'metal_pattern': 'GS98',
                  'kind': 'Rosseland absorption mean over the source frequency interval; scattering and grains excluded',
                  'source_path': str(a.solar_opacity), 'source_sha256': sha(a.solar_opacity),
                  'log10_temperature_K': (np.asarray(table['log_temperature'])/math.log(10)).tolist(),
                  'log10_density_g_cm3': (np.asarray(table['log_density'])/math.log(10)).tolist(),
                  'log10_kappa_cm2_g': np.log10(means).tolist(),
                  'rosseland_weight_fraction_in_frequency_interval': weight_fractions,
                  'integration': 'Trapezoid in h nu / kT using source constants; harmonic absorption mean, normalized over sampled interval.'}
    map_path = a.paper/'photospheric_opacity_map.json'
    map_path.write_text(json.dumps(background, indent=2)+'\n')
    report = {'reference': 'https://www.astroexplorer.org/details/10_1086_304125_fg6',
              'history_sha256': sha(history_path), 'photosphere_csv_sha256': sha(csv_path),
              'opacity_map_sha256': sha(map_path), 'script_sha256': sha(__file__),
              'atmosphere_table_sha256': sha(table_path), 'source_files_sha256': source_hashes,
              'source_nodes_checked': len(nodes), 'history_states': len(output),
              'optical_depth': 2/3, 'density': 'Baryonic density from the TLUSTY source molecular EOS.',
              'interpolation': 'Log T, density and gas pressure in XH, XHe3, log Teff and log g, using complete source cells. No photospheric EOS inversion using the deeper stellar EOS.',
              'runtime_matching_comparison': {'samples': len(actual), 'maximum_relative_error': runtime_error,
                                              'probe_sha256': sha(a.probe)},
              'scope': 'Reconstruction from accepted atmosphere structures; no independent interpolation-accuracy bound at tau=2/3. Fixed initial-composition absorption background is a reference, not the changing total opacity along the track. No early contraction or future cooling added.',
              'first_state': dict(zip(columns, output[0])), 'last_state': dict(zip(columns, output[-1]))}
    (a.paper/'photospheric_evolution_provenance.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ['source_nodes_checked', 'history_states', 'runtime_matching_comparison', 'first_state', 'last_state']}, indent=2))


if __name__ == '__main__':
    main()
