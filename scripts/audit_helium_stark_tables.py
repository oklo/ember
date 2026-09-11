#!/usr/bin/env python3
"""Inspect published helium-line tables without changing stellar opacity.

The finite wavelength ranges are integrated as supplied. A separate 4471-A
wing estimate uses the asymptotic |wavelength offset|**(-5/2) law discussed
by Tremblay et al. (2026); it is a diagnostic, not a renormalization or an
accepted opacity prescription. The 2009 profiles have no Doppler broadening.
"""
import argparse
from collections import defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import tarfile

import numpy as np
from scipy.interpolate import PchipInterpolator


TEMPERATURES = np.array([10000., 20000., 40000.])
DENSITIES = np.r_[10.**np.arange(13., 18., .5), 6.e17]
SIMULATED_LINES = {3820, 3868, 3965, 4026, 4121, 4144, 4169,
                   4388, 4438, 4471, 4713, 4922, 5048}
PUBLISHED_MD5 = {
    'Tremblay26': '701bba709891d0cf84e16850614876e2',
    'Beauchamp25_NLD': '7ece842049d868f3ef2d7ca03f956781',
    'Beauchamp25_LD': '0e7d05c8c57552fa44fb6db292253f73',
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_2026(path, name):
    data = Path(path).read_bytes()
    if hashlib.md5(data).hexdigest() != PUBLISHED_MD5[name]:
        raise ValueError('table differs from the published Zenodo v1 checksum')
    tokens = data.decode().split()
    result, i = {}, 0
    while i < len(tokens):
        n, label, wavelength = tokens[i:i+3]
        n, wavelength = int(n), int(wavelength)
        if label != 'LAMBDA' or wavelength in result or n < 3:
            raise ValueError('invalid or duplicate line header')
        i += 3
        x = np.array(tokens[i:i+n], float)
        i += n
        y = np.array(tokens[i:i+33*n], float).reshape(11, 3, n)
        i += 33*n
        if (not np.all(np.isfinite(x)) or not np.all(np.diff(x) > 0)
                or not np.all(np.isfinite(y)) or not np.all(y >= 0)):
            raise ValueError('invalid wavelength grid or profile')
        result[wavelength] = (x, y)
    if len(result) != 36:
        raise ValueError('expected all 36 published transitions')
    return result


def read_2009(path, original_tar):
    data = gzip.decompress(Path(path).read_bytes())
    rows = data.decode().splitlines()
    if len(rows) != 900:
        raise ValueError('expected 900 CDS rows')
    groups = defaultdict(list)
    for row in rows:
        values = [float(row[k:k+8]) if row[k:k+8].strip() else np.nan
                  for k in range(21, 154, 12)]
        groups[float(row[:5])].append([float(row[8:17]), *values])
    if len(groups) != 10:
        raise ValueError('expected ten electron densities')
    # The original authors tables provide a separate format check on the
    # fixed-width CDS export. Read archive members without extracting paths.
    with tarfile.open(original_tar) as archive:
        for i, (ne, rows) in enumerate(sorted(groups.items()), 2):
            text = archive.extractfile(f'table{i:02d}.txt').read().decode()
            numeric = [line for line in text.splitlines()
                       if line.strip() and not line.lstrip().startswith('#')]
            cds = np.array(rows)
            if len(numeric) != len(cds):
                raise ValueError('original and CDS wavelength counts differ')
            for line, reference in zip(numeric, cds):
                values = np.array([float(v) for v in line.split()])
                if not np.array_equal(values, reference[np.isfinite(reference)]):
                    raise ValueError('original and CDS profile values differ')
    return {ne: np.array(rows) for ne, rows in groups.items()}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in PUBLISHED_MD5:
        p.add_argument('--'+name.lower().replace('_', '-'), type=Path, required=True)
    p.add_argument('--gigosos-cds', type=Path, required=True)
    p.add_argument('--gigosos-original', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    paths = {name: getattr(a, name.lower()) for name in PUBLISHED_MD5}
    inputs = {str(path.resolve()): digest(path) for path in
              [*paths.values(), a.gigosos_cds, a.gigosos_original, Path(__file__)]}
    tables = {name: read_2026(path, name) for name, path in paths.items()}
    unchanged = 0
    for line, (x, y) in tables['Tremblay26'].items():
        xb, yb = tables['Beauchamp25_NLD'][line]
        if not np.array_equal(x, xb):
            raise ValueError('published comparison wavelength grids differ')
        limit = 5 if line == 4713 else 2 if line in SIMULATED_LINES else 11
        if not np.array_equal(y[:limit], yb[:limit]):
            raise ValueError('published nonsimulation entries differ from NLD')
        unchanged += limit*3
    areas = {}
    for name, lines in tables.items():
        summaries = []
        for line, (x, y) in sorted(lines.items()):
            linear = np.trapezoid(y, x, axis=-1)
            cubic = PchipInterpolator(x, y, axis=-1).integrate(x[0], x[-1])
            summary = {'wavelength_A': line, 'samples': len(x),
                       'offset_range_A': [x[0], x[-1]],
                       'linear_area_range': [linear.min(), linear.max()],
                       'pchip_area_range': [cubic.min(), cubic.max()]}
            if line == 4471:
                tail = (abs(x[0])*y[..., 0]+x[-1]*y[..., -1])/1.5
                summary.update(linear_areas=linear.tolist(), pchip_areas=cubic.tolist(),
                               asymptotic_wing_estimates=tail.tolist(),
                               pchip_plus_wing_estimate=(cubic+tail).tolist())
            summaries.append(summary)
        areas[name] = summaries
    groups = read_2009(a.gigosos_cds, a.gigosos_original)
    ion_mass = []
    for log_ne_si, data in sorted(groups.items()):
        x, profiles = data[:, 0], data[:, 1:].reshape(-1, 3, 4)
        for it, temperature in enumerate([5000, 10000, 20000, 40000]):
            y = profiles[:, :, it]
            expected_missing = (temperature == 5000 and log_ne_si > 22.33
                                or temperature == 10000 and log_ne_si > 23.33)
            if expected_missing:
                if not np.all(np.isnan(y)):
                    raise ValueError('unexpected data in a published missing cell')
                continue
            if not np.all(np.isfinite(y)) or not np.all(y >= 0) or not np.all(np.diff(x) > 0):
                raise ValueError('invalid 2009 source profile')
            ion_mass.append({'log10_electron_density_cm3': log_ne_si-6.,
                             'temperature_K': temperature,
                             'reduced_masses_proton_units': [.8, 2., 4.],
                             'finite_range_areas': np.trapezoid(y, x, axis=0).tolist(),
                             'hydrogen_vs_helium_ion_L1': float(np.trapezoid(abs(y[:, 0]-y[:, 1]), x))})
    if any(digest(path) != checksum for path, checksum in inputs.items()):
        raise ValueError('an input changed during the audit')
    report = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
              'sources': ['https://doi.org/10.5281/zenodo.18722143',
                          'https://arxiv.org/abs/2603.04374',
                          'https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A%2BA/503/293'],
              'input_sha256': inputs, 'profiles_per_2026_file': 1188,
              'simulation_profiles': 1188-unchanged,
              'nonsimulation_profiles_exactly_equal_to_NLD': unchanged,
              'temperature_K': TEMPERATURES.tolist(), 'electron_density_cm3': DENSITIES.tolist(),
              'area_diagnostics': areas, 'ion_mass_diagnostics_2009': ion_mass,
              'limitations': ['No opacity or source executable is modified.',
                  'The 2026 simulated perturbers are He II ions; proton dynamics require comparison.',
                  'The simulated profiles omit line dissolution.',
                  'Finite-range areas are not unit-area normalization tests.',
                  'The 4471-A wing estimate assumes the asymptotic power law holds at each boundary.',
                  'No wavelength, temperature or density extrapolation is accepted by this audit.',
                  'The 2009 ion-mass differences also contain simulation sampling noise.']}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'simulation_profiles': 1188-unchanged, 'unchanged_profiles': unchanged,
                      'ion_mass_comparisons': len(ion_mass), 'output': str(a.output)}))


if __name__ == '__main__':
    main()
