"""Wavelength-dependent absorption and scattering by independent spheres.

Optical efficiencies come from the separately prepared LX-MIE executable.
Cross sections are integrated over explicit discrete NUMBER weights, then
divided by the mass of those grains. No abundance, size, material density,
settling efficiency, or out-of-range optical constant is supplied implicitly.
"""
import json
import math
from pathlib import Path
import subprocess
import numpy as np
from prepare_nongrey_sources import digest


def prepared(path):
    result = json.loads(Path(path).read_text())
    if result.get('format') != 1 or digest(result['probe']) != result['probe_sha256']:
        raise ValueError('grain source executable identity differs')
    return result


def optical_constants(receipt, material, wavelength_um):
    if material not in receipt['optical_data']:
        raise ValueError('material absent from pinned optical compilation')
    path = Path(receipt['sources']['lxmie']['source']) / 'compilation' / material
    if digest(path) != receipt['optical_data'][material]['sha256']:
        raise ValueError('optical constants checksum differs')
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] != 3 or len(data) < 2 or not np.isfinite(data).all():
        raise ValueError('invalid optical constants')
    if np.any(data[:, :2] <= 0) or np.any(data[:, 2] < 0) or np.any(np.diff(data[:, 0]) <= 0):
        raise ValueError('invalid optical constants or wavelength ordering')
    wave = np.asarray(wavelength_um, dtype=float)
    if wave.ndim != 1 or not len(wave) or not np.isfinite(wave).all() or np.any(wave <= 0):
        raise ValueError('invalid requested wavelengths')
    if wave.min() < data[0, 0] or wave.max() > data[-1, 0]:
        raise ValueError(f'{material}: requested wavelengths outside published compilation')
    # Linear n and logarithmic positive k against log wavelength. Intervals
    # touching an exact zero use linear k, preserving zero endpoints.
    axis = np.log(data[:, 0]); target = np.log(wave)
    j = np.clip(np.searchsorted(axis, target, side='right') - 1, 0, len(axis) - 2)
    f = (target - axis[j]) / (axis[j + 1] - axis[j])
    n = (1 - f) * data[j, 1] + f * data[j + 1, 1]
    low, high = data[j, 2], data[j + 1, 2]
    k = (1 - f) * low + f * high
    positive = (low > 0) & (high > 0)
    k[positive] = np.exp((1 - f[positive]) * np.log(low[positive]) + f[positive] * np.log(high[positive]))
    return n, k


def efficiencies(probe, n, k, x):
    n, k, x = np.broadcast_arrays(np.asarray(n, float), np.asarray(k, float), np.asarray(x, float))
    if not n.size or not all(np.isfinite(v).all() for v in [n, k, x]) or np.any(n <= 0) or np.any(k < 0) or np.any(x < 1e-5) or np.any(x > 1e4):
        raise ValueError('unsupported Mie query')
    queries = ''.join(f'{a:.17g} {b:.17g} {c:.17g}\n' for a, b, c in zip(n.flat, k.flat, x.flat, strict=True))
    run = subprocess.run([str(probe)], input=queries, text=True, capture_output=True, timeout=120)
    if run.returncode:
        raise ValueError('Mie source failed: ' + run.stderr.strip())
    rows = np.array([list(map(float, line.split())) for line in run.stdout.splitlines()])
    if rows.shape != (n.size, 4) or not np.isfinite(rows).all() or np.any(rows[:, :3] < 0) or np.any(np.abs(rows[:, 3]) > 1):
        raise ValueError('invalid Mie source response')
    if np.any(np.abs(rows[:, 0] + rows[:, 1] - rows[:, 2]) > 1e-12 * np.maximum(rows[:, 2], 1e-300)):
        raise ValueError('Mie extinction does not equal absorption plus scattering')
    return rows.reshape((*n.shape, 4))


def spectrum(receipt, material, density_g_cm3, radii_um, number_weights, wavelength_um):
    radii, weights, wave = (np.asarray(v, dtype=float) for v in [radii_um, number_weights, wavelength_um])
    if radii.ndim != 1 or not len(radii) or weights.shape != radii.shape or not np.isfinite(radii).all() or not np.isfinite(weights).all() or np.any(radii <= 0) or np.any(weights < 0) or not weights.sum() > 0:
        raise ValueError('invalid explicit grain size distribution')
    if not math.isfinite(density_g_cm3) or density_g_cm3 <= 0:
        raise ValueError('invalid grain material density')
    n, k = optical_constants(receipt, material, wave)
    q = efficiencies(receipt['probe'], n[None, :], k[None, :], 2 * math.pi * radii[:, None] / wave[None, :])
    weights = weights / weights.sum(); radius_cm = radii * 1e-4
    mass = np.sum(weights * (4 * math.pi / 3) * density_g_cm3 * radius_cm**3)
    area = weights[:, None] * math.pi * radius_cm[:, None]**2
    absorption = np.sum(area * q[:, :, 0], axis=0) / mass
    scattering = np.sum(area * q[:, :, 1], axis=0) / mass
    scattering_g = np.sum(area * q[:, :, 1] * q[:, :, 3], axis=0) / mass
    g = np.divide(scattering_g, scattering, out=np.zeros_like(scattering), where=scattering > 0)
    return {'wavelength_um': wave, 'absorption_cm2_per_g_dust': absorption,
            'scattering_cm2_per_g_dust': scattering, 'asymmetry_parameter': g,
            'transport_scattering_cm2_per_g_dust': scattering - scattering_g}


def mixture(spectra, dust_mass_fractions):
    """Opacity per TOTAL mass, including gas and dust; no gas opacity added.

    The caller supplies a spectrum for every retained condensate. A transport
    coefficient is reported for diagnosis, not substituted into a radiative
    transfer calculation that requires the full scattering phase function.
    """
    if not spectra or set(spectra) != set(dust_mass_fractions):
        raise ValueError('every retained condensate requires optical data')
    fractions = list(dust_mass_fractions.values())
    if any(not math.isfinite(v) or v < 0 for v in fractions) or sum(fractions) > 1:
        raise ValueError('invalid condensate mass fractions')
    wave = next(iter(spectra.values()))['wavelength_um']
    absorption = np.zeros_like(wave); scattering = np.zeros_like(wave); scattering_g = np.zeros_like(wave)
    for material, table in spectra.items():
        if not np.array_equal(table['wavelength_um'], wave):
            raise ValueError('grain wavelength grids differ')
        fraction = dust_mass_fractions[material]
        absorption += fraction * table['absorption_cm2_per_g_dust']
        scattering += fraction * table['scattering_cm2_per_g_dust']
        scattering_g += fraction * table['scattering_cm2_per_g_dust'] * table['asymmetry_parameter']
    g = np.divide(scattering_g, scattering, out=np.zeros_like(scattering), where=scattering > 0)
    return {'wavelength_um': wave, 'absorption_cm2_per_g_total': absorption,
            'scattering_cm2_per_g_total': scattering, 'asymmetry_parameter': g,
            'transport_scattering_cm2_per_g_total': scattering - scattering_g}


def condensate_mass_fractions(chemistry, element_masses):
    """Convert equilibrium formula-unit densities into fractions of total mass.

    Element masses must use the atmosphere's abundance convention. An
    independent chemistry conservation audit is required before this call.
    Every nonzero condensate is retained, including species with no optical
    data; the caller must diagnose missing material coverage explicitly.
    """
    if chemistry['status'] or chemistry['mode'] != 'equilibrium':
        raise ValueError('converged equilibrium chemistry required')
    if any(r['flag'] or not all(r['element_conserved']) for r in chemistry['rows']):
        raise ValueError('unconverged or unconserved chemistry row')
    elements = chemistry['elements']
    weights = np.array([0. if e == 'e-' else element_masses[e] for e in elements])
    if not np.isfinite(weights).all() or np.any(weights < 0) or any(weights[i] <= 0 for i, e in enumerate(elements) if e != 'e-'):
        raise ValueError('invalid elemental mass convention')
    gas = np.asarray([r['gas'] for r in chemistry['rows']], float)
    solid = np.asarray([r['condensed'] for r in chemistry['rows']], float)
    if not all(np.isfinite(v).all() and np.all(v >= 0) for v in [gas, solid]):
        raise ValueError('invalid chemical number density')
    gas_mass = np.asarray([s['stoichiometry'] for s in chemistry['gas_species']]) @ weights
    solid_mass = np.asarray([s['stoichiometry'] for s in chemistry['condensates']]) @ weights
    condensed = solid * solid_mass[None, :]
    total = gas @ gas_mass + condensed.sum(axis=1)
    if not np.isfinite(total).all() or np.any(total <= 0):
        raise ValueError('invalid total material density')
    return {s['symbol']: condensed[:, j] / total for j, s in enumerate(chemistry['condensates'])
            if np.any(condensed[:, j] > 0)}
