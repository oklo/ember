#!/usr/bin/env python3
"""Check grain absorption/scattering against Optool and the Rayleigh limit."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import numpy as np
from grain_opacity import efficiencies, mixture, optical_constants, prepared, spectrum
from prepare_nongrey_sources import digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path); p.add_argument('work', type=Path)
    p.add_argument('--report', type=Path)
    a = p.parse_args(); receipt = prepared(a.source)
    if digest(receipt['optool']) != receipt['optool_sha256']:
        raise ValueError('independent Optool executable changed')
    a.work.mkdir(parents=True, exist_ok=False); work = a.work.resolve()
    # Fixed before comparison: include weak/strong absorption, lossless
    # scattering, and size parameters on both sides of the Rayleigh regime.
    cases = [(1.5, .1, .001), (1.5, .1, .1), (1.5, .1, 2),
             (1.5, .1, 20), (1.5, 0, 2), (3, 4, .1), (3, 4, 1), (3, 4, 10)]
    checks = []
    for index, (n, k, x) in enumerate(cases):
        radius = x / (2 * math.pi); density = 3.
        path = work / f'case-{index:02d}'; path.mkdir()
        command = [receipt['optool'], f'{n}:{k}:{density}', '-mie', '-a', f'{radius:.17g}',
                   '-l', '1', '-s', '1800', '-o', str(path)]
        with (path / 'run.log').open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=30)
        lines = [line for line in (path / 'dustkapscatmat.dat').read_text().splitlines() if line.strip() and not line.startswith('#')]
        if list(map(int, lines[:3])) != [0, 1, 1800]:
            raise ValueError('unexpected independent opacity format')
        wave, absorption, scattering, g = map(float, lines[3].split())
        factor = 4 * density * radius * 1e-4 / 3
        other = np.array([absorption * factor, scattering * factor, (absorption + scattering) * factor, g])
        actual = efficiencies(receipt['probe'], n, k, x)
        # Optool's printed table has seven significant figures. The allowed
        # difference also covers the independent numerical algorithms.
        allowance = 1e-4 * np.abs(other) + 5e-7
        selected = [1, 3] if k == 0 else [0, 1, 2, 3]
        if np.any(np.abs(actual[selected] - other[selected]) > allowance[selected]):
            raise ValueError(f'independent Mie comparison failed: {index}, {actual}, {other}')
        if k == 0 and (actual[0] != 0 or actual[1] != actual[2]):
            raise ValueError('lossless grain must scatter without absorbing')
        checks.append({'n': n, 'k': k, 'size_parameter': x, 'lxmie': actual.tolist(),
                       'optool': other.tolist(), 'absolute_difference': np.abs(actual - other).tolist(),
                       'independently_compared_columns': selected})
    if np.any(efficiencies(receipt['probe'], 1., 0., [1e-5, 1., 1e4])):
        raise ValueError('no-contrast grains have nonzero optical cross section')
    for invalid in ['1.5 0.1\n', '1.5 -0.1 1\n', '1.5 0.1 0\n']:
        r = subprocess.run([receipt['probe']], input=invalid, text=True, capture_output=True, timeout=10)
        if not r.returncode:
            raise ValueError('invalid or incomplete grain query accepted')
    x = .001; m = complex(1.5, .1); polarizability = (m*m - 1) / (m*m + 2)
    expected = np.array([4*x*polarizability.imag, 8/3*x**4*abs(polarizability)**2])
    actual = efficiencies(receipt['probe'], m.real, m.imag, x)
    rayleigh_error = np.abs(actual[:2] / expected - 1)
    if np.max(rayleigh_error) > 1e-4 or abs(actual[3]) > 1e-5:
        raise ValueError('Rayleigh absorption/scattering limit failed')
    wave = np.geomspace(.2, 100, 120)
    first = spectrum(receipt, 'Al2O3.dat', 4., [.01], [1.], wave)
    second = spectrum(receipt, 'Al2O3.dat', 4., [.1], [1.], wave)
    combined = spectrum(receipt, 'Al2O3.dat', 4., [.01, .1], [1., 2.], wave)
    # Independent mass accounting: number weights 1:2 become mass weights
    # 1:2000 for radii differing by ten. Mixing must not use number weights
    # directly to average a per-mass opacity or scattering asymmetry.
    mixed = mixture({'a': first, 'b': second}, {'a': 1/2001, 'b': 2000/2001})
    errors = {}
    for quantity in ['absorption', 'scattering', 'transport_scattering']:
        errors[quantity] = float(np.max(np.abs(combined[quantity + '_cm2_per_g_dust'] / mixed[quantity + '_cm2_per_g_total'] - 1)))
    if max(errors.values()) > 1e-12:
        raise ValueError('grain mass weighting or scattering weighting failed')
    zero = mixture({'a': first}, {'a': 0.})
    if np.any(zero['absorption_cm2_per_g_total']) or np.any(zero['scattering_cm2_per_g_total']):
        raise ValueError('zero retained grain mass does not recover zero grain opacity')
    try:
        optical_constants(receipt, 'Al2O3.dat', [.1])
    except ValueError:
        pass
    else:
        raise ValueError('unsupported optical wavelength was silently extrapolated')
    report = {'scope': 'Spherical particle optics and mass normalization only; no coupled atmosphere or stellar result',
        'source_receipt': str(a.source.resolve()), 'source_receipt_sha256': digest(a.source),
        'probe_sha256': receipt['probe_sha256'], 'optool_sha256': receipt['optool_sha256'],
        'columns': ['Qabs', 'Qsca', 'Qext', 'scattering_asymmetry'], 'independent_cases': checks,
        'acceptance': {'independent_relative': 1e-4, 'independent_absolute': 5e-7, 'rayleigh_relative': 1e-4},
        'independent_angular_resolution': 1800,
        'angular_refinement_reason': 'Optool computes g by angular quadrature. Its default 180-point grid gives a spurious negative g in the small-particle control; resolving that integral retains the original comparison tolerance.',
        'lossless_control': 'Optool imposes Qabs >= 1e-4 times the raw Qext in optool.f90. Compare its Qsca and g only for k=0. Check LX-MIE absorption and extinction against the exact lossless identity separately; no absorption floor is imported.',
        'rayleigh_relative_errors': rayleigh_error.tolist(), 'size_distribution_relative_errors': errors,
        'zero_grain_limit_exact': True, 'rejects_unsupported_wavelength': True,
        'scripts_sha256': {name: digest(Path(__file__).with_name(name)) for name in
            ['grain_mie_probe.cpp', 'grain_opacity.py', 'prepare_grain_sources.py', 'audit_grain_mie.py']},
        'passed': True}
    (work / 'audit.json').write_text(json.dumps(report, indent=2) + '\n')
    if a.report:
        a.report.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'independent_cases': len(checks), 'max_rayleigh_relative_error': float(max(rayleigh_error)), 'passed': True}))


if __name__ == '__main__':
    main()
