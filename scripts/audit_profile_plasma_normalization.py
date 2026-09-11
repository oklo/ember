#!/usr/bin/env python3
"""Bound the plasma-cutoff normalization effect on one fixed stellar profile.

Complete ionization of H and He plus Z/A<=1/2 for GS98 metals bounds the
electron density. The classical electron plasma frequency then supplies an
upper cutoff for this diagnostic. Opacity spectra, conductivity, composition
and stellar structure are otherwise held fixed. This does not bound the full
refractive-index correction or the change in a newly evolved stellar model.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess

from scipy.integrate import quad
from run_evolution_snapshot import input_data


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('profile', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--probe', type=Path, required=True)
    p.add_argument('--eos-family', type=Path, required=True)
    p.add_argument('--opacity-directory', type=Path, required=True)
    a = p.parse_args()
    rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(a.profile.open())]
    columns = ['mass_g', 'radius_cm', 'density_g_cm3', 'temperature_K',
               'luminosity_erg_s', 'X', 'Y3', 'Y4']
    request = ''.join(' '.join(format(row[k], '.17g') for k in columns)+'\n' for row in rows)
    command = [str(a.probe), '--gs98', '--eos-family', str(a.eos_family),
               '--opacity-directory', str(a.opacity_directory)]
    material = input_data(['metal:'+str(a.eos_family), '--opacity-directory', str(a.opacity_directory)])
    inputs = {str(path.resolve()): digest(path) for path in [a.profile, a.probe, Path(__file__), *material]}
    result = subprocess.run(command, input=request, text=True, capture_output=True, check=True)
    physics = [json.loads(line) for line in result.stdout.splitlines()]
    if len(physics) != len(rows):
        raise ValueError('profile diagnostics changed row count')
    records = []
    for i, (row, local) in enumerate(zip(rows, physics, strict=True)):
        T, rho = row['temperature_K'], row['density_g_cm3']
        if abs(sum(row[k] for k in ['X', 'Y3', 'Y4'])-.98) > 1e-12:
            raise ValueError('profile is not the fixed Z=.02 mixture')
        ye = row['X']+2*row['Y3']/3+row['Y4']/2+.02/2
        ne = rho*ye/1.66053906660e-24
        omega = math.sqrt(4*math.pi*(4.80320471257e-10)**2*ne/9.1093837015e-28)
        cutoff = 1.054571817e-27*omega/(1.380649e-16*T)
        def weight(u):
            if u == 0:
                return 0.
            return u**4*math.exp(-u)/math.expm1(-u)**2
        removed, error = quad(weight, 0, cutoff, epsabs=1e-12, epsrel=1e-12)
        removed /= 4*math.pi**4/15
        radiative_fraction = local[6]/local[4]
        records.append(dict(zone=i, mass_fraction=row['mass_g']/rows[-1]['mass_g'],
                            temperature_K=T, density=rho, upper_free_electrons_per_baryon=ye,
                            upper_classical_cutoff_u=cutoff, upper_removed_rosseland_weight=removed,
                            original_radiative_opacity=local[4], original_combined_opacity=local[6],
                            radiative_fraction_of_diffusive_conductance=radiative_fraction,
                            upper_radiative_opacity_relative_increase=1/(1-removed)-1,
                            upper_combined_opacity_relative_increase=1/(1-radiative_fraction*removed)-1))
    keys = ['upper_classical_cutoff_u', 'upper_removed_rosseland_weight',
            'upper_radiative_opacity_relative_increase', 'upper_combined_opacity_relative_increase']
    report = dict(scope=__doc__, input_sha256=inputs, command=command,
                  probe_response_sha256=hashlib.sha256(result.stdout.encode()).hexdigest(),
                  maxima={key: max(records, key=lambda r: r[key]) for key in keys},
                  central=records[0], zones=records, accepted_for_stellar_opacity=False,
                  limits='Continuous classical cutoff; the finite TOPS frequency-bin rounding is not included. Cooler zones use a complete-ionization upper bound, not an assumed ionization state. No claim about evolutionary response follows from a fixed-profile conductance change.')
    if any(digest(Path(path)) != value for path, value in inputs.items()):
        raise ValueError('an input changed during the diagnostic')
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({key: report['maxima'][key][key] for key in keys}))


if __name__ == '__main__':
    main()
