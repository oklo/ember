#!/usr/bin/env python3
"""Compare the TOPS cutoff inference with its reported free-electron density.

The mean atomic mass is estimated from the returned helium number and mass
fractions and the NIST standard helium atomic weight. This is a classical,
collisionless diagnostic, not an adopted dielectric model. The conduction
comparison uses Ember's existing baryonic GS98 mixture approximation at the
same X and Z; it does not construct a stellar structure.
"""
import argparse
import json
import math
from pathlib import Path
import re
import subprocess

from audit_tops_spectral_means import Integrals, source, digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['work', 'spectral_report', 'conduction_probe', 'output']:
        p.add_argument(name, type=Path)
    a = p.parse_args()
    table = Path('data/conduction/condtab21wd_metals.dat')
    root = a.work/'cutoff-on'
    inputs = {str(path.resolve()): digest(path) for path in
              [root/'receipt.json', root/'request.json', root/'source.txt',
               a.spectral_report, a.conduction_probe, table, Path(__file__),
               Path(__file__).with_name('audit_tops_spectral_means.py'),
               Path('scripts/benchmark_conduction_profile.cpp')]}
    spectra, means = source(root)
    receipt = json.loads((root/'receipt.json').read_text())
    report = json.loads(a.spectral_report.read_text())
    for name, checksum in report['input_sha256'].items():
        if digest(Path(name)) != checksum:
            raise ValueError('spectral comparison input changed')
        inputs[name] = checksum
    parts = (root/'source.txt').read_text().split(
        'No. Fraction Mass Fraction  At. No.  Chem. Sym.  Mat ID.\n')[1].split('Temperature grid')[0]
    helium = [row.split() for row in parts.splitlines() if re.search(r'\sHe\s', row)]
    if len(helium) != 1 or receipt['Z'] != .02:
        raise ValueError('expected one helium entry and the GS98 Z=.02 control')
    number, mass = map(float, helium[0][:2])
    atomic_weight = 4.002602
    mean_atomic_mass = atomic_weight*number/mass
    rows = []
    records = []
    for state in report['states']:
        temperature, rho = state['temperature_keV'], state['density']
        key = temperature, rho
        T = temperature*1.602176634e-9/1.380649e-16
        ne = rho/(mean_atomic_mass*1.66053906660e-24)*means[key]['free_electrons_per_ion']
        u = 1.054571817e-27*math.sqrt(4*math.pi*4.80320471257e-10**2*ne/9.1093837015e-28)/(temperature*1.602176634e-9)
        integral = Integrals(spectra[key], temperature, 16)
        above = integral.above(u)
        step = (4*math.pi**4/15)/above[1]
        refractive = (4*math.pi**4/15)/integral.refractive(u)[1]
        x = receipt['X']
        rows.append(f'1 1 {rho:.17g} {T:.17g} 1 {x:.17g} 0 {1-.02-x:.17g}\n')
        inferred = state['calculations'][-1]['inferred_cutoff_u']
        records.append(dict(temperature_K=T, density=rho, hydrogen=x, helium3=0.,
                            estimated_electron_density=ne, classical_cutoff_u=u,
                            inferred_cutoff_u=inferred, inferred_cutoff_relative_difference=inferred/u-1,
                            reported_conditional_rosseland=means[key]['rosseland'],
                            independently_predicted_conditional_rosseland=above[0]/above[1],
                            full_normalization_classical_step_rosseland=step,
                            full_normalization_classical_refractive_rosseland=refractive))
    request = a.work/(a.output.stem+'.conduction.txt')
    request.write_text(''.join(rows))
    r = subprocess.run([str(a.conduction_probe.resolve()), str(table.resolve()),
                        str(request.resolve()), '1'], text=True, capture_output=True, check=True)
    lines = r.stdout.splitlines()
    if len(lines) != len(records)+1 or not lines[-1].startswith('checksum '):
        raise ValueError('unexpected conduction response')
    for record, line in zip(records, lines[:-1], strict=True):
        values = list(map(float, line.split()))
        if len(values) != 5 or not all(map(math.isfinite, values)) or values[0] <= 0:
            raise ValueError('invalid conductive opacity')
        kc = values[0]
        original = record['reported_conditional_rosseland']
        combined = 1/(1/original+1/kc)
        record.update(conductive_opacity=kc, original_radiative_fraction=combined/original,
                      original_combined_opacity=combined,
                      maximum_combined_increase_if_radiation_vanishes=kc/combined-1)
        for name in ['step', 'refractive']:
            kr = record[f'full_normalization_classical_{name}_rosseland']
            total = 1/(1/kr+1/kc)
            record[f'{name}_combined_opacity'] = total
            record[f'{name}_combined_relative_change'] = total/combined-1
    inputs[str(request.resolve())] = digest(request)
    result = dict(scope=__doc__, accepted_for_stellar_opacity=False,
                  mean_atomic_mass_from_helium=mean_atomic_mass,
                  helium_atomic_weight=atomic_weight,
                  atomic_weight_reference='https://physics.nist.gov/PhysRefData/Handbook/Tables/heliumtable1_a.htm',
                  assumptions='Classical plasma frequency from reported free electrons; standard helium mass and rounded source fractions. No relativistic, collisional or bound-electron dispersion. Conductive mixture is the existing Ember approximation.',
                  input_sha256=inputs, records=records)
    if any(digest(Path(name)) != checksum for name, checksum in inputs.items()):
        raise ValueError('diagnostic input changed')
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(records, indent=2))


if __name__ == '__main__':
    main()
