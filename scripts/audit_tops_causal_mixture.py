"""Apply the causal conductivity diagnostic to a verified TOPS mixture.

The conductivity sum counts bound and free electrons. The collisionless
comparison uses only free electrons; the two electron densities are retained
separately. This is not an accepted dense-opacity table.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess

import numpy as np

from audit_tops_spectral_means import source, digest
from audit_tops_refractive_conventions import means as component_means
from audit_tops_causal_helium import rosseland, C, HBAR, KEV, KB, ME, E
from causal_conductivity import ContinuedConductivity


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('work', 'classical_check', 'conduction_probe', 'mass_probe', 'scratch', 'output'):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    if a.output.exists() or a.scratch.exists():raise FileExistsError('use fresh diagnostic paths')
    a.scratch.mkdir(parents=True)
    check = json.loads(a.classical_check.read_text())
    inputs = dict(check['input_sha256'])
    paths = [a.classical_check, a.conduction_probe, a.mass_probe, Path(__file__),
             Path('scripts/conduction_composition_probe.cpp'),
             Path('/tmp/ember-opacity-mass-basis-v1.cpp'),
             Path('include/ember/composition.hpp'), Path('include/ember/gs98_mixture.hpp')]
    paths += [Path(__file__).with_name(n) for n in ('audit_tops_spectral_means.py',
              'audit_tops_refractive_conventions.py', 'audit_tops_causal_helium.py', 'causal_conductivity.py')]
    inputs.update({str(v.resolve()): digest(v) for v in paths})
    for name, h in inputs.items():
        if digest(Path(name)) != h:raise ValueError('source dependency changed')
    root = a.work/'cutoff-on'
    spectra, source_means = source(root)
    receipt = json.loads((root/'receipt.json').read_text())
    text = (root/'source.txt').read_text().split(
        'No. Fraction Mass Fraction  At. No.  Chem. Sym.  Mat ID.\n')[1].split('Temperature grid')[0]
    composition = [line.split() for line in text.splitlines() if len(line.split()) == 5]
    number_sum = sum(float(r[0]) for r in composition)
    total_charge = sum(float(r[0])*float(r[2]) for r in composition)/number_sum
    if abs(number_sum-1) > 1e-4:raise ValueError('source element numbers are not normalized')
    hscale, hescale, zscale = map(float, subprocess.check_output([str(a.mass_probe)], text=True).split())
    X, Z = receipt['X'], receipt['Z']
    mass_scale = 1/(X/hscale+(1-X-Z)/hescale+Z/zscale)
    xb, zb = mass_scale*X/hscale, mass_scale*Z/zscale
    # This is the inverse of ElementalOpacity's baryonic-to-atomic conversion.
    if max(abs(xb*hscale/mass_scale-X), abs(zb*zscale/mass_scale-Z),
           abs((1-xb-zb)*hescale/mass_scale-(1-X-Z))) > 1e-14:
        raise ValueError('mass-basis round trip failed')
    states = []
    lines = []
    for old in check['records']:
        key = next(k for k in spectra if k[1] == old['density'] and
                   abs(k[0]*KEV/KB/old['temperature_K']-1) < 1e-12)
        states.append((key, old))
        lines.append(f'1 1 {key[1]/mass_scale:.17g} {old["temperature_K"]:.17g} 1 {xb:.17g} 0 {1-xb-zb:.17g}\n')
    request = a.scratch/'conduction.txt';request.write_text(''.join(lines))
    response = subprocess.run([str(a.conduction_probe), 'data/conduction/condtab21wd_metals.dat',
                               str(request), format(zb, '.17g')], text=True, capture_output=True)
    (a.scratch/'conduction.out').write_text(response.stdout)
    (a.scratch/'conduction.err').write_text(response.stderr);response.check_returncode()
    lines = response.stdout.splitlines()
    if len(lines) != len(states)+1 or not lines[-1].startswith('checksum '):raise ValueError('invalid conduction response')
    records = []
    for (key, old), line in zip(states, lines[:-1], strict=True):
        data = spectra[key]
        free_charge = source_means[key]['free_electrons_per_ion']
        bound_fraction = 1-free_charge/total_charge
        if bound_fraction < -1e-5:raise ValueError('free-electron count exceeds all electrons')
        free_ne = old['estimated_electron_density']
        total_ne = free_ne*max(1., total_charge/free_charge)
        scale = key[0]*KEV/HBAR
        total_plasma = math.sqrt(4*math.pi*E**2*total_ne/ME)/scale
        free_plasma = old['classical_cutoff_u']
        D = C*key[1]/scale
        kc_b = float(line.split()[0]);kc_atomic = kc_b/mass_scale
        baseline = component_means(data, key[0], free_plasma, 16)[
            'absorption_divided_by_n_scattering_unchanged']
        baseline_coarse = component_means(data, key[0], free_plasma, 8)[
            'absorption_divided_by_n_scattering_unchanged']
        variants = []
        for fraction in (.5, .75, 1.):
            calculations = []
            for intervals, refinement in ((256, 1), (512, 1), (512, 2)):
                model = ContinuedConductivity(data[:, 0]/key[0], data[:, 2], D,
                                              total_plasma, fraction*free_plasma, refinement)
                result = rosseland(model, data, key[0], intervals)
                # The helper's ancillary total-electron collisionless estimate
                # is not the physical free-electron reference used here.
                result.pop('classical_component_rosseland')
                calculations.append(result)
            final = calculations[-1]
            combined = lambda k: 1/(1/k+1/kc_atomic)
            variants.append({'join_over_free_plasma': fraction, 'collision_u': model.collision_u,
                'sum_relative_error': model.sum_relative_error,
                'rosseland_atomic': final['rosseland'],
                'radiative_relative_to_free_collisionless': final['rosseland']/baseline-1,
                'combined_relative_to_free_collisionless': combined(final['rosseland'])/combined(baseline)-1,
                'frequency_integration_relative_change': calculations[1]['rosseland']/calculations[0]['rosseland']-1,
                'spectral_refinement_relative_change': final['rosseland']/calculations[1]['rosseland']-1,
                'inverse_rosseland_fraction_negative_epsilon1': final['inverse_rosseland_fraction_negative_epsilon1'],
                'calculations': calculations})
        records.append({'temperature_keV': key[0], 'density_atomic': key[1],
                        'density_baryonic': key[1]/mass_scale, 'electron_density_free': free_ne,
                        'electron_density_total': total_ne, 'bound_electron_fraction': max(0., bound_fraction),
                        'conductive_opacity_baryonic': kc_b, 'classical_component_rosseland_atomic': baseline,
                        'classical_quadrature_relative_change': baseline/baseline_coarse-1, 'variants': variants})
        (a.scratch/'partial.json').write_text(json.dumps(records, indent=2)+'\n')
        print(json.dumps({'temperature_keV': key[0], 'density_atomic': key[1], 'variants': variants}), flush=True)
    error = max(abs(v[k]) for r in records for v in r['variants'] for k in
                ('frequency_integration_relative_change', 'spectral_refinement_relative_change'))
    for name, h in inputs.items():
        if digest(Path(name)) != h:raise ValueError('source changed during the calculation')
    report = {'scope': __doc__, 'accepted_for_stellar_opacity': False, 'records': records,
              'mass_scale_atomic_over_baryonic': mass_scale, 'hydrogen_baryonic': xb, 'metals_baryonic': zb,
              'full_mean_charge_per_ion': total_charge, 'maximum_numerical_relative_change': error,
              'numerical_relative_criterion': 1e-4, 'numerical_checks_passed': error <= 1e-4,
              'assumptions': ['The total electronic sum includes bound electrons; only free electrons set the collisionless reference cutoff.',
                 'The measured high-frequency spectrum is retained; the low-frequency part is a continuous Drude model constrained by the total sum.',
                 'Join variation is a model sensitivity, not a global physical error bound. Collective scattering is held fixed.',
                 'Composition and density conversion follows the runtime ElementalOpacity convention; isotope and GS98 proxy approximations remain.'],
              'input_sha256': inputs}
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'report': str(a.output), 'numerical_checks_passed': error <= 1e-4}), flush=True)


if __name__ == '__main__':main()
