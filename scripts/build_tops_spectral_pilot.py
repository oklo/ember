#!/usr/bin/env python3
"""Build a small continuous spectral-opacity table and test its interpolation.

Finite-temperature transverse electron dispersion, absorption divided by its
refractive index, unchanged scattering and full Rosseland normalization follow
the individually checked hot-state calculation. A log-bilinear pilot measures
the sampling needed; it is not the runtime stellar interpolation or a complete
composition family. The previous physical opacity tables remain unchanged.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from audit_tops_spectral_means import source, digest
from audit_tops_electron_dispersion import electron_moments, mean, controls, KEV, KB, HBAR, ME


def interpolate(nodes, temperatures, densities, temperature, density):
    ti = min(max(int(np.searchsorted(temperatures, temperature, side='right'))-1, 0), len(temperatures)-2)
    ri = min(max(int(np.searchsorted(densities, density, side='right'))-1, 0), len(densities)-2)
    t0, t1 = temperatures[ti:ti+2]; r0, r1 = densities[ri:ri+2]
    if not t0 <= temperature <= t1 or not r0 <= density <= r1:
        raise ValueError('interpolation would extrapolate')
    ft = math.log(temperature/t0)/math.log(t1/t0)
    fr = math.log(density/r0)/math.log(r1/r0)
    return math.exp(sum(wt*wr*math.log(nodes[t, r])
                        for t, wt in ((t0, 1-ft), (t1, ft))
                        for r, wr in ((r0, 1-fr), (r1, fr))))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('plan', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('use a new report path')
    spec = json.loads(a.plan.read_text())
    paths = [a.plan, Path(__file__)]
    paths += [Path(__file__).with_name(n) for n in (
        'audit_tops_spectral_means.py', 'audit_tops_electron_dispersion.py', 'import_tops_composition.py')]
    inputs = {str(path.resolve()): digest(path) for path in paths}
    limits = controls()
    records = []
    seen = set()
    for directory in spec['reused_sources']+[v['work'] for v in spec['requests']]:
        root = Path(directory)
        inputs.update({str(path.resolve()): digest(path) for path in
                       (root/'receipt.json', root/'request.json', root/'source.txt')})
        spectra, grey = source(root)
        receipt = json.loads((root/'receipt.json').read_text())
        if receipt['X'] != 0 or receipt['Z'] != .02:
            raise ValueError('pilot expects the verified X=0, Z=0.02 mixture')
        composition = (root/'source.txt').read_text().split(
            'No. Fraction Mass Fraction  At. No.  Chem. Sym.  Mat ID.\n')[1].split('Temperature grid')[0]
        helium = [v.split() for v in composition.splitlines() if len(v.split()) == 5 and v.split()[3] == 'He']
        if len(helium) != 1:
            raise ValueError('missing source helium fraction')
        number, mass = map(float, helium[0][:2])
        atomic_mass = 4.002602*number/mass
        for key, data in spectra.items():
            if key in seen:
                raise ValueError('duplicate state in pilot sources')
            seen.add(key)
            t, rho = key
            ne = rho/(atomic_mass*1.66053906660e-24)*grey[key]['free_electrons_per_ion']
            classical_cutoff = HBAR*math.sqrt(4*math.pi*4.80320471257e-10**2*ne/ME)/(t*KEV)
            electron = electron_moments(t*KEV/KB, ne)
            cutoff = classical_cutoff*math.sqrt(electron['plasma_frequency_squared_ratio'])
            values = [float(mean(data, t, cutoff, electron['vstar_squared'], order)) for order in (8, 16)]
            error = abs(values[0]/values[1]-1)
            if error > 1e-7:
                raise ValueError('frequency quadrature did not converge')
            records.append({'temperature_keV': t, 'density_atomic_g_cm3': rho,
                            'rosseland_atomic_cm2_g': values[1],
                            'reported_conditional_rosseland_atomic_cm2_g': grey[key]['rosseland'],
                            'quadrature_relative_change': error, 'electron': electron,
                            'source': str(root), 'free_electron_density_cm3': ne})
    tt, rr = spec['node_temperatures_keV'], spec['node_densities_g_cm3']
    anchors = {(t, r) for t in tt for r in rr}
    if not anchors <= seen:
        raise ValueError('missing interpolation anchor')
    columns = ('rosseland_atomic_cm2_g', 'reported_conditional_rosseland_atomic_cm2_g')
    nodes = {name: {(r['temperature_keV'], r['density_atomic_g_cm3']): r[name]
                    for r in records if (r['temperature_keV'], r['density_atomic_g_cm3']) in anchors}
             for name in columns}
    # Recover the completed reference values without changing its pinned input.
    reference = Path('docs/results/tops_contraction_electron_dispersion_v1.json')
    inputs[str(reference.resolve())] = digest(reference)
    old = json.loads(reference.read_text())
    reference_error = max(abs(nodes[columns[0]][r['temperature_keV'], r['density']]/
                              r['radiative_opacities']['transverse_dispersion']-1)
                          for r in old['records'])
    if reference_error > 1e-12:
        raise ValueError('reused physical means differ from the completed reference')
    comparisons = []
    for row in records:
        key = row['temperature_keV'], row['density_atomic_g_cm3']
        row['role'] = 'anchor' if key in anchors else 'independent_check'
        if key in anchors:
            continue
        comparisons.append({'temperature_keV': key[0], 'density_atomic_g_cm3': key[1],
            'relative_interpolation_errors': {
                name: interpolate(nodes[name], tt, rr, *key)/row[name]-1 for name in columns}})
    if len(comparisons) != 9:
        raise ValueError('pilot needs all nine independent states')
    maximum = {name: max(abs(r['relative_interpolation_errors'][name]) for r in comparisons)
               for name in columns}
    if any(digest(Path(name)) != checksum for name, checksum in inputs.items()):
        raise ValueError('source changed during the calculation')
    report = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
              'limit_checks': limits, 'reused_reference_relative_error': reference_error,
              'records': records, 'comparisons': comparisons,
              'maximum_relative_interpolation_errors': maximum,
              'radiative_interpolation_relative_criterion': spec['radiative_interpolation_relative_criterion'],
              'pilot_interpolation_passed': maximum[columns[0]] <= spec['radiative_interpolation_relative_criterion'],
              'assumptions': ['Atomic mass basis retained; no conversion to runtime baryonic tables yet.',
                  'Collisionless finite-temperature electron model; hot-state collision sensitivity was evaluated separately.',
                  'This interpolation check uses the stated log-bilinear pilot, not Ember runtime derivatives.'],
              'input_sha256': inputs}
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'output': str(a.output), 'states': len(records), 'checks': len(comparisons),
                      'maximum_relative_interpolation_errors': maximum,
                      'pilot_interpolation_passed': report['pilot_interpolation_passed']}), flush=True)


if __name__ == '__main__':
    main()
