#!/usr/bin/env python3
"""Compare frequency groups and full spectra at specified stellar compositions.

This checks the group representation, not table interpolation or complete
dielectric physics. The source mixture, coordinates and electron counts must
agree before either transport result is compared.
"""
import argparse
import json
import math
from pathlib import Path

from audit_tops_spectral_means import digest, source
from audit_tops_electron_dispersion import KEV, KB, HBAR, ME, electron_moments, mean
from audit_tops_group_transport import constant_control
from tops_groups import source_groups, transport


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('plan', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('keep completed comparisons immutable')
    spec = json.loads(a.plan.read_text())
    inputs = dict(spec['seed_input_sha256'])
    code = ['audit_tops_group_compositions.py', 'audit_tops_group_transport.py',
            'tops_groups.py', 'audit_tops_spectral_means.py',
            'audit_tops_electron_dispersion.py', 'import_tops_composition.py',
            'fetch_tops_group_plan.py', 'fetch_tops_spectral_plan.py']
    for path in [a.plan, *(Path(__file__).with_name(n) for n in code)]:
        inputs[str(path.resolve())] = digest(path)

    def pin(root):
        for name in ['receipt.json', 'request.json', 'source.txt', 'recipe.json']:
            path = root/name
            inputs[str(path.resolve())] = digest(path)

    def same_mixture(root, job):
        receipt = json.loads((root/'receipt.json').read_text())
        if any(receipt[k] != job[k] for k in ('X', 'Z', 'metals')):
            raise ValueError('comparison mixture mismatch')

    for name, h in inputs.items():
        if digest(Path(name)) != h:
            raise ValueError('changed control input: '+name)
    rows = []
    control = constant_control()
    for job in spec['mixtures']:
        for kind in ('spectral', 'group'):
            path = Path(job[kind+'_plan_path'])
            if json.loads(path.read_text()) != job[kind+'_plan']:
                raise ValueError('materialized source plan differs from specification')
            inputs[str(path.resolve())] = digest(path)
        spectra, electrons, source_paths = {}, {}, {}
        for request in job['spectral_plan']['requests']:
            root = Path(request['work']); same_mixture(root, job)
            data, grey = source(root); pin(root)
            text = (root/'source.txt').read_text()
            comp = text.split('No. Fraction Mass Fraction  At. No.  Chem. Sym.  Mat ID.\n')[1].split('Temperature grid')[0]
            helium = [line.split() for line in comp.splitlines()
                      if len(line.split()) == 5 and line.split()[3] == 'He']
            if len(helium) != 1:
                raise ValueError('missing spectral helium normalization')
            number, mass = map(float, helium[0][:2])
            atomic_mass = 4.002602*number/mass
            for key, values in data.items():
                if key in spectra:
                    raise ValueError('duplicate reference spectrum')
                spectra[key] = values
                electrons[key] = key[1]/(atomic_mass*1.66053906660e-24)*grey[key]['free_electrons_per_ion']
                source_paths[key] = str(root)
        seen = set()
        for request in job['group_plan']['requests']:
            root = Path(request['work']); same_mixture(root, job)
            boundaries, data, grey = source_groups(root); pin(root)
            for key, values in data.items():
                if key in seen or key not in spectra:
                    raise ValueError('duplicate or unmatched group state')
                seen.add(key)
                ne = electrons[key]
                if grey[key]['electron_density_cm3'] != ne:
                    raise ValueError('group and spectral electron counts differ')
                e = electron_moments(key[0]*KEV/KB, ne)
                cutoff = HBAR*math.sqrt(4*math.pi*4.80320471257e-10**2*ne/ME)/(key[0]*KEV)
                cutoff *= math.sqrt(e['plasma_frequency_squared_ratio'])
                full = [float(mean(spectra[key], key[0], cutoff, e['vstar_squared'], n))
                        for n in (8, 16)]
                grouped = [transport(boundaries, values, key[0], ne, n)
                           for n in (16, 32)]
                full_error = abs(full[0]/full[1]-1)
                group_error = abs(grouped[0]['rosseland_component_estimate']/grouped[1]['rosseland_component_estimate']-1)
                if full_error > 1e-7 or group_error > 1e-6:
                    raise ValueError('transport quadrature did not converge')
                value = grouped[1]['rosseland_component_estimate']
                rows.append({'X': job['X'], 'Z': job['Z'],
                    'temperature_keV': key[0], 'density_atomic_g_cm3': key[1],
                    'free_electron_density_cm3': ne, 'full_spectrum_source': source_paths[key],
                    'group_source': str(root), 'reference_rosseland_atomic_cm2_g': full[1],
                    'group_rosseland_atomic_cm2_g': value,
                    'component_relative_error': value/full[1]-1,
                    'full_spectrum_quadrature_relative_change': full_error,
                    'group_quadrature_relative_change': group_error,
                    'group_result': grouped[1]})
        if seen != set(spectra):
            raise ValueError('unmatched full-spectrum state')
        print(json.dumps({'X': job['X'], 'Z': job['Z'], 'states': len(seen),
            'maximum_relative_error': max(abs(r['component_relative_error']) for r in rows
                                         if r['X'] == job['X'] and r['Z'] == job['Z'])}), flush=True)
    maximum = max(abs(r['component_relative_error']) for r in rows)
    for name, h in inputs.items():
        if digest(Path(name)) != h:
            raise ValueError('input changed during comparison: '+name)
    result = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
        'constant_component_control': control, 'states': len(rows),
        'mixtures': len(spec['mixtures']), 'maximum_relative_error': maximum,
        'relative_comparison_criterion': spec['relative_comparison_criterion'],
        'group_representation_check_passed': maximum <= spec['relative_comparison_criterion'],
        'records': rows, 'input_sha256': inputs,
        'limitations': ['The group scattering fraction is inferred from group means, not independently supplied.',
            'These comparisons test frequency grouping at sampled source states, not interpolation between them.',
            'Both calculations retain the same finite-temperature free-electron dispersion and atomic spectra.',
            'Collisional and bound-electron dielectric changes are outside this numerical comparison.']}
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('records', 'input_sha256', 'limitations')}), flush=True)


if __name__ == '__main__':
    main()
