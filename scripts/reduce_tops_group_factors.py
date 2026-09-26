#!/usr/bin/env python3
"""Reduce verified TOPS groups with adaptive quadrature and exact result reuse.

Export the group transport mean, its dimensionless refractive ratio F, and
native uncut grey opacity times F. The last construction preserves the native
grey averaging convention; it is not an independent accuracy check. Outputs
remain candidates, with the original uncut recovery diagnostic retained.
"""
import argparse
from collections import Counter
import json
import math
from pathlib import Path

from audit_tops_spectral_means import digest
from audit_tops_electron_dispersion import KEV, KB
from reduce_tops_group_plan import uncut_mean
from tops_groups import source_groups, transport

QUADRATURE_CRITERION = 1e-6
UNCUT_RECOVERY_CRITERION = 1e-4


def verify(inputs):
    for name, expected in inputs.items():
        if digest(Path(name)) != expected:
            raise ValueError('input changed: ' + name)


def add_inputs(inputs, extra):
    for name, value in extra.items():
        if name in inputs and inputs[name] != value:
            raise ValueError('inconsistent input identity: ' + name)
        inputs[name] = value


def adaptive_state(boundaries, groups, t, rho, grey, source):
    previous = transport(boundaries, groups, t, grey['electron_density_cm3'], 16)
    previous_uncut = uncut_mean(boundaries, groups, t, 16)
    trials = []
    for order in (32, 64, 128, 256):
        result = transport(boundaries, groups, t, grey['electron_density_cm3'], order)
        value = result['rosseland_component_estimate']
        uncut = uncut_mean(boundaries, groups, t, order)
        if not all(math.isfinite(x) and x > 0 for x in (value, uncut)):
            raise ValueError('nonpositive or nonfinite group mean')
        change = abs(previous['rosseland_component_estimate'] / value - 1)
        uncut_change = abs(previous_uncut / uncut - 1)
        trials.append({'order': order, 'transport_relative_change': change,
                       'uncut_relative_change': uncut_change})
        if max(change, uncut_change) <= QUADRATURE_CRITERION:
            return {'temperature_keV': t, 'density_atomic_g_cm3': rho,
                    'rosseland_atomic_cm2_g': value,
                    'uncut_recombined_rosseland': uncut,
                    'uncut_source_rosseland': grey['rosseland'],
                    'uncut_recovery_relative_error': uncut / grey['rosseland'] - 1,
                    'quadrature_relative_change': change,
                    'uncut_quadrature_relative_change': uncut_change,
                    'quadrature_order': order, 'quadrature_trials': trials,
                    'free_electron_density_cm3': grey['electron_density_cm3'],
                    'cutoff_u': result['cutoff_u'], 'source': str(source)}
        previous, previous_uncut = result, uncut
    raise ValueError(f'quadrature did not converge at T={t}, rho={rho}: {trials}')


def validate_coordinates(records, job):
    tt = sorted({r['temperature_keV'] for r in records})
    rr = sorted({r['density_atomic_g_cm3'] for r in records})
    if tt != job['temperatures_keV'] or len(rr) != len(job['densities_g_cm3']):
        raise ValueError('source axes differ from plan')
    if any(abs(x / y - 1) > 5e-5 for x, y in zip(rr, job['densities_g_cm3'], strict=True)):
        raise ValueError('source densities differ from plan')
    if len(records) != len(tt) * len(rr):
        raise ValueError('incomplete source rectangle')


def table_text(records, tt, rr, x, z, field, description):
    lines = [f'EMBER_OPACITY_TABLE 2 1 {len(tt)} {len(rr)} {description}',
             ' '.join(format(math.log10(r), '.17g') for r in rr),
             ' '.join(format(math.log10(t * KEV / KB), '.17g') for t in tt),
             f'{x:.17g} {z:.17g}']
    prefixes = []
    for t in tt:
        row = sorted(r for temp, r in records if temp == t)
        if len(row) < 4 or row != rr[:len(row)]:
            raise ValueError('source support is not a contiguous density prefix')
        prefixes.append(len(row))
        values = [records[t, r][field] for r in row]
        if any(not math.isfinite(v) or v <= 0 for v in values):
            raise ValueError('invalid table value')
        lines.append(str(len(row)) + ' ' + ' '.join(format(math.log10(v), '.17g') for v in values))
    return '\n'.join(lines) + '\n', prefixes


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('plan', type=Path)
    p.add_argument('output_directory', type=Path)
    p.add_argument('report', type=Path)
    p.add_argument('--reuse-report', type=Path, action='append', default=[])
    a = p.parse_args()
    if a.output_directory.exists() or a.report.exists():
        raise FileExistsError('use new candidate output paths')
    plan = json.loads(a.plan.read_text())
    seed = Path(plan['baseline'])
    expected = json.loads((seed / 'receipt.json').read_text())
    paths = [a.plan, seed / 'receipt.json', Path(__file__),
             *(Path(__file__).with_name(n) for n in ('tops_groups.py',
               'reduce_tops_group_plan.py', 'audit_tops_electron_dispersion.py',
               'audit_tops_spectral_means.py', 'import_tops_composition.py',
               'fetch_tops_group_plan.py'))]
    inputs = {str(v.resolve()): digest(v) for v in paths}
    cached = {}
    for path in a.reuse_report:
        old = json.loads(path.read_text())
        if (old['X'], old['Z']) != (expected['X'], expected['Z']):
            raise ValueError('reused report has a different composition')
        add_inputs(inputs, old['input_sha256'])
        add_inputs(inputs, {str(path.resolve()): digest(path)})
        for record in old['records']:
            source = str(Path(record['source']).resolve())
            key = record['temperature_keV'], record['density_atomic_g_cm3']
            if key in cached.setdefault(source, {}):
                raise ValueError('duplicate cached source state')
            if record['quadrature_relative_change'] > QUADRATURE_CRITERION:
                raise ValueError('cached group quadrature fails')
            cached[source][key] = record
    verify(inputs)
    jobs = list(plan['requests'])
    for retained in plan.get('retained_requests', []):
        path = Path(retained['plan'])
        if digest(path) != retained['plan_sha256']:
            raise ValueError('retained plan changed')
        jobs.append(json.loads(path.read_text())['requests'][retained['request_index']])
        add_inputs(inputs, {str(path.resolve()): digest(path)})
    records = {}
    reused = 0
    for job in jobs:
        root = Path(job['work'])
        receipt = json.loads((root / 'receipt.json').read_text())
        if any(receipt[k] != expected[k] for k in ('X', 'Z', 'metals')):
            raise ValueError('source mixture differs from plane')
        source_inputs = {str((root / n).resolve()): digest(root / n)
                         for n in ('receipt.json', 'request.json', 'source.txt', 'recipe.json')}
        retained = cached.get(str(root.resolve()))
        if retained is not None:
            if any(inputs.get(name) != h for name, h in source_inputs.items()):
                raise ValueError('cached source is not fully identified')
            part = list(retained.values())
            reused += len(part)
        else:
            boundaries, groups, grey = source_groups(root)
            part = [adaptive_state(boundaries, data, t, rho, grey[t, rho], root)
                    for (t, rho), data in groups.items()]
        validate_coordinates(part, job)
        add_inputs(inputs, source_inputs)
        for original in part:
            record = dict(original)
            key = record['temperature_keV'], record['density_atomic_g_cm3']
            if key in records:
                raise ValueError('duplicate state in candidate plane')
            factor = record['rosseland_atomic_cm2_g'] / record['uncut_recombined_rosseland']
            if not math.isfinite(factor) or factor < 1 - 1e-12:
                raise ValueError('invalid refractive ratio')
            record['refractive_ratio'] = factor
            record['native_uncut_times_ratio_atomic_cm2_g'] = record['uncut_source_rosseland'] * factor
            records[key] = record
        print(json.dumps({'states': len(records), 'reused': reused, 'source': str(root)}), flush=True)
    tt = sorted({t for t, r in records})
    rr = sorted({r for t, r in records})
    if len(tt) < 4 or len(rr) < 4:
        raise ValueError('insufficient runtime axes')
    texts = {}
    for field, name, description in (
            ('rosseland_atomic_cm2_g', 'group_transport.dat', 'ATOMIC refractive group transport'),
            ('refractive_ratio', 'factor.dat', 'DIMENSIONLESS refractive ratio; not standalone opacity'),
            ('native_uncut_times_ratio_atomic_cm2_g', 'normalized_transport.dat',
             'ATOMIC native uncut grey opacity times group refractive ratio')):
        texts[name], prefixes = table_text(records, tt, rr, expected['X'], expected['Z'], field, description)
    verify(inputs)
    a.output_directory.mkdir(parents=True)
    for name, contents in texts.items():
        (a.output_directory / name).write_text(contents)
    maximum = max(abs(r['uncut_recovery_relative_error']) for r in records.values())
    result = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
              'X': expected['X'], 'Z': expected['Z'], 'states': len(records),
              'reused_states': reused, 'newly_integrated_states': len(records) - reused,
              'temperatures_keV': tt, 'densities_atomic_g_cm3': rr,
              'density_prefix_sizes': prefixes,
              'quadrature_relative_criterion': QUADRATURE_CRITERION,
              'maximum_quadrature_relative_change': max(r['quadrature_relative_change'] for r in records.values()),
              'new_quadrature_orders': dict(Counter(r['quadrature_order'] for r in records.values() if 'quadrature_order' in r)),
              'uncut_recovery_relative_criterion': UNCUT_RECOVERY_CRITERION,
              'maximum_uncut_recovery_relative_error': maximum,
              'uncut_recovery_check_passed': maximum <= UNCUT_RECOVERY_CRITERION,
              'normalization_is_an_independent_accuracy_check': False,
              'limitations': ['Group scattering fraction remains a narrow-group approximation.',
                              'Quadrature convergence does not establish frequency-group accuracy.',
                              'Cold dense states are not validated by hot full-spectrum comparisons.',
                              'Interpolation and completed stellar input coverage require separate checks.'],
              'output_sha256': {str((a.output_directory / n).resolve()): digest(a.output_directory / n) for n in texts},
              'input_sha256': inputs, 'records': list(records.values())}
    a.report.parent.mkdir(parents=True, exist_ok=True)
    a.report.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in (
        'records', 'input_sha256', 'temperatures_keV', 'densities_atomic_g_cm3', 'density_prefix_sizes')}), flush=True)


if __name__ == '__main__':
    main()
