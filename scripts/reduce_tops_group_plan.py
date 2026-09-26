#!/usr/bin/env python3
"""Reduce verified TOPS groups into a separate single-composition opacity table.

Check recovery of the source's uncut grey mean before applying the refractive
transport calculation. Every output value comes from a source state. This
operation does not accept the table's interpolation or install stellar input.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from audit_tops_spectral_means import digest
from audit_tops_electron_dispersion import KEV, KB, RW
from tops_groups import source_groups, transport


def uncut_mean(boundaries, groups, temperature, order=32):
    lo, hi = boundaries[:-1]/temperature, boundaries[1:]/temperature
    nodes, weights = leggauss(order)
    u = (lo[:, None]+hi[:, None])/2+(hi-lo)[:, None]/2*nodes
    wr = u**4*np.exp(-u)/(-np.expm1(-u))**2
    integral = np.sum((hi-lo)[:, None]/2*weights*wr/groups[:, 1, None])
    return float(RW/integral)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('plan', type=Path)
    p.add_argument('table', type=Path)
    p.add_argument('report', type=Path)
    a = p.parse_args()
    if a.table.exists() or a.report.exists():
        raise FileExistsError('use new candidate output paths')
    plan = json.loads(a.plan.read_text())
    seed = Path(plan['baseline'])
    expected = json.loads((seed/'receipt.json').read_text())
    inputs = {str(v.resolve()):digest(v) for v in
              [a.plan, seed/'receipt.json', Path(__file__),
               *(Path(__file__).with_name(n) for n in ['tops_groups.py',
                   'audit_tops_electron_dispersion.py', 'audit_tops_spectral_means.py',
                   'import_tops_composition.py', 'fetch_tops_group_plan.py'])]}
    records = {}
    jobs = list(plan['requests'])
    for retained in plan.get('retained_requests', []):
        path = Path(retained['plan'])
        if digest(path) != retained['plan_sha256']:
            raise ValueError('retained source plan changed')
        old = json.loads(path.read_text())
        jobs.append(old['requests'][retained['request_index']])
        inputs[str(path.resolve())] = digest(path)
    for job in jobs:
        root = Path(job['work'])
        receipt = json.loads((root/'receipt.json').read_text())
        if any(receipt[k] != expected[k] for k in ('X', 'Z', 'metals')):
            raise ValueError('group mixture differs from the candidate plane')
        boundaries, groups, grey = source_groups(root)
        if sorted({k[0] for k in groups}) != job['temperatures_keV']:
            raise ValueError('group temperatures differ from the plan')
        actual_rho = sorted({k[1] for k in groups})
        if len(actual_rho) != len(job['densities_g_cm3']) or any(
                abs(x/y-1) > 5e-5 for x,y in zip(actual_rho,job['densities_g_cm3'],strict=True)):
            raise ValueError('group densities differ from the plan')
        for name in ['receipt.json', 'request.json', 'source.txt', 'recipe.json']:
            path = root/name; inputs[str(path.resolve())] = digest(path)
        for key, data in groups.items():
            if key in records:
                raise ValueError('duplicate group state in candidate plane')
            t, rho = key
            pair = [transport(boundaries, data, t, grey[key]['electron_density_cm3'], n)
                    for n in (16, 32)]
            value = pair[1]['rosseland_component_estimate']
            quadrature = abs(pair[0]['rosseland_component_estimate']/value-1)
            if not math.isfinite(value) or value <= 0 or quadrature > 1e-6:
                raise ValueError('invalid or unconverged group transport')
            uncut = uncut_mean(boundaries, data, t)
            records[key] = {'temperature_keV':t, 'density_atomic_g_cm3':rho,
                'rosseland_atomic_cm2_g':value,
                'uncut_recombined_rosseland':uncut,
                'uncut_source_rosseland':grey[key]['rosseland'],
                'uncut_recovery_relative_error':uncut/grey[key]['rosseland']-1,
                'quadrature_relative_change':quadrature,
                'free_electron_density_cm3':grey[key]['electron_density_cm3'],
                'cutoff_u':pair[1]['cutoff_u'], 'source':str(root)}
            if len(records) % 256 == 0:
                print(json.dumps({'reduced_states':len(records)}),flush=True)
    tt = sorted({t for t,r in records})
    rr = sorted({r for t,r in records})
    if len(tt) < 4 or len(rr) < 4:
        raise ValueError('insufficient runtime axes')
    lines = [f'EMBER_OPACITY_TABLE 2 1 {len(tt)} {len(rr)} ATOMIC groups; finite-temperature free-electron dispersion; atomic mass basis',
             ' '.join(format(math.log10(r),'.17g') for r in rr),
             ' '.join(format(math.log10(t*KEV/KB),'.17g') for t in tt),
             f'{expected["X"]:.17g} {expected["Z"]:.17g}']
    prefixes = []
    for t in tt:
        row = sorted(r for temp,r in records if temp == t)
        if len(row) < 4 or row != rr[:len(row)]:
            raise ValueError('source support is not a contiguous density prefix')
        prefixes.append(len(row))
        lines.append(str(len(row))+' '+' '.join(format(math.log10(records[t,r]['rosseland_atomic_cm2_g']),'.17g') for r in row))
    maximum = max(abs(r['uncut_recovery_relative_error']) for r in records.values())
    passed = maximum <= 1e-4
    for name,h in inputs.items():
        if digest(Path(name)) != h:
            raise ValueError('input changed during reduction')
    a.table.parent.mkdir(parents=True,exist_ok=True)
    a.table.write_text('\n'.join(lines)+'\n')
    result = {'scope':__doc__, 'accepted_for_stellar_opacity':False,
        'X':expected['X'], 'Z':expected['Z'], 'states':len(records),
        'temperatures_keV':tt, 'densities_atomic_g_cm3':rr, 'density_prefix_sizes':prefixes,
        'uncut_recovery_relative_criterion':1e-4,
        'maximum_uncut_recovery_relative_error':maximum,
        'uncut_recovery_check_passed':passed,
        'maximum_quadrature_relative_change':max(r['quadrature_relative_change'] for r in records.values()),
        'table':str(a.table.resolve()),'table_sha256':digest(a.table),
        'input_sha256':inputs,'records':list(records.values())}
    a.report.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('records','input_sha256','temperatures_keV','densities_atomic_g_cm3','density_prefix_sizes')}),flush=True)
    if not passed:
        raise ValueError('frequency groups do not recover the source uncut mean')


if __name__ == '__main__':
    main()
