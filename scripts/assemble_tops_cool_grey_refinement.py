#!/usr/bin/env python3
"""Refine a cold grey diagnostic using means accompanying verified groups.

Retain every supported original grey value in the common cold rectangle.
Check all overlapping means, and use only new density coordinates from the
group reduction. No duplicate grey requests or interpolated source values.
This diagnostic assembly is not an accepted stellar opacity family.
"""
import argparse
import json
import math
from pathlib import Path

from audit_refractive_opacity_family_v2 import read_table
from audit_tops_electron_dispersion import KEV, KB
from fetch_tops_composition import digest
from reduce_tops_group_factors import add_inputs, verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('base_report', 'factor_report', 'reference_report', 'output', 'report'):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    if a.output.exists() or a.report.exists():
        raise FileExistsError('preserve completed diagnostic assemblies')
    inputs = {}
    def load(path):
        r = json.loads(path.read_text())
        add_inputs(inputs, {str(path.resolve()): digest(path)})
        add_inputs(inputs, r['input_sha256'])
        add_inputs(inputs, r['output_sha256'])
        return r
    base, factor, reference = map(load, (a.base_report, a.factor_report, a.reference_report))
    for name in ('assemble_tops_cool_grey_refinement.py', 'audit_refractive_opacity_family_v2.py',
                 'audit_tops_electron_dispersion.py', 'reduce_tops_group_factors.py',
                 'fetch_tops_composition.py'):
        path = Path(__file__).with_name(name)
        add_inputs(inputs, {str(path.resolve()): digest(path)})
    verify(inputs)
    x, z = factor['X'], factor['Z']
    if any((r['X'], r['Z']) != (x, z) for r in (base, reference)):
        raise ValueError('composition mismatch')
    table = read_table(Path(next(iter(base['output_sha256']))))
    if table['X'] != [x] or table['Z'] != z:
        raise ValueError('base table composition mismatch')
    tt = factor['temperatures_keV']
    if any(t >= .025 or t not in base['temperatures_keV'] for t in tt):
        raise ValueError('expected a cold native-temperature diagnostic')
    base_rr = [r for r in base['densities_atomic_g_cm3'] if r <= max(factor['densities_atomic_g_cm3'])]
    rr = sorted(set(base_rr) | set(factor['densities_atomic_g_cm3']))
    new_rr = set(rr)-set(base_rr)
    refs = [r for r in reference['records'] if r['temperature_keV'] in tt]
    if any(abs(rho/r['density_atomic_g_cm3']-1) <= 5e-5 for rho in new_rr for r in refs):
        raise ValueError('new grey density coincides with a withheld reference')
    logs = {}
    for t in tt:
        i = base['temperatures_keV'].index(t)
        if abs(table['logT'][i]-math.log10(t*KEV/KB)) > 1e-12:
            raise ValueError('base temperature metadata mismatch')
        for j, rho in enumerate(base_rr):
            if abs(table['logrho'][j]-math.log10(rho)) > 1e-12:
                raise ValueError('base density metadata mismatch')
            if j < len(table['planes'][x][i]):
                logs[t, rho] = table['planes'][x][i][j]
    retained = dict(logs)
    overlaps = []
    added = 0
    for row in factor['records']:
        key = row['temperature_keV'], row['density_atomic_g_cm3']
        v = row['uncut_source_rosseland']
        if not math.isfinite(v) or v <= 0:
            raise ValueError('invalid source mean')
        if key in logs:
            error = abs(10**logs[key]/v-1)
            if error > 1e-12: raise ValueError('group and grey sources disagree')
            overlaps.append(error)
        elif key[1] in base_rr:
            raise ValueError('group source would fill an excluded original state')
        else:
            logs[key] = math.log10(v)
            added += 1
    if any(logs[k] != v for k, v in retained.items()):
        raise ValueError('original source log value changed')
    lines = [f'EMBER_OPACITY_TABLE 2 1 {len(tt)} {len(rr)} TOPS cold uncut grey diagnostic; group-accompanying density refinement',
             ' '.join(format(math.log10(r), '.17g') for r in rr),
             ' '.join(format(math.log10(t*KEV/KB), '.17g') for t in tt), f'{x:.17g} {z:.17g}']
    prefixes = []
    for t in tt:
        present = sorted(r for temp, r in logs if temp == t)
        if len(present) < 4 or present != rr[:len(present)]:
            raise ValueError('source support does not form a density prefix')
        prefixes.append(len(present))
        lines.append(str(len(present))+' '+' '.join(format(logs[t,r], '.17g') for r in present))
    verify(inputs)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text('\n'.join(lines)+'\n')
    check = read_table(a.output)
    for (t, rho), value in logs.items():
        if check['planes'][x][tt.index(t)][rr.index(rho)] != value:
            raise ValueError('serialized value changed')
    result = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
              'X': x, 'Z': z, 'temperatures_keV': tt, 'densities_atomic_g_cm3': rr,
              'density_prefix_sizes': prefixes, 'states': len(logs),
              'original_values_retained_exactly': len(retained), 'new_source_values': added,
              'overlap_comparisons': len(overlaps), 'maximum_overlap_relative_error': max(overlaps),
              'reference_density_points_added': 0, 'input_sha256': inputs,
              'output_sha256': {str(a.output.resolve()): digest(a.output)}}
    a.report.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('scope','input_sha256','temperatures_keV','densities_atomic_g_cm3','density_prefix_sizes')}))


if __name__ == '__main__':
    main()
