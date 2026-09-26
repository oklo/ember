#!/usr/bin/env python3
"""Measure a metal atmosphere response from completed sources at matched nodes.

The reference values must exist in the retained atmosphere table. Reduced-Z
columns are independently revalidated. Missing reference or reduced-Z states
remain absent; no missing source point is filled by interpolation.
"""
import argparse
from datetime import datetime, timezone
import itertools
import json
import math
from pathlib import Path

from assemble_nongrey_grid import load_continuation
from audit_atmosphere_retention import read
from prepare_nongrey_sources import digest
from write_scientific_result import write_result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('reference', type=Path)
    p.add_argument('plan', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('report', type=Path)
    a = p.parse_args()
    if a.output.exists() or a.report.exists():
        raise FileExistsError('preserve existing response tables and reports')
    plan = json.loads(a.plan.read_text())
    header, _, reference = read(a.reference)
    z0 = sum(map(float, header['metals']))
    z1 = plan['source_Z']
    if (abs(z0-plan['reference_Z']) > 1e-14 or not 0 <= z1 < z0 < 1
            or header['basis'] != ['baryon_mass']):
        raise ValueError('invalid reference or source composition')
    axes = [plan['hydrogen'], [math.log10(v) for v in plan['teff_K']], plan['log_g']]
    for axis in axes:
        if len(axis) < 2 or sorted(set(axis)) != axis or not all(map(math.isfinite, axis)):
            raise ValueError('invalid response axis')
    if not 0 <= plan['maximum_helium3'] < .01:
        raise ValueError('this approximation requires a declared small isotope range')
    paths = [Path(__file__).resolve(), a.reference.resolve(), a.plan.resolve()]
    hashes = {str(path): digest(path) for path in paths}
    # Additional reference sources supply actual He3=0 nodes only. They do
    # not extend the isotope coverage of the atmosphere used by evolution.
    def same_composition(spec, z):
        return (abs(sum(spec['metals'])-z) <= 1e-14
                and spec['tau'] == float(header['tau'][0])
                and spec['approximation'] == header['approximation'][0]
                and all(abs(v/z-r) <= 1e-12 for v, r in
                        zip(spec['metals'], metal_ratios, strict=True)))

    def add_reference(key, value):
        previous = reference.get(key, (0,))
        if previous[0] and previous != value:
            raise ValueError('conflicting reference source coordinate')
        if value[0]:
            reference[key] = value

    records = {}
    metal_ratios = [float(v)/z0 for v in header['metals']]
    for path in plan.get('reference_tables', []):
        extra_header, _, extra = read(Path(path))
        if any(extra_header[k] != header[k] for k in
               ['approximation', 'basis', 'tau', 'metals']):
            raise ValueError('additional reference table physics differs')
        hashes[str(Path(path).resolve())] = digest(path)
        for key, value in extra.items():
            if key[1] == 0:
                add_reference(key, value)
    for path in plan.get('reference_continuations', []):
        key, state, spec, _, dependencies = load_continuation(Path(path))
        if key[1] != 0 or not same_composition(spec, z0):
            raise ValueError('reference source composition or physics differs')
        add_reference((key[0], 0., math.log10(key[2]), key[3]),
                      (1, math.log10(state['T']), math.log10(state['Pgas'])))
        hashes.update(dependencies)
    for path in plan['continuations']:
        key, state, spec, _, dependencies = load_continuation(Path(path))
        if (key[1] != 0 or abs(sum(spec['metals'])-z1) > 1e-14
                or spec['tau'] != float(header['tau'][0])
                or spec['approximation'] != header['approximation'][0]
                or any(abs(v/z1-r) > 1e-12 for v, r in zip(spec['metals'], metal_ratios, strict=True))):
            raise ValueError('source composition, matching depth or approximation differs')
        raw_key = (key[0], math.log10(key[2]), key[3])
        if raw_key in records:
            raise ValueError('duplicate response source coordinate')
        records[raw_key] = dict(state=state, path=path)
        hashes.update(dependencies)
    lines = ['EMBER_METAL_ATMOSPHERE_RESPONSE 1',
             'source "Measured GS98 metal response from validated LTE gas atmospheres"',
             'approximation "Interpolate log matching-temperature and gas-pressure changes in Z; fixed GS98 ratios; metal and small-He3 responses assumed separable"',
             'basis baryon_mass', f'reference_Z {plan["reference_Z"]:.17g}',
             f'source_Z {z1:.17g}', f'maximum_helium3 {plan["maximum_helium3"]:.17g}',
             f'tau {float(header["tau"][0]):.17g}']
    for name, axis in zip(['hydrogen', 'log_teff', 'log_g'], axes):
        lines.append(f'{name} {len(axis)} '+' '.join(format(v, '.17g') for v in axis))
    lines.append('data')
    comparisons, missing = [], []
    for key in itertools.product(*axes):
        reference_key = (key[0], 0., key[1], key[2])
        ref = reference.get(reference_key, (0,))
        if not ref[0] or key not in records:
            lines.append('0')
            missing.append(dict(coordinates=key, missing_reference=not bool(ref[0]),
                                missing_reduced_Z=key not in records))
            continue
        low = records[key]['state']
        delta = [math.log(low['T'])-math.log(10)*ref[1],
                 math.log(low['Pgas'])-math.log(10)*ref[2]]
        lines.append('1 '+' '.join(format(v, '.17g') for v in delta))
        comparisons.append(dict(coordinates=[key[0], 10**key[1], key[2]],
                                source=records[key]['path'], delta_logT=delta[0],
                                delta_logPgas=delta[1],
                                reference_T=10**ref[1], reference_Pgas=10**ref[2],
                                source_T=low['T'], source_Pgas=low['Pgas']))
    if not comparisons:
        raise ValueError('no matched reference/source states')
    for path, expected in hashes.items():
        if digest(path) != expected:
            raise ValueError('input changed during calibration: '+path)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text('\n'.join(lines)+'\n')
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  outcome='measured_unselected_response', selected_for_evolution=False,
                  input_sha256=hashes, output_sha256=digest(a.output),
                  response=str(a.output.resolve()), comparisons=comparisons, missing=missing,
                  limitations=['Calibrated at source nodes only; intermediate-Z/H and small-He3 assumptions require stellar sensitivity checks.',
                               'Reference table physics and its previously measured interpolation uncertainty remain unchanged.',
                               'This file does not select or evolve a stellar model.'])
    write_result(a.report, report)
    print(len(comparisons), 'matched source states;', len(missing), 'missing')


if __name__ == '__main__':
    main()
