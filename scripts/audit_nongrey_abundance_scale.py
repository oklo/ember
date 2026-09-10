#!/usr/bin/env python3
"""Test the native opacity source's abundance-scale invariance at fixed T/rho.

All abundances, including trace placeholders, are multiplied by one common
factor. Baryonic fractions and element masses are unchanged. This is an
offline source diagnostic; it neither installs a different normalization nor
substitutes a small hydrogen abundance for the zero-hydrogen endpoint.
"""
import argparse
import json
import math
from pathlib import Path

from generate_nongrey_grid import opacity_inputs, execute, sequence
from nongrey_opacity import validate_table
from prepare_nongrey_sources import digest, data_digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['prepared', 'specification', 'work', 'output']:
        p.add_argument(name, type=Path)
    p.add_argument('--temperature', type=float, required=True)
    p.add_argument('--hydrogen', type=float, required=True)
    p.add_argument('--helium3', type=float, required=True)
    p.add_argument('--scales', type=float, nargs='+', default=[1., .15, 3.])
    a = p.parse_args()
    if a.work.exists() or a.output.exists():
        raise FileExistsError('use fresh source work and report paths')
    if a.temperature <= 0 or not math.isfinite(a.temperature):
        raise ValueError('invalid temperature')
    if len(a.scales) < 2 or a.scales[0] != 1. or len(set(a.scales)) != len(a.scales):
        raise ValueError('start with the unscaled control, followed by unique scales')
    if any(v <= 0 or not math.isfinite(v) for v in a.scales):
        raise ValueError('invalid abundance scale')
    prepared = json.loads(a.prepared.read_text())
    spec = json.loads(a.specification.read_text())
    if 'depletion' in prepared or 'initialization_only' in prepared:
        raise ValueError('canonical gas opacity source required')
    dependencies = {str(p.resolve()): digest(p) for p in
                    [a.prepared, a.specification, Path(prepared['synspec']),
                     *map(Path, prepared['line_lists'])]}
    if digest(prepared['synspec']) != prepared['executables']['synspec']:
        raise ValueError('source executable changed')
    if any(digest(p) != h for p, h in zip(prepared['line_lists'], prepared['line_list_sha256'], strict=True)):
        raise ValueError('line list changed')
    if data_digest(Path(prepared['synple'])/'data') != prepared['data_sha256']:
        raise ValueError('source material inputs changed')
    a.work.mkdir(parents=True)
    (a.work/'request.json').write_text(json.dumps(vars(a), default=str, indent=2)+'\n')
    reference = None
    results = []
    for i, scale in enumerate(a.scales):
        directory = a.work/f'scale-{i:03d}'
        abundance, masses = opacity_inputs(directory, prepared, spec, a.hydrogen,
                                            a.helium3, a.temperature)
        scaled = [v*scale for v in abundance]
        if any(v <= 0 or v > 1e6 for v in scaled):
            raise ValueError('scaled value invokes a different native input convention')
        lines = (directory/'fort.5').read_text().splitlines()
        for j, value in enumerate(scaled):
            row = lines[5+j].split()
            if len(row) != 3 or float(row[1]) != abundance[j]:
                raise ValueError('unrecognized native element input')
            row[1] = format(value, '.17e')
            lines[5+j] = ' '.join(row)
        (directory/'fort.5').write_text('\n'.join(lines)+'\n')
        execute(prepared['synspec'], directory, ['fort.63', 'fort.29'])
        table = validate_table(directory/'fort.63', scaled, [a.temperature], sequence(spec['log_density']))
        if reference is None:
            reference = table
        if any(table[k] != reference[k] for k in ['shape', 'frequency', 'log_temperature', 'log_density', 'flags']):
            raise ValueError('comparison changed source axes or opacity processes')
        differences = {}
        for key in ['log_electron_density', 'log_opacity']:
            values = [v-r for v, r in zip(table[key], reference[key], strict=True)]
            differences[key] = {'max_abs_log_difference': max(map(abs, values)),
                                'relative_difference_range': [math.expm1(min(values)), math.expm1(max(values))]}
        original_mass = sum(v*m for v, m in zip(abundance, masses, strict=True))
        scaled_mass = sum(v*m for v, m in zip(scaled, masses, strict=True))
        mass_error = max(abs(v*m/scaled_mass-r*m/original_mass)
                         for v, r, m in zip(scaled, abundance, masses, strict=True))
        record = {'scale': scale, 'mass_fraction_max_abs_difference': mass_error,
                  'shape': table['shape'], 'differences': differences,
                  'work': str(directory.resolve()), 'completed_sha256': digest(directory/'completed.json')}
        results.append(record)
        print(json.dumps(record), flush=True)
    if any(digest(p) != h for p, h in dependencies.items()):
        raise ValueError('source inputs changed during diagnostic')
    report = {'scope': __doc__, 'XH': a.hydrogen, 'X3': a.helium3, 'temperature_K': a.temperature,
              'results': results, 'input_files_sha256': dependencies, 'script_sha256': digest(__file__),
              'interpretation': 'Agreement here would test numerical normalization at this isotherm only; zero-hydrogen chemistry and TLUSTY atmosphere normalization still require separate verification.'}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
