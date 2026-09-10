#!/usr/bin/env python3
"""Summarize recorded F77 center-stability transitions, retaining their brackets.

The diagnostic adds only a print call to a separate source copy. Original
output bytes must agree with an unmodified-source control where supplied.
The first nonconvective central flag and the gradient crossing are separate
events because F77 has ten-percent convective-flag hysteresis.
"""
import argparse
import difflib
import hashlib
import json
import math
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def records(path):
    rows = []
    keys = ['model', 'age_yr', 'central_X', 'central_Y3', 'central_convective_flag',
            'connected_nonconvective_core_mass_fraction', 'central_gradient_ratio',
            'central_T_K', 'central_rho', 'mesh_points']
    for line in path.read_text().splitlines():
        if not line.startswith('EMBER_CORE '):
            continue
        values = list(map(float, line.split()[1:]))
        if len(values) != len(keys) or not all(math.isfinite(v) for v in values):
            raise ValueError('invalid central diagnostic record')
        row = dict(zip(keys, values))
        for name in ['model', 'central_convective_flag', 'mesh_points']:
            row[name] = int(row[name])
        rows.append(row)
    if not rows or any(b['model'] <= a['model'] or b['age_yr'] <= a['age_yr']
                       for a, b in zip(rows, rows[1:])):
        raise ValueError('empty or nonmonotonic F77 diagnostic history')
    return rows


def transitions(rows, quantity, threshold):
    return [{'direction': 'down' if quantity(b) < threshold else 'up',
             'before': a, 'after': b}
            for a, b in zip(rows, rows[1:])
            if (quantity(a) < threshold) != (quantity(b) < threshold)]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('work', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    inputs = {str(f.resolve()): sha(f) for f in a.work.rglob('*') if f.is_file()}
    runs = []
    for name, control, opacity in [('diagnostic', 'baseline', 'Ferguson; 6000-model control'),
                                   ('ajr-20000', 'ajr-baseline', 'AJR1983; historical opacity option'),
                                   ('ferguson-20000', None, 'Ferguson; extended current preset')]:
        receipt = json.loads((a.work/(name+'.receipt.json')).read_text())
        if receipt['returncode'] not in [0, 2]:
            raise ValueError('unexpected F77 exit status')
        output = a.work/(name+'.out')
        rows = records(output)
        agreement = None
        if control:
            original = a.work/(control+'.out')
            original_receipt = json.loads((a.work/(control+'.receipt.json')).read_text())
            stripped = b''.join(line for line in output.read_bytes().splitlines(keepends=True)
                                if not line.startswith(b'EMBER_CORE '))
            agreement = not original_receipt['returncode'] and original.read_bytes() == stripped
            if not agreement:
                raise ValueError('diagnostic changed original F77 output')
        # Retain early transients too; identify later main-sequence onset from
        # the actual sequence instead of an interpolated or guessed age.
        runs.append({'name': name, 'opacity': opacity, 'receipt': receipt,
                     'run_completed_without_error': receipt['returncode'] == 0,
                     'history_scope': 'Prints occur only after accepted solves; a later failure does not establish the requested endpoint.',
                     'unmodified_output_bytes_identical': agreement,
                     'recorded_models': len(rows), 'first': rows[0], 'final': rows[-1],
                     'central_flag_transitions': transitions(rows, lambda r: r['central_convective_flag'], .5),
                     'central_gradient_unity_crossings': transitions(rows, lambda r: r['central_gradient_ratio'], 1.)})
    source = a.work/'baseline.f'
    diagnostic = a.work/'diagnostic.f'
    patch = ''.join(difflib.unified_diff(source.read_text().splitlines(keepends=True),
                                       diagnostic.read_text().splitlines(keepends=True),
                                       fromfile='baseline.f', tofile='diagnostic.f'))
    if any(sha(Path(path)) != checksum for path, checksum in inputs.items()):
        raise ValueError('recorded input or output changed during summary')
    report = {'scope': __doc__, 'runs': runs, 'diagnostic_patch': patch,
              'compiler_commands': [['/opt/homebrew/bin/gfortran', '-O2', name+'.f', '-o', name]
                                    for name in ['baseline', 'diagnostic']],
              'decks': {f.name: f.read_text() for f in a.work.glob('*.in')},
              'reproduction': 'Compile the preserved source and patched copy with the listed commands; run with each recorded deck from a directory containing the pinned scvh/ and ferguson/ inputs. Compare original output after removing EMBER_CORE lines.',
              'limitations': 'These reproduce the currently retained F77 source with explicitly named opacity options. They do not recover an unarchived historical executable or deck, determine a converged transition age, or transfer that age to Ember. The core diagnostic includes conductive transport through the F77 effective opacity.',
              'input_files_sha256': inputs, 'script_sha256': sha(Path(__file__))}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps([{'name': r['name'], 'central_flag_transitions': r['central_flag_transitions'],
                       'final': r['final']} for r in runs]))


if __name__ == '__main__':
    main()
