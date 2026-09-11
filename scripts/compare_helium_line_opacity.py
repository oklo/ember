#!/usr/bin/env python3
"""Compare the He I 4471-A omission control with unchanged source isotherms.

This measures sensitivity of sampled absorption means. It does not accept
the omitted-line opacity or establish a bound on an atmosphere solution.
Reference files correspond, in order, to the first diagnostic temperatures.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from compare_nongrey_opacity import compare
from nongrey_opacity import read_table
from prepare_nongrey_sources import digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('diagnostic', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--references', type=Path, nargs='+', required=True)
    a = p.parse_args()
    receipt_path = a.diagnostic/'completed.json'
    receipt = json.loads(receipt_path.read_text())
    if receipt.get('accepted_for_stellar_opacity') is not False:
        raise ValueError('explicit omission diagnostic required')
    prepared_path = a.diagnostic/'prepared.json'
    prepared = json.loads(prepared_path.read_text())
    if (not prepared.get('initialization_only')
            or digest(prepared_path) != receipt['prepared_sha256']):
        raise ValueError('diagnostic source identity differs')
    models = receipt['models']
    if not 1 <= len(a.references) <= len(models):
        raise ValueError('invalid number of original isotherms')
    inputs = {str(f.resolve()): digest(f) for f in
              [receipt_path, prepared_path, Path(__file__),
               Path(__file__).with_name('compare_nongrey_opacity.py'), *a.references]}
    for model in models:
        path = Path(model['table'])
        if digest(path) != model['table_sha256']:
            raise ValueError('diagnostic table changed')
        inputs[str(path.resolve())] = model['table_sha256']
    rows = []
    for model, old in zip(models, a.references):
        new = Path(model['table'])
        result = compare(new, old)
        candidate, reference = read_table(new), read_table(old)
        if candidate['log_electron_density'] != reference['log_electron_density']:
            raise ValueError('omission changed the source electron densities')
        x = np.array(candidate['log_opacity'], float).reshape(candidate['shape'])
        y = np.array(reference['log_opacity'], float).reshape(reference['shape'])
        changed = np.any(x != y, axis=(1, 2))
        wavelength = 2.997925e18/np.array(candidate['frequency'])
        result.update(negative_original_line_contribution_samples=int(np.count_nonzero(y < x)),
                      changed_samples=int(np.count_nonzero(x != y)),
                      affected_wavelength_range_A=([float(wavelength[changed].min()),
                          float(wavelength[changed].max())] if changed.any() else None),
                      tables={'candidate': str(new), 'reference': str(old)})
        rows.append(result)
    if any(digest(f) != checksum for f, checksum in inputs.items()):
        raise ValueError('an input changed during comparison')
    report = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
              'XH': receipt['recipe']['XH'], 'X3': receipt['recipe']['X3'],
              'source_recipe': receipt, 'input_sha256': inputs, 'controls': rows,
              'finite_diagnostic_temperatures_without_reference_K':
                  [r['temperature_K'] for r in models[len(a.references):]]}
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'comparisons': len(rows), 'maximum_Rosseland_relative_change':
        max(abs(r['means_relative_change'][1]) for c in rows for r in c['rows'])}))


if __name__ == '__main__':
    main()
