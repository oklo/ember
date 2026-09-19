#!/usr/bin/env python3
"""Check completed cached EOS overlap states without waiting for whole cubes.

Missing isotherms are not failures and are not accepted as completed sources.
This check cannot replace the final complete-cube merge or EOS accuracy audit.
"""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path

from metal_eos_composition import mixture
from merge_metal_eos_density_sources import overlap_difference


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['previous', 'addition', 'output']:
        p.add_argument(name, type=Path)
    a = p.parse_args()
    roots = [a.previous, a.addition]
    specs = [json.loads((root/'specification.json').read_text()) for root in roots]
    for key in ['hydrogen', 'helium3', 'probe_sha256', 'source_archive_sha256', 'source_options']:
        if specs[0][key] != specs[1][key]:
            raise ValueError('source identities differ')
    if specs[0]['logQ'][-1] != specs[1]['logQ'][0]:
        raise ValueError('source density overlap differs')
    inputs = {str((root/'specification.json').resolve()): sha((root/'specification.json').read_bytes()) for root in roots}
    for path in [Path(__file__), Path(__file__).with_name('merge_metal_eos_density_sources.py'),
                 Path(__file__).with_name('metal_eos_composition.py')]:
        inputs[str(path.resolve())] = sha(path.read_bytes())
    records = []
    for jt, t in enumerate(specs[1]['logT']):
        if t not in specs[0]['logT']:
            raise ValueError('source temperatures do not overlap exactly')
        it = specs[0]['logT'].index(t)
        for plane, (x, y) in enumerate((x, y) for x in specs[0]['hydrogen'] for y in specs[0]['helium3']):
            paths = [root/f'plane-{plane:03d}'/f'temperature-{k:03d}.json.gz' for root, k in zip(roots, [it, jt])]
            if not all(path.exists() for path in paths):
                continue
            m = mixture(x, y); response = []
            for path, spec, index in zip(paths, specs, [-1, 0], strict=True):
                blob = path.read_bytes(); inputs[str(path.resolve())] = sha(blob)
                saved = json.loads(gzip.decompress(blob))
                request = ' '.join(map(str, m['eps']))+'\n'+' '.join(map(str, spec['source_options']))+'\n'+''.join(
                    f'{math.log(m["source_mass_scale"])+math.log(10)*(q+1.5*(t-6)):.17g} {math.log(10)*t:.17g}\n'
                    for q in spec['logQ'])
                if (saved['input'] != request
                        or saved['input_sha256'] != sha((spec['probe_sha256']+request).encode())
                        or len(saved['data']) != len(spec['logQ'])
                        or any(len(row) != 22 or not all(map(math.isfinite, row)) for row in saved['data'])):
                    raise ValueError('source cache differs from its requested state')
                response.append(saved['data'][index])
            left, right = response
            difference = overlap_difference(left, right)
            records.append(dict(hydrogen=x, helium3=y, logT=t, source_flags=[left[0], right[0]],
                                source_iterations=[left[1], right[1]], maximum_scaled_difference=difference))
    converged = [r for r in records if r['maximum_scaled_difference'] is not None]
    if not converged:
        raise ValueError('no completed converged overlap')
    maximum = max(r['maximum_scaled_difference'] for r in converged)
    if any(sha(Path(path).read_bytes()) != value for path, value in inputs.items()):
        raise ValueError('a completed source input changed during the comparison')
    report = dict(scope=__doc__, input_sha256=inputs, checked=len(records), converged=len(converged),
                  different_iteration_counts=sum(r['source_iterations'][0] != r['source_iterations'][1] for r in converged),
                  maximum_scaled_physical_difference=maximum, criterion=1e-8, passed=maximum <= 1e-8,
                  records=records, completed_family_accepted=False)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ['scope', 'input_sha256', 'records']}))
    if not report['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
