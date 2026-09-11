#!/usr/bin/env python3
"""Join EOS density ranges at one shared density for each temperature.

Existing source rows are retained exactly. The shared coordinate checks the
effect of starting a new source process at high density; each mixture's
converged responses must agree before its merged source is written. This
prepares raw material data, not an accepted thermodynamic table or a stellar restart.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
import math
from pathlib import Path

from metal_eos_composition import mixture


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(value):
    return (json.dumps(value, separators=(',', ':'), allow_nan=False)+'\n').encode()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['previous', 'addition', 'output']:
        p.add_argument(name, type=Path)
    p.add_argument('--jobs', type=int, default=2)
    a = p.parse_args()
    if not 1 <= a.jobs <= 4 or a.output.exists():
        raise ValueError('use a fresh output directory and one to four workers')
    specs = [root/'specification.json' for root in [a.previous, a.addition]]
    old_spec, extra_spec = [json.loads(path.read_text()) for path in specs]
    options = old_spec.get('source_options', [3, 1, -2])
    if options not in [[3, 1, -2], [3, 223, -2]]:
        raise ValueError('unsupported source electron treatment')
    for key in ['hydrogen', 'helium3', 'logT', 'probe_sha256', 'source_archive_sha256']:
        if old_spec[key] != extra_spec[key]:
            raise ValueError('source composition, temperature or executable differs')
    if (extra_spec.get('source_options', [3, 1, -2]) != options
            or old_spec.get('source_radiation_included', True) != (options == [3, 1, -2])
            or extra_spec.get('source_radiation_included', True) != (options == [3, 1, -2])):
        raise ValueError('source physics or radiation treatment differs')
    old_q, extra_q = old_spec['logQ'], extra_spec['logQ']
    if len(old_q) < 5 or len(extra_q) < 5:
        raise ValueError('source density axes are too short')
    q = old_q+extra_q[1:]
    step = old_q[1]-old_q[0]
    if (not step > 0 or old_q[-1] != extra_q[0]
            or any(not math.isfinite(v) for v in q)
            or any(abs((right-left)/step-1) > 1e-10 for left, right in zip(q, q[1:]))):
        raise ValueError('density grids must have exactly one shared endpoint and equal spacing')
    inputs = {str(path.resolve()): digest(path) for path in [*specs, Path(__file__)]}
    keys = [(x, y) for x in old_spec['hydrogen'] for y in old_spec['helium3']]
    # Refuse an unfinished generator before creating a candidate directory.
    for i in range(len(keys)):
        for root in [a.previous, a.addition]:
            path = root/f'plane-{i:03d}'/'source.json.gz'
            inputs[str(path.resolve())] = digest(path)
    spec = {k: old_spec[k] for k in ['hydrogen', 'helium3', 'logT',
                                   'probe_sha256', 'source_archive_sha256']}
    spec['logQ'] = q
    if options != [3, 1, -2]:
        spec.update(source_options=options, source_radiation_included=False)
    spec['density_merge_sources'] = {str(path.resolve()): digest(path) for path in specs}
    a.output.mkdir(parents=True)
    (a.output/'specification.json').write_text(json.dumps(spec, indent=2)+'\n')

    def merge(job):
        i, (x, y) = job
        paths = [root/f'plane-{i:03d}'/'source.json.gz' for root in [a.previous, a.addition]]
        old, extra = [json.loads(gzip.decompress(path.read_bytes())) for path in paths]
        m = mixture(x, y)
        for raw, axis in [(old, old_q), (extra, extra_q)]:
            if (any(raw[k] != v for k, v in m.items()) or raw['logT'] != spec['logT']
                    or raw['logQ'] != axis or raw['options'] != options
                    or raw['probe_sha256'] != spec['probe_sha256']
                    or raw['source_archive_sha256'] != spec['source_archive_sha256']
                    or raw['version'] != 'FreeEOS 3.0.0'
                    or len(raw['data']) != len(spec['logT'])*len(axis)
                    or any(len(row) != 22 or not all(map(math.isfinite, row)) for row in raw['data'])):
                raise ValueError('source plane identity, dimensions or responses differ')
        rows, overlap = [], []
        for it, t in enumerate(spec['logT']):
            retained = old['data'][it*len(old_q):(it+1)*len(old_q)]
            added = extra['data'][it*len(extra_q):(it+1)*len(extra_q)]
            left, right = retained[-1], added[0]
            converged = left[0] == right[0] == 0
            difference = (max(abs(v-w)/max(1., abs(v), abs(w))
                              for v, w in zip(left[1:], right[1:], strict=True)) if converged else None)
            if converged and difference > 1e-8:
                raise ValueError(f'independent density overlap differs at X={x}, X3={y}, logT={t}')
            overlap.append({'logT': t, 'source_flags': [left[0], right[0]],
                            'maximum_scaled_difference': difference})
            rows.extend(retained)
            rows.extend(added[1:])
        if not any(r['maximum_scaled_difference'] is not None for r in overlap):
            raise ValueError('no converged source overlap connects the density ranges')
        raw = {**old, 'logQ': q, 'data': rows,
               'density_merge_source_sha256': [inputs[str(path.resolve())] for path in paths]}
        folder = a.output/f'plane-{i:03d}'
        folder.mkdir()
        (folder/'mixture.json').write_text(json.dumps(m, indent=2)+'\n')
        output = folder/'source.json.gz'
        output.write_bytes(gzip.compress(encode(raw), mtime=0))
        record = {'hydrogen': x, 'helium3': y, 'retained_rows': len(old['data']),
                  'added_rows': len(spec['logT'])*(len(extra_q)-1),
                  'source_failures_for_masking': sum(row[0] != 0 for row in rows),
                  'overlap': overlap, 'source_sha256': digest(output)}
        (folder/'merge_receipt.json').write_text(json.dumps(record, indent=2)+'\n')
        print(f'merged X={x:.4g}, X3={y:.4g}', flush=True)
        return record

    with ThreadPoolExecutor(a.jobs) as pool:
        records = list(pool.map(merge, enumerate(keys)))
    if any(digest(Path(path)) != checksum for path, checksum in inputs.items()):
        raise ValueError('a source changed during assembly')
    (a.output/'merge_manifest.json').write_text(json.dumps(
        {'scope': __doc__, 'input_sha256': inputs, 'planes': records}, indent=2)+'\n')


if __name__ == '__main__':
    main()
