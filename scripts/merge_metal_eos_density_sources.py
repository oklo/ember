#!/usr/bin/env python3
"""Join EOS density ranges at one shared density for each temperature.

Existing source rows are retained exactly. The shared coordinate checks the
effect of starting a new source process at high density; each mixture's
converged responses must agree before its merged source is written. This
prepares raw material data, not an accepted thermodynamic table or a stellar restart.
The addition may cover a contiguous upper-temperature subset. Unrequested
cold dense states are then explicitly null under a checked coverage declaration.
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


def overlap_difference(left, right):
    if left[0] != 0 or right[0] != 0:
        return None
    # Source iteration count (column 1) is not a physical response.
    return max(abs(v-w)/max(1., abs(v), abs(w))
               for v, w in zip(left[2:], right[2:], strict=True))


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
    for key in ['hydrogen', 'helium3', 'probe_sha256', 'source_archive_sha256']:
        if old_spec[key] != extra_spec[key]:
            raise ValueError('source composition, temperature or executable differs')
    old_t, extra_t = old_spec['logT'], extra_spec['logT']
    if (len(extra_t) < 5 or extra_t[0] not in old_t
            or old_t[old_t.index(extra_t[0]):] != extra_t
            or old_spec.get('source_coverage') or extra_spec.get('source_coverage')):
        raise ValueError('addition must use a contiguous upper-temperature subset of complete source cubes')
    if (extra_spec.get('source_options', [3, 1, -2]) != options
            or old_spec.get('source_radiation_included', True) != (options == [3, 1, -2])
            or extra_spec.get('source_radiation_included', True) != (options == [3, 1, -2])):
        raise ValueError('source physics or radiation treatment differs')
    precisions = [s.get('precision_fallback') for s in [old_spec, extra_spec]]
    if all(precisions) and precisions[0] != precisions[1]:
        raise ValueError('source precision fallbacks differ')
    precision = next((value for value in precisions if value), None)
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
    if old_t != extra_t:
        spec['source_coverage'] = dict(kind='hot_density_extension',
                                      original_logQ_max=old_q[-1], minimum_added_logT=extra_t[0])
    if options != [3, 1, -2]:
        spec.update(source_options=options, source_radiation_included=False)
    if precision:
        spec['precision_fallback'] = precision
    spec['density_merge_sources'] = {str(path.resolve()): digest(path) for path in specs}
    a.output.mkdir(parents=True)
    (a.output/'specification.json').write_text(json.dumps(spec, indent=2)+'\n')

    def merge(job):
        i, (x, y) = job
        paths = [root/f'plane-{i:03d}'/'source.json.gz' for root in [a.previous, a.addition]]
        old, extra = [json.loads(gzip.decompress(path.read_bytes())) for path in paths]
        m = mixture(x, y)
        for raw, axis, thermal, source_spec in [(old, old_q, old_t, old_spec),
                                               (extra, extra_q, extra_t, extra_spec)]:
            if (any(raw[k] != v for k, v in m.items()) or raw['logT'] != thermal
                    or raw['logQ'] != axis or raw['options'] != options
                    or raw['probe_sha256'] != spec['probe_sha256']
                    or raw['source_archive_sha256'] != spec['source_archive_sha256']
                    or raw['version'] != 'FreeEOS 3.0.0'
                    or len(raw['data']) != len(thermal)*len(axis)
                    or any(len(row) != 22 or not all(map(math.isfinite, row)) for row in raw['data'])):
                raise ValueError('source plane identity, dimensions or responses differ')
            if raw.get('precision_fallback') != source_spec.get('precision_fallback'):
                raise ValueError('source plane precision differs from its specification')
        rows, overlap = [], []
        for it, t in enumerate(spec['logT']):
            retained = old['data'][it*len(old_q):(it+1)*len(old_q)]
            if t not in extra_t:
                rows.extend(retained)
                rows.extend([None]*(len(extra_q)-1))
                continue
            jt = extra_t.index(t)
            added = extra['data'][jt*len(extra_q):(jt+1)*len(extra_q)]
            left, right = retained[-1], added[0]
            difference = overlap_difference(left, right)
            if difference is not None and difference > 1e-8:
                raise ValueError(f'independent density overlap differs at X={x}, X3={y}, logT={t}')
            overlap.append({'logT': t, 'source_flags': [left[0], right[0]],
                            'source_iteration_counts': [left[1], right[1]],
                            'maximum_scaled_difference': difference})
            rows.extend(retained)
            rows.extend(added[1:])
        if not any(r['maximum_scaled_difference'] is not None for r in overlap):
            raise ValueError('no converged source overlap connects the density ranges')
        raw = {**old, 'logQ': q, 'data': rows,
               'density_merge_source_sha256': [inputs[str(path.resolve())] for path in paths]}
        if 'source_coverage' in spec:
            raw['source_coverage'] = spec['source_coverage']
        if precision:
            raw['precision_fallback'] = precision
            overrides = []
            for source, axis, thermal in [(old, old_q, old_t), (extra, extra_q, extra_t)]:
                for override in source.get('precision_fallback_isotherms', []):
                    it = override['temperature_index']
                    if (not 0 <= it < len(thermal) or override['logT'] != thermal[it]
                            or not source.get('precision_fallback')):
                        raise ValueError('invalid source precision isotherm')
                    # Preserve the input fingerprint of the complete source
                    # isotherm and identify exactly which rows were retained.
                    overrides.append({**override, 'temperature_index': old_t.index(thermal[it]),
                                      'source_temperature_index': it,
                                      'source_logQ_range': [axis[0], axis[-1]],
                                      'retained_logQ_range': [axis[1] if source is extra else axis[0], axis[-1]]})
            raw['precision_fallback_isotherms'] = overrides
        folder = a.output/f'plane-{i:03d}'
        folder.mkdir()
        (folder/'mixture.json').write_text(json.dumps(m, indent=2)+'\n')
        output = folder/'source.json.gz'
        output.write_bytes(gzip.compress(encode(raw), mtime=0))
        record = {'hydrogen': x, 'helium3': y, 'retained_rows': len(old['data']),
                  'added_rows': len(extra_t)*(len(extra_q)-1),
                  'absent_source_rows': sum(row is None for row in rows),
                  'source_failures_for_masking': sum(row is not None and row[0] != 0 for row in rows),
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
