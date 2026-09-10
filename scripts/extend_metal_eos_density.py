#!/usr/bin/env python3
"""Calculate only missing low-density FreeEOS states, retaining old raw rows.

The new family uses the same source executable, composition and temperature
axes. Each new isotherm starts a fresh source process. Old source rows are
copied exactly; they are never replaced by a second calculation. Import and
independent thermodynamic checks are required before selecting the result.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess

from metal_eos_composition import mixture


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(value):
    return (json.dumps(value, separators=(',', ':'), allow_nan=False)+'\n').encode()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['probe', 'previous_manifest', 'work']:
        parser.add_argument(name, type=Path)
    parser.add_argument('--minimum-log-q', type=float, default=-3.)
    parser.add_argument('--jobs', type=int, default=4)
    args = parser.parse_args()
    if not 1 <= args.jobs <= 8 or not math.isfinite(args.minimum_log_q):
        raise ValueError('invalid extension request')
    previous = json.loads(args.previous_manifest.read_text())
    probe_sha = digest(args.probe)
    if probe_sha != previous['probe_sha256']:
        raise ValueError('source executable differs from the original family')
    old_q = previous['logQ']; ts = previous['logT']
    step = (old_q[-1]-old_q[0])/(len(old_q)-1)
    count = round((old_q[0]-args.minimum_log_q)/step)
    if count <= 0 or abs(old_q[0]-count*step-args.minimum_log_q) > 1e-12:
        raise ValueError('lower density limit must extend the existing uniform grid')
    added_q = [old_q[0]-(count-i)*step for i in range(count)]
    spec = {key: previous[key] for key in
            ['hydrogen', 'helium3', 'logT', 'probe_sha256', 'source_archive_sha256']}
    spec.update(logQ=added_q+old_q,
                previous_manifest=str(args.previous_manifest.resolve()),
                previous_manifest_sha256=digest(args.previous_manifest),
                generator_sha256=digest(Path(__file__)))
    args.work.mkdir(parents=True, exist_ok=True)
    saved = args.work/'specification.json'
    if saved.exists() and json.loads(saved.read_text()) != spec:
        raise ValueError('changed inputs require a new work directory')
    saved.write_text(json.dumps(spec, indent=2)+'\n')
    records = {(r['hydrogen'], r['helium3']): r for r in previous['planes']}
    keys = [(x, y) for x in spec['hydrogen'] for y in spec['helium3']]
    if set(keys) != set(records) or len(keys) != len(previous['planes']):
        raise ValueError('incomplete or repeated source composition')

    def calculate(job):
        index, (x, y) = job
        folder = args.work/f'plane-{index:03d}'; folder.mkdir(exist_ok=True)
        record = records[x, y]
        source = args.previous_manifest.parent.parent/record['source']
        if digest(source) != record['source_sha256']:
            raise ValueError('old raw source checksum mismatch')
        old = json.loads(gzip.decompress(source.read_bytes()))
        m = mixture(x, y)
        if any(old[k] != v for k, v in m.items()):
            raise ValueError('composition recipe differs from the original source')
        if old['logT'] != ts or old['logQ'] != old_q or old['probe_sha256'] != probe_sha:
            raise ValueError('source coordinates or executable differ')
        if len(old['data']) != len(ts)*len(old_q):
            raise ValueError('incomplete old source rows')
        (folder/'mixture.json').write_text(json.dumps(m, indent=2)+'\n')
        rows = []; source_failures = 0
        for it, t in enumerate(ts):
            request = ' '.join(map(str, m['eps']))+'\n3 1 -2\n'+''.join(
                f'{math.log(m["source_mass_scale"])+math.log(10)*(q+1.5*(t-6)):.17g} {math.log(10)*t:.17g}\n'
                for q in added_q)
            fingerprint = hashlib.sha256((probe_sha+request).encode()).hexdigest()
            cache = folder/f'temperature-{it:03d}.json.gz'
            if cache.exists():
                entry = json.loads(gzip.decompress(cache.read_bytes()))
                if entry['input_sha256'] != fingerprint or entry['input'] != request:
                    raise ValueError('cached low-density request differs')
                added = entry['data']
            else:
                result = subprocess.run([str(args.probe.resolve())], input=request,
                                        text=True, capture_output=True, timeout=600)
                added = [list(map(float, line.split())) for line in result.stdout.splitlines()]
                if result.returncode:
                    (folder/f'temperature-{it:03d}.failure.log').write_text(result.stdout+'\n'+result.stderr)
                    raise ValueError(f'source failed at XH={x}, X3={y}, logT={t}')
                for row in added:
                    if len(row) != 22 or not all(math.isfinite(v) for v in row):
                        raise ValueError('invalid direct source response')
                    row[2] /= m['source_mass_scale']
                    for k in [5, 6, 9, 10, 11, 12, 13, 16, 17]:
                        row[k] *= m['source_mass_scale']
                entry = {'input_sha256': fingerprint, 'input': request,
                         'stderr': result.stderr, 'data': added}
                cache.write_bytes(gzip.compress(encode(entry), mtime=0))
            if len(added) != count or any(len(r) != 22 or not all(math.isfinite(v) for v in r) for r in added):
                raise ValueError('incomplete or nonfinite new source rows')
            source_failures += sum(r[0] != 0 for r in added)
            rows.extend(added)
            rows.extend(old['data'][it*len(old_q):(it+1)*len(old_q)])
        new = {**old, 'logQ': spec['logQ'], 'data': rows,
               'retained_source_sha256': record['source_sha256']}
        output = folder/'source.json.gz'
        output.write_bytes(gzip.compress(encode(new), mtime=0))
        evidence = {'hydrogen': x, 'helium3': y, 'source_sha256': digest(output),
                    'retained_rows': len(old['data']), 'new_rows': len(ts)*count,
                    'new_source_failures_for_masking': source_failures,
                    'retained_source_sha256': record['source_sha256']}
        (folder/'extension_receipt.json').write_text(json.dumps(evidence, indent=2)+'\n')
        print(json.dumps(evidence), flush=True)
        return evidence

    with ThreadPoolExecutor(args.jobs) as pool:
        results = list(pool.map(calculate, enumerate(keys)))
    if digest(args.probe) != probe_sha or digest(args.previous_manifest) != spec['previous_manifest_sha256']:
        raise ValueError('source inputs changed during generation')
    (args.work/'extension_manifest.json').write_text(json.dumps(
        {'scope': __doc__, 'specification': spec, 'planes': results}, indent=2)+'\n')


if __name__ == '__main__':
    main()
