#!/usr/bin/env python3
"""Calculate hotter FreeEOS rows while retaining every existing source state.

The source executable, mixture recipes and density/composition coordinates
remain fixed. New rows require import, retained-potential checks and independent
thermodynamic comparisons before this family can be selected for evolution.
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
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def encode(value):
    return (json.dumps(value, separators=(',', ':'), allow_nan=False)+'\n').encode()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['probe', 'previous_manifest', 'work']:
        parser.add_argument(name, type=Path)
    parser.add_argument('--maximum-log-t', type=float, default=7.35)
    parser.add_argument('--jobs', type=int, default=4)
    args = parser.parse_args()
    if not 1 <= args.jobs <= 8 or not math.isfinite(args.maximum_log_t):
        raise ValueError('invalid extension request')
    previous = json.loads(args.previous_manifest.read_text())
    checksum = digest(args.probe)
    if checksum != previous['probe_sha256']:
        raise ValueError('source executable changed')
    ts, qs = previous['logT'], previous['logQ']
    step = ts[1]-ts[0]
    count = round((args.maximum_log_t-ts[-1])/step)
    if count <= 0 or abs(ts[-1]+count*step-args.maximum_log_t) > 1e-12:
        raise ValueError('maximum temperature must extend the uniform source grid')
    if any(abs((b-a)/step-1) > 1e-10 for a,b in zip(ts,ts[1:])):
        raise ValueError('nonuniform original temperature grid')
    added_ts = [ts[-1]+(i+1)*step for i in range(count)]
    spec = {k: previous[k] for k in ['hydrogen','helium3','logQ','probe_sha256','source_archive_sha256']}
    spec.update(logT=ts+added_ts, previous_manifest=str(args.previous_manifest.resolve()),
                previous_manifest_sha256=digest(args.previous_manifest), generator_sha256=digest(__file__))
    args.work.mkdir(parents=True, exist_ok=True)
    path = args.work/'specification.json'
    if path.exists() and json.loads(path.read_text()) != spec:
        raise ValueError('changed inputs require a new work directory')
    path.write_text(json.dumps(spec, indent=2)+'\n')
    records = {(r['hydrogen'],r['helium3']): r for r in previous['planes']}
    keys = [(x,y) for x in spec['hydrogen'] for y in spec['helium3']]
    if set(keys) != set(records) or len(keys) != len(previous['planes']):
        raise ValueError('incomplete or repeated source composition')

    def calculate(job):
        index, (x,y) = job
        folder = args.work/f'plane-{index:03d}'
        folder.mkdir(exist_ok=True)
        record = records[x,y]
        source = args.previous_manifest.parent.parent/record['source']
        if digest(source) != record['source_sha256']:
            raise ValueError('old source checksum mismatch')
        old = json.loads(gzip.decompress(source.read_bytes()))
        m = mixture(x,y)
        if any(old[k] != v for k,v in m.items()):
            raise ValueError('source composition recipe changed')
        if (old['logT'] != ts or old['logQ'] != qs or old['probe_sha256'] != checksum
                or old['source_archive_sha256'] != spec['source_archive_sha256']
                or old['version'] != 'FreeEOS 3.0.0' or old['options'] != [3,1,-2]
                or len(old['data']) != len(ts)*len(qs)):
            raise ValueError('old source coordinates, physics or row count differ')
        (folder/'mixture.json').write_text(json.dumps(m,indent=2)+'\n')
        retained = len(old['data'])
        rows = old['data'].copy()
        failures = 0
        for i,t in enumerate(added_ts):
            scale = m['source_mass_scale']
            request = ' '.join(map(str,m['eps']))+'\n3 1 -2\n'+''.join(
                f'{math.log(scale)+math.log(10)*(q+1.5*(t-6)):.17g} {math.log(10)*t:.17g}\n'
                for q in qs)
            fingerprint = hashlib.sha256((checksum+request).encode()).hexdigest()
            cache = folder/f'temperature-{len(ts)+i:03d}.json.gz'
            if cache.exists():
                entry = json.loads(gzip.decompress(cache.read_bytes()))
                if entry['input_sha256'] != fingerprint or entry['input'] != request:
                    raise ValueError('cached source request changed')
            else:
                result = subprocess.run([str(args.probe.resolve())], input=request, text=True,
                                        capture_output=True, timeout=600)
                added = [list(map(float,line.split())) for line in result.stdout.splitlines()]
                if (result.returncode or len(added) != len(qs)
                        or any(len(r) != 22 or not all(math.isfinite(v) for v in r) for r in added)):
                    (folder/f'temperature-{len(ts)+i:03d}.failure.log').write_text(result.stdout+'\n'+result.stderr)
                    raise ValueError(f'invalid hot source response at XH={x}, X3={y}, logT={t}')
                for row in added:
                    row[2] /= scale
                    for k in [5,6,9,10,11,12,13,16,17]:
                        row[k] *= scale
                entry = {'input_sha256':fingerprint, 'input':request, 'stderr':result.stderr, 'data':added}
                cache.write_bytes(gzip.compress(encode(entry),mtime=0))
            added = entry['data']
            if len(added) != len(qs) or any(len(r) != 22 or not all(math.isfinite(v) for v in r) for r in added):
                raise ValueError('incomplete hot source rows')
            failures += sum(r[0] != 0 for r in added)
            rows.extend(added)
        if rows[:retained] != old['data'] or digest(source) != record['source_sha256']:
            raise ValueError('old source states changed during extension')
        result = {**old, 'logT':spec['logT'], 'data':rows, 'retained_source_sha256':record['source_sha256']}
        output = folder/'source.json.gz'
        output.write_bytes(gzip.compress(encode(result),mtime=0))
        receipt = {'hydrogen':x, 'helium3':y, 'source_sha256':digest(output),
                   'retained_rows':retained, 'new_rows':len(added_ts)*len(qs),
                   'new_source_failures_for_masking':failures, 'retained_source_sha256':record['source_sha256']}
        (folder/'extension_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        print(json.dumps(receipt),flush=True)
        return receipt

    with ThreadPoolExecutor(args.jobs) as pool:
        results = list(pool.map(calculate,enumerate(keys)))
    if digest(args.probe) != checksum or digest(args.previous_manifest) != spec['previous_manifest_sha256']:
        raise ValueError('source inputs changed during extension')
    (args.work/'extension_manifest.json').write_text(json.dumps(
        {'scope':__doc__, 'specification':spec, 'planes':results},indent=2)+'\n')


if __name__ == '__main__':
    main()
