#!/usr/bin/env python3
"""Measure dense opacity interpolation errors after including conduction.

Reuses the pinned opacity source and runtime queries, then evaluates the
selected conductive transport at those same material states. Removing all
radiation gives an additional fixed-state diagnostic, not an opacity accuracy
bound or an evolutionary error bound. No stellar input is accepted here.
"""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess

from import_tops_composition import read


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['source_audit', 'conduction_probe', 'conduction_table', 'output']:
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    source = json.loads(args.source_audit.read_text())
    rawpath = args.source_audit.with_suffix('.probe.json.gz')
    if digest(rawpath) != source['raw_probe_sha256']:
        raise ValueError('saved opacity queries changed')
    inputs = {**source['input_sha256'], **{str(p.resolve()): digest(p) for p in
              [args.source_audit, rawpath, args.conduction_probe, args.conduction_table, Path(__file__)]}}
    if any(digest(p) != checksum for p, checksum in inputs.items()):
        raise ValueError('source or runtime dependency changed')
    raw = json.loads(gzip.decompress(rawpath.read_bytes()))
    manifests = [Path(p) for p in source['input_sha256'] if p.endswith('/manifest.json')]
    if len(manifests) != 2:
        raise ValueError('expected source-node and independent-source manifests')
    records = []; probe_records = []; batch = 0
    for kind, path in zip(['nodes', 'heldouts'], manifests, strict=True):
        for entry in json.loads(path.read_text())['requests']:
            work = Path(entry['work']); receipt = json.loads((work/'receipt.json').read_text())
            temperatures, densities, cells, excluded = read(work/'source.txt', receipt,
                                                           dimensions=tuple(receipt['dimensions']))
            material = [(t, rho) for t in temperatures
                        if math.log10(t*1e3*1.602176634e-12/1.380649e-16) >= 5.7 for rho in densities]
            previous = raw[batch]; batch += 1
            rows = [json.loads(line) for line in previous['stdout'].splitlines()]
            if previous['Z'] != entry['Z'] or len(rows) != len(material):
                raise ValueError('opacity source ordering differs')
            states = []
            for (t, rho), row in zip(material, rows, strict=True):
                expected = [entry['X'], t*1e3*1.602176634e-12/1.380649e-16, rho]
                if row['query'] != expected:
                    raise ValueError('opacity source coordinate differs')
                if row['covered']:
                    if (t, rho) in excluded:
                        raise ValueError('opacity query admits substituted source')
                    states.append((row, cells[t, rho]))
            for y3 in [0., .12]:
                request = ''.join(' '.join(format(v, '.17g') for v in
                                  [entry['X'], entry['Z'], y3, row['query'][1], row['query'][2]])+'\n'
                                  for row, _ in states)
                response = subprocess.run([str(args.conduction_probe), str(args.conduction_table)],
                                          input=request, text=True, capture_output=True, check=True)
                conductivity = list(map(float, response.stdout.splitlines()))
                if len(conductivity) != len(states) or any(not math.isfinite(k) or k <= 0 for k in conductivity):
                    raise ValueError('invalid conduction response')
                probe_records.append(dict(kind=kind, X=entry['X'], Z=entry['Z'], Y3=y3,
                                          input=request, stdout=response.stdout))
                for (row, direct), kc in zip(states, conductivity, strict=True):
                    interpolation = row['kappa']
                    combined = 1/(1/interpolation+1/kc)
                    reference = 1/(1/direct+1/kc)
                    records.append(dict(kind=kind, X=entry['X'], Z=entry['Z'], Y3=y3,
                                        temperature_K=row['query'][1], density=row['query'][2],
                                        radiative_source=direct, radiative_interpolated=interpolation,
                                        conductive_opacity=kc,
                                        radiative_relative_difference=interpolation/direct-1,
                                        combined_relative_difference=combined/reference-1,
                                        remove_all_radiation_relative_increase=kc/combined-1))
    if batch != len(source['comparisons']):
        raise ValueError('source comparison count differs')
    def summary(selected):
        return dict(states=len(selected), maxima={key: max(selected, key=lambda r: abs(r[key])) for key in
                    ['radiative_relative_difference', 'combined_relative_difference',
                     'remove_all_radiation_relative_increase']})
    result = dict(scope=__doc__, input_sha256=inputs, accepted_for_evolution=False,
                  all_states=summary(records), independent=summary([r for r in records if r['kind'] == 'heldouts']),
                  hydrogen_burning_temperatures=summary([r for r in records if 1e6 <= r['temperature_K'] <= 2e7]),
                  composition_note='GS98 baryonic X/Z with Y3=0 and 0.12. The opacity source treats helium as He4; its existing isotope approximation is retained.')
    archive = args.output.with_suffix('.queries.json.gz')
    archive.write_bytes(gzip.compress(json.dumps(dict(queries=probe_records, records=records),
                                               separators=(',', ':'), allow_nan=False).encode(), mtime=0))
    result.update(raw_queries=str(archive), raw_queries_sha256=digest(archive))
    if any(digest(p) != checksum for p, checksum in inputs.items()):
        raise ValueError('an input changed during the diagnostic')
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: {name: r[name] for name, r in result[k]['maxima'].items()}
                      for k in ['all_states', 'independent', 'hydrogen_burning_temperatures']}))


if __name__ == '__main__':
    main()
