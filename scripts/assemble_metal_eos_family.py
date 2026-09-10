#!/usr/bin/env python3
"""Assemble compatible imported EOS planes without replacing existing data.

The output has its own family manifest. Existing potential and source files
are reused only if byte-identical. Source physics and material axes must agree,
and every composition corner must be present.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import shutil


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('output_family', type=Path)
    p.add_argument('source_manifests', type=Path, nargs='+')
    a = p.parse_args()
    if a.output_family.exists():
        raise FileExistsError('use a new family manifest')
    common, planes, copies, parents = None, {}, {}, []
    for manifest in a.source_manifests:
        record = json.loads(manifest.read_text())
        identity = {k: record[k] for k in ['logT', 'logQ', 'probe_sha256', 'source_archive_sha256']}
        if common is None:
            common = identity
        elif identity != common:
            raise ValueError('source physics or material axes differ')
        parents.append({'manifest': str(manifest.resolve()), 'sha256': sha(manifest)})
        root = manifest.parent.parent
        if len(record['planes']) != len(record['hydrogen'])*len(record['helium3']):
            raise ValueError('incomplete input family')
        for row in record['planes']:
            key = row['hydrogen'], row['helium3']
            if key in planes:
                raise ValueError('overlapping composition planes; select disjoint input families')
            saved = dict(row)
            for field in ['potential', 'source']:
                source = root / row[field]
                if sha(source) != row[field+'_sha256']:
                    raise ValueError('source manifest checksum mismatch')
                name = Path(row[field]).name
                relative = Path('sources') / name if field == 'source' else Path(name)
                target = a.output_family.parent / relative
                if target in copies:
                    raise ValueError('duplicate output filename')
                if target.exists() and sha(target) != row[field+'_sha256']:
                    raise ValueError('refusing to replace different existing data')
                copies[target] = source
                saved[field] = str(relative)
            planes[key] = saved
    xs = sorted({x for x, y in planes})
    ys = sorted({y for x, y in planes})
    if len(xs) < 2 or len(ys) < 2 or set(planes) != set(itertools.product(xs, ys)):
        raise ValueError('incomplete output composition rectangle')
    provenance = a.output_family.parent / 'sources' / (a.output_family.stem+'_manifest.json')
    if provenance.exists():
        raise FileExistsError('output provenance already exists')
    lines = ['EMBER_METAL_HELMHOLTZ 1',
             'hydrogen '+str(len(xs))+' '+' '.join(map(str, xs)),
             'helium3 '+str(len(ys))+' '+' '.join(map(str, ys))]
    ordered = [planes[key] for key in itertools.product(xs, ys)]
    lines.extend(json.dumps(row['potential']) for row in ordered)
    for target, source in copies.items():
        if not target.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
    a.output_family.write_text('\n'.join(lines)+'\n')
    provenance.write_text(json.dumps({**common, 'hydrogen': xs, 'helium3': ys,
        'planes': ordered, 'family': a.output_family.name, 'family_sha256': sha(a.output_family),
        'parent_manifests': parents, 'assembly_script_sha256': sha(__file__)}, indent=2)+'\n')
    print(json.dumps({'family': str(a.output_family), 'planes': len(ordered),
                      'hydrogen': xs, 'helium3': ys, 'sha256': sha(a.output_family)}))


if __name__ == '__main__':
    main()
