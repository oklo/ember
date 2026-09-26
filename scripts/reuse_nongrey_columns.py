#!/usr/bin/env python3
"""Populate a NEW atmosphere work directory with verified identical columns.

This reuses completed calculations at unchanged physical coordinates and source
resolution. It does not interpolate, transform compositions, or run a solver.
The normal generator can fill the remaining columns with --reuse-opacity.
"""
import argparse
from datetime import datetime, timezone
import gzip
import itertools
import json
from pathlib import Path
import shutil

from assemble_helium_fraction_atmospheres import NUMERICS, load_plane
from assemble_nongrey_grid import physical_identity
from generate_nongrey_grid import input_fingerprint, temperatures, sequence
from prepare_nongrey_sources import digest


def reuse(prepared_path, specification_path, sources, output):
    prepared_path, specification_path = map(Path, [prepared_path, specification_path])
    prepared = json.loads(prepared_path.read_text())
    spec = json.loads(specification_path.read_text())
    output = Path(output)
    if output.exists():
        raise FileExistsError('source reuse requires a new directory')
    if len(spec['hydrogen']) != 1 or len(spec['helium3']) != 1:
        raise ValueError('one target composition plane required')
    h, he3 = spec['hydrogen'][0], spec['helium3'][0]
    identity = dict(physical_identity(spec, prepared),
                    numerical_settings={k: spec.get(k) for k in NUMERICS})
    selected, pins, opacity = {}, {}, None
    for source in sources:
        root = Path(source).resolve(strict=True)
        _, rows, own_identity, own_pins, own_spec = load_plane(
            root, h, he3/(1-h-sum(spec['metals'])))
        if (own_identity != identity or own_spec['approximation'] != spec['approximation']
                or temperatures(own_spec) != temperatures(spec)
                or sequence(own_spec['log_density']) != sequence(spec['log_density'])):
            raise ValueError('source reuse would change physics or numerical resolution')
        own_opacity = root/'plane-000/opacity/fort.63'
        if opacity is not None and digest(opacity) != digest(own_opacity):
            raise ValueError('identical opacity bytes required for source reuse')
        opacity = own_opacity
        pins.update(own_pins)
        for row in rows:
            t, g = row['coordinates'][2:]
            if t in spec['teff_K'] and g in spec['log_g']:
                if (t, g) in selected:
                    raise ValueError('duplicate accepted coordinate: choose one source explicitly')
                selected[t, g] = root, row
    if not selected:
        raise ValueError('no identical accepted source coordinates')
    if any(digest(p) != hsh for p, hsh in pins.items()):
        raise ValueError('source changed during reuse validation')
    output.mkdir(parents=True)
    (output/'specification.json').write_text(json.dumps(spec, indent=2)+'\n')
    (output/'provenance.json').write_text(json.dumps(prepared, indent=2)+'\n')
    target_opacity = output/'plane-000/opacity/fort.63'
    target_opacity.parent.mkdir(parents=True)
    # Use independent bytes: the generator may later rewrite its merged table.
    # A hard link here would put the original source archive at risk.
    shutil.copy2(opacity, target_opacity)
    reused = []
    for (t, g), (root, row) in selected.items():
        old = Path(row['validation']).parent
        directory = output/'plane-000'/f"model-{spec['teff_K'].index(t):03d}-{spec['log_g'].index(g):03d}"
        directory.mkdir()
        record = dict(row['source_record'])
        receipt = json.loads((old/'completed.json').read_text())
        names = set(receipt['outputs']) | {'fort.5', 'fort.15', 'tas', 'ember-masses.dat',
                                         'physics.json', 'opacity.sha256'}
        if 'initial_structure' in record:
            names.add('fort.8')
        for name in names:
            src = old/name
            if src.is_file():
                shutil.copy2(src, directory/name)
            else:
                (directory/name).write_bytes(gzip.decompress(src.with_name(name+'.gz').read_bytes()))
        shutil.copy2(old/'completed.json', directory/'completed.json')
        for kind in ['log', 'convergence', 'atmosphere_input', 'element_masses', 'parameters', 'initial_structure']:
            if kind not in record:
                continue
            name = Path(record[kind]).name
            shutil.copy2(root/record[kind], directory/name)
            record[kind] = str((directory/name).relative_to(output))
        if input_fingerprint(prepared['executables']['tlusty'], directory) != receipt['input_sha256']:
            raise ValueError('copied completed input fingerprint changed')
        if any(digest(directory/n) != hsh for n, hsh in receipt['outputs'].items()):
            raise ValueError('copied completed output hash changed')
        record['reused_source'] = str(old)
        (directory/'validated.json').write_text(json.dumps(record, indent=2)+'\n')
        reused.append(dict(coordinates=row['coordinates'], original=str(old), copy=str(directory)))
    # Independently reparse the new directory through the same importer that
    # will later assemble atmosphere tables. No newly solved source is claimed.
    _, copied, _, _, _ = load_plane(output, h, he3/(1-h-sum(spec['metals'])))
    if len(copied) != len(selected):
        raise ValueError('reused source count differs after independent import')
    for row in copied:
        if row['state'] != selected[tuple(row['coordinates'][2:])][1]['state']:
            raise ValueError('repackaged matching state changed')
    if any(digest(p) != hsh for p, hsh in pins.items()):
        raise ValueError('original source changed during repackaging')
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  reused_columns=reused, missing_coordinates=[list(k) for k in itertools.product(
                      spec['teff_K'], spec['log_g']) if k not in selected],
                  solver_calls=0, matching_states_bit_identical=True,
                  input_sha256={**pins, str(prepared_path): digest(prepared_path),
                      str(specification_path): digest(specification_path), __file__: digest(__file__)})
    (output/'reuse_manifest.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('prepared'); p.add_argument('specification'); p.add_argument('output')
    p.add_argument('--source', action='append', required=True)
    a = p.parse_args()
    r = reuse(a.prepared, a.specification, a.source, a.output)
    print(json.dumps({'reused_columns': len(r['reused_columns']),
                      'missing_columns': len(r['missing_coordinates']), 'solver_calls': 0}))
