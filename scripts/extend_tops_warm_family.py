#!/usr/bin/env python3
"""Extend a warm opacity family with retained native low-hydrogen TOPS planes.

Preserve the existing axes and composition rows. Every added value must have
an original, checksum-verified source cell without density substitution.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import shlex

from import_tops_composition import read


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-manifest', type=Path, required=True)
    parser.add_argument('--warm-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('preserve existing output')
    family = json.loads(args.source_manifest.read_text())
    manifest = args.warm_manifest.read_text().splitlines()
    header = shlex.split(manifest[0])
    if header[:2] != ['EMBER_OPACITY_MIXTURE', '1'] or header[3] != 'logRho':
        raise ValueError('native-density opacity family required')
    if len(manifest) != int(header[2]) + 1:
        raise ValueError('incomplete family manifest')
    inputs = {str(p.resolve()): digest(p) for p in
              [args.source_manifest, args.warm_manifest, Path(__file__),
               Path(__file__).with_name('import_tops_composition.py')]}
    products = []
    unchanged = added = 0
    # Validate and assemble all source planes before creating the output.
    for entry in manifest[1:]:
        z, filename = shlex.split(entry)
        z = float(z)
        table = args.warm_manifest.parent / filename
        inputs[str(table.resolve())] = digest(table)
        lines = table.read_text().splitlines()
        nx, nt, nr = map(int, lines[0].split()[:3])
        density = list(map(float, lines[1].split()))
        temperature = list(map(float, lines[2].split()))
        if len(lines) != 3 + nx * (nt + 1) or len(density) != nr or len(temperature) != nt:
            raise ValueError('inconsistent existing table dimensions')
        xmin = float(lines[3].split()[0])
        missing = sorted((p for p in family['planes'] if p['Z'] == z and p['X'] < xmin),
                         key=lambda p: p['X'])
        if len({p['X'] for p in missing}) != len(missing):
            raise ValueError('duplicate source composition')
        body = []
        for plane in missing:
            source = args.source_manifest.parent / plane['file']
            tt, rr, values, excluded = read(source, plane)
            inputs[str(source.resolve())] = digest(source)
            if len(tt) != nt or len(rr) < nr:
                raise ValueError('source dimensions do not cover the warm table')
            if any(abs(math.log10(t * 1e3 * 1.602176634e-12 / 1.380649e-16) - q) > 1e-13
                   for t, q in zip(tt, temperature)):
                raise ValueError('source temperature coordinates differ')
            if any(abs(math.log10(r) - q) > 1e-13 for r, q in zip(rr[:nr], density)):
                raise ValueError('source density coordinates differ')
            if any((t, r) in excluded for t in tt for r in rr[:nr]):
                raise ValueError('added plane contains a substituted density')
            body.append(f'{plane["X"]:.17g} {z:.17g}')
            body.extend(' '.join(format(math.log10(values[t, r]), '.17g') for r in rr[:nr])
                        for t in tt)
            added += nt * nr
        contents = '\n'.join([
            f'{nx + len(missing)} {nt} {nr} Retained native TOPS ATOMIC; added low-hydrogen planes',
            lines[1], lines[2], *body, *lines[3:]]) + '\n'
        products.append((z, f'Z{z:.3f}.dat', contents))
        unchanged += nx * nt * nr
    if len({name for _, name, _ in products}) != len(products):
        raise ValueError("metallicities collide in output filenames")
    args.output.mkdir(parents=True)
    for _, filename, contents in products:
        (args.output / filename).write_text(contents)
    output = args.output / 'opacity_warm.dat'
    output.write_text(manifest[0] + '\n' + ''.join(
        f'{z:.17g} "{name}"\n' for z, name, _ in products))
    report = dict(input_sha256=inputs, unchanged_states=unchanged, added_states=added,
                  output_sha256={p.name: digest(p) for p in
                                 [output, *[args.output / name for _, name, _ in products]]},
                  scope='Retained native source cells; no extrapolation or density substitutions.')
    (args.output / 'assembly.json').write_text(json.dumps(report, indent=2) + '\n')
    print(f'{added} source values added; {unchanged} existing values retained')


if __name__ == '__main__':
    main()
