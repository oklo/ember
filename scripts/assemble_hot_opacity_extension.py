#!/usr/bin/env python3
"""Add original lower-hydrogen TOPS hot planes to an unchanged opacity family.

Cool TOPS and AESOPUS files are copied unchanged. The high-temperature source
axes and every old hydrogen plane must agree exactly. This assembly check
does not replace independent opacity, derivative or stellar-profile checks.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import shutil


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def planes(path):
    lines = path.read_text().splitlines()
    nx, nt, nr = map(int, lines[0].split()[:3])
    if len(lines) != 3 + nx*(nt+1):
        raise ValueError('invalid opacity table length')
    rows = [lines[3+i*(nt+1):3+(i+1)*(nt+1)] for i in range(nx)]
    xs = [float(row[0].split()[0]) for row in rows]
    if xs != sorted(set(xs)) or xs[0] < 0:
        raise ValueError('invalid hydrogen coordinates')
    return (nt, nr, lines[1:3]), xs, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['previous', 'imported', 'output']:
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('use a fresh output directory')
    replacements = {}; evidence = []
    for label in ['low', 'high']:
        name = f'tops_gs98_mixture_{label}.dat'
        if (args.previous/name).read_bytes() != (args.imported/name).read_bytes():
            raise ValueError('TOPS metallicity manifest changed')
    manifest = args.previous/'tops_gs98_mixture_high.dat'
    for line in manifest.read_text().splitlines()[1:]:
        z, name = shlex.split(line)
        if Path(name).name != name:
            raise ValueError('expected a table within its family directory')
        old, new = args.previous/name, args.imported/name
        old_axes, old_x, old_rows = planes(old)
        new_axes, new_x, new_rows = planes(new)
        added = len(new_x)-len(old_x)
        if added <= 0 or old_axes != new_axes or new_x[added:] != old_x or new_rows[added:] != old_rows:
            raise ValueError('extension must preserve every old source cell and add only lower hydrogen')
        replacements[name] = new
        evidence.append({'metallicity': float(z), 'added_hydrogen': new_x[:added],
                         'previous_sha256': digest(old), 'extended_sha256': digest(new)})
    shutil.copytree(args.previous, args.output)
    for name, source in replacements.items():
        shutil.copy2(source, args.output/name)
    report = {'scope': __doc__, 'previous_directory': str(args.previous.resolve()),
              'imported_directory': str(args.imported.resolve()), 'extensions': evidence,
              'copied_unchanged_sha256': {str(p.relative_to(args.previous)): digest(p)
                  for p in args.previous.rglob('*') if p.is_file() and p.name not in replacements},
              'script_sha256': digest(Path(__file__))}
    (args.output/'hot_extension_manifest.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
