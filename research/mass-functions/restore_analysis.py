"""Restore missing source scripts to their original workspace layout; never overwrite."""
from pathlib import Path
import argparse
import json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True, help='Existing Ember checkout')
    args = parser.parse_args()
    root = args.root.resolve()
    if not (root / 'CMakeLists.txt').is_file():
        parser.error('--root must contain Ember CMakeLists.txt')
    bundle = Path(__file__).resolve().parent
    workspaces = json.loads((bundle / 'provenance/workspaces.json').read_text())
    restored = skipped = 0
    for group, workspace in workspaces.items():
        for src in sorted((bundle / 'analysis' / group).rglob('*')):
            if not src.is_file() or src.suffix not in ('.py', '.md'):
                continue
            dst = root / 'out' / workspace / src.relative_to(bundle / 'analysis' / group)
            dst.parent.mkdir(parents=True, exist_ok=True)
            try:
                with dst.open('x') as f:
                    f.write(src.read_text().replace('@EMBER_ROOT@', str(root)))
                restored += 1
            except FileExistsError:
                skipped += 1
    print(f'Restored {restored} missing files; left {skipped} existing files unchanged.')
    print('Scientific input data and frozen executables must be supplied separately.')


if __name__ == '__main__':
    main()
