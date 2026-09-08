#!/usr/bin/env python3
"""Verify checksums and reproduce every composition table from archived sources.

Run from the repository root. Python standard library only; no network or
Fortran needed. Independent EOS reference queries have a separate generator.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    eos=Path('data/eos')
    spec=json.loads((eos/'sources/freeeos300_composition_manifest.json').read_text())
    with tempfile.TemporaryDirectory(prefix='ember-composition-') as temporary:
        output=Path(temporary)
        for plane in spec['planes']:
            source=eos/'sources'/plane['source'];table=eos/plane['table']
            if sha(source)!=plane['source_sha256'] or sha(table)!=plane['table_sha256']:
                raise ValueError(f'checksum mismatch for X={plane["X"]}')
            target=output/table.name
            subprocess.run([sys.executable,'scripts/import_freeeos_potential.py',str(source),str(target),
                            '--hydrogen',str(plane['X'])],check=True)
            if target.read_bytes()!=table.read_bytes():raise ValueError(f'non-reproducible table: {table}')
        opacity=Path('data/opacity')
        subprocess.run([sys.executable,'scripts/import_tops_composition.py',
                        str(opacity/'sources/tops_composition_manifest.json'),str(output)],check=True)
        for label in ['low','high']:
            name=f'tops_gs98_composition_z020_{label}.dat'
            if (output/name).read_bytes()!=(opacity/name).read_bytes():
                raise ValueError(f'non-reproducible opacity: {name}')
    print('All four EOS planes and both multi-X opacity rectangles reproduce byte for byte.')


if __name__=='__main__':main()
