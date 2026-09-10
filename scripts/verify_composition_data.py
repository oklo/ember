#!/usr/bin/env python3
"""Verify checksums and reproduce every composition table from archived sources.

Run from the repository root. Python standard library only; no network or
Fortran needed. Independent EOS reference queries have a separate generator.
"""
import hashlib
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--extended',action='store_true')
    args=parser.parse_args()
    eos=Path('data/eos')
    spec=json.loads((eos/('sources/freeeos300_extended_manifest.json' if args.extended else 'sources/freeeos300_composition_manifest.json')).read_text())
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
        if args.extended:
            for script,source,target in [
                ('import_tops_mixtures.py',opacity/'sources/tops_gs98_mixture_manifest.json',opacity),
                ('import_aesopus_mixtures.py',opacity/'sources/aesopus21_gs98_mixtures.zip',opacity),
                ('import_conduction.py',Path('data/conduction/sources'),Path('data/conduction'))]:
                destination=output/script;destination.mkdir()
                subprocess.run([sys.executable,'scripts/'+script,str(source),str(destination)],check=True)
                for table in destination.iterdir():
                    if table.read_bytes()!=(target/table.name).read_bytes():
                        raise ValueError(f'non-reproducible extended input: {table.name}')
    print(f'All {len(spec["planes"])} EOS planes and selected opacity/conductivity inputs reproduce byte for byte.')


if __name__=='__main__':main()
