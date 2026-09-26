#!/usr/bin/env python3
"""Restore pinned atmosphere data and line lists while retaining the selected executables.

Run under a controller with an elapsed-time and scratch-space bound. The
prepared receipt is written only after all restored data match their hashes.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

from prepare_nongrey_sources import ASSETS, LINES, SOURCES, data_digest, digest, fetch, run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prepared',type=Path)
    parser.add_argument('scratch',type=Path)
    args=parser.parse_args()
    assert not args.scratch.exists()
    root=args.scratch.resolve(); root.mkdir()
    prepared=json.loads(args.prepared.read_text())
    assert prepared['opacity_method']=='sampling'
    for name in ('tlusty','synspec'):
        assert digest(prepared[name])==prepared['executables'][name]
    for name,h in prepared['patches'].items():assert digest(SOURCES/name)==h
    for name in ['synple.tar.gz']+[n+'.gz' for n in LINES]:
        assert ASSETS[name][1]==prepared['inputs'][name]['sha256']
    downloads=root/'downloads';downloads.mkdir()
    archive=fetch(downloads,'synple.tar.gz')
    source=root/'synple';source.mkdir()
    with tarfile.open(archive) as tar:
        members=[m for m in tar if m.isfile() and
                 ('/data/' in m.name or m.name.endswith('/synspec/list2bin.f'))]
        assert members
        tar.extractall(source,members=members,filter='data')
    converter_source,=source.glob('*/synspec/list2bin.f')
    synple=converter_source.parents[1]
    assert data_digest(synple/'data')==prepared['data_sha256']
    converter=root/'list2bin'
    run(['gfortran','-O2','-std=legacy','-o',str(converter),str(converter_source)],
        root,root/'converter-build.log')
    line_root=root/'lines';line_root.mkdir()

    def restore(item):
        name,expected=item
        compressed=fetch(downloads,name+'.gz')
        directory=line_root/name;directory.mkdir()
        ascii_path=directory/'ascii'
        with gzip.open(compressed,'rb') as src,ascii_path.open('xb') as dst:
            shutil.copyfileobj(src,dst)
        with ascii_path.open('rb') as src:
            run([str(converter)],directory,directory/'conversion.log',src)
        assert 'lines included' in (directory/'conversion.log').read_text()
        binary=directory/'fort.12'
        assert digest(binary)==expected,name
        ascii_path.unlink() # Reproducible intermediate created by this job.
        print(json.dumps(dict(line_list=name,binary_bytes=binary.stat().st_size,sha256=expected)),flush=True)
        return str(binary)

    with ThreadPoolExecutor(max_workers=2) as pool:
        lines=list(pool.map(restore,zip(LINES,prepared['line_list_sha256'],strict=True)))
    prepared.update(synple=str(synple),line_lists=lines,
                    restoration=dict(parent_prepared=str(args.prepared.resolve()),
                        parent_prepared_sha256=digest(args.prepared),
                        restorer_sha256=digest(__file__),converter_sha256=digest(converter),
                        converter_source_sha256=digest(converter_source),
                        compiler=subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],
                        selected_executables_unchanged=True,
                        unused_dense_synspec_available=Path(prepared['dense_synspec']).is_file()))
    (root/'prepared.json').write_text(json.dumps(prepared,indent=2)+'\n')
    print(root/'prepared.json',flush=True)


if __name__=='__main__':
    main()
