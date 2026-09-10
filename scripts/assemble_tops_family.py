#!/usr/bin/env python3
"""Assemble a distinct candidate TOPS family from pinned source manifests.

The output preserves provenance and original sources. Composition accuracy
must subsequently be measured against independently calculated heldouts.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

from import_tops_composition import read


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def assemble(manifests, output):
    if output.exists():raise FileExistsError('use a new candidate-family manifest')
    planes = {}
    for manifest in manifests:
        for record in json.loads(manifest.read_text())['planes']:
            if record.get('validation_role') == 'heldout':
                raise ValueError('heldout source needs an explicit, separately documented promotion to a candidate node')
            key = (record['Z'],record['X'])
            if key in planes:raise ValueError('duplicate candidate composition')
            source = (manifest.parent/record['file']).resolve(strict=True)
            request = (manifest.parent/record['request']).resolve(strict=True)
            if sha(request) != record['request_sha256']:raise ValueError('changed source request')
            read(source,record)
            planes[key] = {**record, 'file': os.path.relpath(source,output.parent.resolve()),
                           'request': os.path.relpath(request,output.parent.resolve()),
                           'parent_manifest': os.path.relpath(manifest.resolve(),output.parent.resolve())}
    if len({z for z,x in planes}) < 2:raise ValueError('at least two source metallicities required')
    result = {'scope': 'Candidate composition family only. Source identity validated; independent interpolation acceptance is separate.',
              'parent_manifests_sha256': {os.path.relpath(p.resolve(),output.parent.resolve()):sha(p) for p in manifests},
              'planes': [planes[k] for k in sorted(planes)]}
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x') as f:f.write(json.dumps(result,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('manifests',type=Path,nargs='+')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();assemble(a.manifests,a.output)


if __name__=='__main__':main()
