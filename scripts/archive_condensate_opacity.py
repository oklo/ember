#!/usr/bin/env python3
"""Validate and archive a completed offline depleted-opacity plane."""
import argparse
import gzip
import json
from pathlib import Path
from generate_nongrey_grid import composition, temperatures, sequence, input_fingerprint
from nongrey_opacity import validate_table
from prepare_nongrey_sources import digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('work',type=Path);p.add_argument('archive',type=Path)
    a=p.parse_args()
    provenance=json.loads((a.work/'provenance.json').read_text())
    spec=provenance['specification'];prepared=provenance['prepared']
    abundance,_=composition(provenance['XH'],provenance['X3'],spec['metals'])
    ts=temperatures(spec);rs=sequence(spec['log_density'])
    validate_table(a.work/'opacity.bin',abundance,ts,rs)
    files=[a.work/'provenance.json',a.work/'opacity.bin']
    for i,t in enumerate(ts):
        directory=a.work/f'temperature-{i:03d}'
        saved=json.loads((directory/'completed.json').read_text())
        if saved['input_sha256']!=input_fingerprint(prepared['executables']['synspec'],directory):
            raise ValueError('opacity source input identity changed')
        if not {'fort.63','fort.29','run.log'}.issubset(saved['outputs']):
            raise ValueError('incomplete source receipt')
        for name,expected in saved['outputs'].items():
            if digest(directory/name)!=expected:raise ValueError('opacity source output changed')
        validate_table(directory/'fort.63',abundance,[t],rs)
        names=['completed.json','condensate-abundances.dat','condensates.sha256',
               'ember-condensates.cfg','ember-masses.dat','fort.2','fort.5',
               'fort.55','physics.json','tas',*saved['outputs']]
        if (directory/'recovery.json').exists():names.append('recovery.json')
        if (directory/'reused.json').exists():names.append('reused.json')
        files += [directory/name for name in sorted(set(names))]
    a.archive.mkdir(parents=True,exist_ok=True)
    hashes={}
    for path in files:
        target=a.archive/path.relative_to(a.work)
        target=target.with_name(target.name+'.gz')
        target.parent.mkdir(parents=True,exist_ok=True)
        content=gzip.compress(path.read_bytes(),compresslevel=9,mtime=0)
        if target.exists() and target.read_bytes()!=content:
            raise ValueError('archive already contains different source bytes')
        target.write_bytes(content)
        hashes[str(target.relative_to(a.archive))]=digest(target)
    kind=('Original equilibrium-depleted gas opacity; zero grain opacity' if provenance['mode']=='equilibrium'
          else 'Original gas-only opacity control; no condensate depletion')
    record={'kind':kind+'; no atmosphere acceptance implied','mode':provenance['mode'],
            'XH':provenance['XH'],'X3':provenance['X3'],
            'opacity_sha256':digest(a.work/'opacity.bin'),'files_sha256':hashes}
    (a.archive/'archive.json').write_text(json.dumps(record,indent=2)+'\n')
    print(a.archive,len(hashes),'source artifacts')


if __name__=='__main__':main()
