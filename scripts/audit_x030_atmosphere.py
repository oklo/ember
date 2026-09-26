#!/usr/bin/env python3
"""Review the twelve X=.30 sources and exact retention in the private atmosphere extension."""
import hashlib
import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from assemble_nongrey_grid import load_continuation
from atmosphere_capacity_identity import CapacityEquivalence


def read(path):
    lines = path.read_text().splitlines(); at = lines.index('data')
    labels = ['hydrogen', 'helium3', 'log_teff', 'log_g']; axes = {}
    header = []
    for line in lines[:at]:
        key, *fields = line.split()
        if key in labels:
            assert int(fields[0]) == len(fields)-1
            axes[key] = [float(x) for x in fields[1:]]
        elif key != 'source': header.append(line)
    order = [axes[k] for k in labels]
    assert len(lines)-at-1 == math.prod(map(len, order))
    return header, axes, dict(zip(itertools.product(*order), lines[at+1:]))


def main():
    output = Path('docs/results/fable_x030_atmosphere_review_v1.json')
    if output.exists(): raise ValueError('preserve completed review')
    sealed = Path('docs/research/fable/results/fable-atmosphere-x030-assembly-v1')
    hashes = {}
    def pin(path):
        p = Path(path); h = hashlib.sha256()
        with p.open('rb') as f:
            for block in iter(lambda: f.read(1024*1024), b''): h.update(block)
        hashes[str(p)] = h.hexdigest(); return h.hexdigest()
    ready = json.loads((sealed/'READY.json').read_text())
    for name, expected in ready['files_sha256'].items(): assert pin(sealed/name) == expected
    manifest = json.loads((sealed/'ember-fable-nongrey-x030-candidate-v1.manifest.json').read_text())
    old = Path('data/atmosphere/nongrey_gs98_z020_exhaustion_t6400_g600_v1.dat')
    new = sealed/'ember-fable-nongrey-x030-candidate-v1.dat'
    assert pin(new) == manifest['table_sha256']; pin(old)
    oh, oa, before = read(old); nh, na, after = read(new)
    assert oh == nh and oa == na
    retained = [key for key, row in before.items() if row != '0']
    assert all(after[key] == before[key] for key in retained)
    added = {key for key in before if before[key] == '0' and after[key] != '0'}
    assert len(retained) == 300 and len(added) == 12
    expected = {(0.3, y, min(oa['log_teff'], key=lambda v: abs(v-math.log10(t))), g)
                for y, t, g in itertools.product((0., .12), (5200, 5400, 5600), (5.15, 5.4))}
    assert added == expected
    capacity = CapacityEquivalence(manifest['capacity_replays']); hashes.update(capacity.inputs)
    sources = []
    for entry in manifest['extension']:
        x, y, t, g = entry['coordinates']
        if x != .3 or t not in (5200, 5400, 5600): continue
        key, state, spec, record, inputs = load_continuation(entry['work'], manifest['physics'], capacity)
        assert list(key) == entry['coordinates']
        assert spec == entry['source_specification'] and state == entry['state']
        stored = after[(x, y, min(oa['log_teff'], key=lambda v: abs(v-math.log10(t))), g)].split()
        assert stored[0] == '1'
        assert float(stored[1]) == math.log10(state['T'])
        assert float(stored[2]) == math.log10(state['Pgas'])
        hashes.update(inputs)
        sources.append(dict(coordinates=key, state=state, work=entry['work']))
        print(key, 'raw source and stored values pass', flush=True)
    assert len(sources) == 12
    pin(__file__)
    output.write_text(json.dumps(dict(created_utc=datetime.now(timezone.utc).isoformat(),
        outcome='passed_raw_sources_and_retention', retained_unchanged=len(retained),
        added_sources=sources, sha256=hashes,
        scope='Source validity and exact assembly only; composition interpolation and stellar response remain separate.'), indent=2)+'\n')


if __name__ == '__main__': main()
