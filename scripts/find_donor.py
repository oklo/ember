#!/usr/bin/env python3
"""Choose the nearest accepted donor column for solve_column.py --initial-from.

usage: find_donor.py TARGET_COLUMN SEARCH_ROOT [SEARCH_ROOT ...]
Candidates are directories below the search roots whose result.json is accepted, with specification.json and
provenance.json. A candidate must have the target's source physics (assemble_nongrey_grid.physical_identity) and must
not share the target's own (Teff, log g, H): a donor is a neighbouring state, never a previous solve of the same state.
Distance: |dTeff|/100 K + |dlog g|/0.1 + |dH|/0.005. Prints JSON {donor, distance, teff_K, log_g, hydrogen}.
"""
import json, os, sys
from pathlib import Path
here = Path(__file__).resolve().parent
for p in [here, *here.parents]:
    if (p / 'scripts' / 'assemble_nongrey_grid.py').exists(): sys.path.insert(0, str(p / 'scripts')); break
    if (p / 'assemble_nongrey_grid.py').exists(): sys.path.insert(0, str(p)); break
from assemble_nongrey_grid import physical_identity

def coords(spec): return spec['teff_K'][0], spec['log_g'][0], spec['hydrogen'][0]

def main():
    target = Path(sys.argv[1]).resolve(); tspec = json.loads((target / 'specification.json').read_text())
    tid = physical_identity(tspec, json.loads((target / 'provenance.json').read_text())); tt, tg, tx = coords(tspec)
    best = None
    for root in sys.argv[2:]:
        for dirpath, dirnames, filenames in os.walk(root):
            dirnames[:] = [d for d in dirnames if d != 'attempts']
            if 'result.json' not in filenames or 'specification.json' not in filenames or 'provenance.json' not in filenames: continue
            d = Path(dirpath)
            try:
                if not json.loads((d / 'result.json').read_text()).get('accepted') or not (d / 'fort.7').exists(): continue
                spec = json.loads((d / 'specification.json').read_text())
                if len(spec['teff_K']) != 1 or physical_identity(spec, json.loads((d / 'provenance.json').read_text())) != tid: continue
            except Exception:
                continue
            t, g, x = coords(spec)
            if (t, g, x) == (tt, tg, tx) or d.resolve() == target: continue
            dist = abs(t - tt) / 100 + abs(g - tg) / 0.1 + abs(x - tx) / 0.005
            if best is None or dist < best[0]: best = (dist, str(d), t, g, x)
    if best is None: print(json.dumps({'donor': None})); return 1
    print(json.dumps(dict(donor=best[1], distance=best[0], teff_K=best[2], log_g=best[3], hydrogen=best[4])))

if __name__ == '__main__':
    sys.exit(main())
