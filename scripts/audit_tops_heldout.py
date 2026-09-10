#!/usr/bin/env python3
"""Compare an independent TOPS composition with bracketing source planes.

This measures composition interpolation at identical source T/rho coordinates,
excluding every substituted cell. It is not a physical opacity uncertainty.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

from import_tops_composition import read


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['family_manifest', 'heldout_manifest', 'output']:
        p.add_argument(name, type=Path)
    a = p.parse_args()
    family = json.loads(a.family_manifest.read_text())['planes']
    heldouts = json.loads(a.heldout_manifest.read_text())['planes']
    comparisons = []
    for target in heldouts:
        same_z = sorted((v for v in family if v['Z'] == target['Z']), key=lambda v: v['X'])
        if any(v['X'] == target['X'] for v in same_z):
            raise ValueError('heldout is a runtime node')
        lower = [v for v in same_z if v['X'] < target['X']]
        upper = [v for v in same_z if v['X'] > target['X']]
        if not lower or not upper:
            raise ValueError('heldout is not bracketed')
        lo, hi = lower[-1], upper[0]
        records = [(a.family_manifest, lo), (a.family_manifest, hi),
                   (a.heldout_manifest, target)]
        sources = []
        for manifest, record in records:
            if sha(manifest.parent / record['request']) != record['request_sha256']:
                raise ValueError('changed source request')
            sources.append(read(manifest.parent / record['file'], record))
        tt, rr = sources[0][:2]
        if any(t != tt or r != rr for t, r, _, _ in sources):
            raise ValueError('unaligned source grids')
        weight = (target['X'] - lo['X']) / (hi['X'] - lo['X'])
        rows = []
        for t in tt:
            for rho in rr:
                if any((t, rho) in s[3] for s in sources):
                    continue
                left, right, truth = [s[2][t, rho] for s in sources]
                predicted = math.exp((1-weight)*math.log(left) + weight*math.log(right))
                rows.append({'T_keV': t, 'density_g_cm3': rho, 'heldout': truth,
                             'interpolated': predicted, 'relative_difference': predicted/truth-1})
        if not rows:
            raise ValueError('no common original cells')
        hot = [r for r in rows if r['T_keV'] >= .025]
        comparisons.append({'X': target['X'], 'Z': target['Z'], 'sources': [lo, hi, target],
            'original_cell_count': len(rows), 'excluded_cell_count': len(tt)*len(rr)-len(rows),
            'worst': max(rows, key=lambda r: abs(r['relative_difference'])),
            'hot_rectangle_worst': max(hot, key=lambda r: abs(r['relative_difference'])) if hot else None,
            'rows': rows})
    report = {'scope': __doc__, 'manifest_sha256': {str(path): sha(path)
              for path in [a.family_manifest, a.heldout_manifest]},
              'audit_script_sha256': sha(Path(__file__)), 'comparisons': comparisons}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps([{k: v for k, v in c.items() if k not in ['rows', 'sources']}
                      for c in comparisons], indent=2))


if __name__ == '__main__':
    main()
