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


KEV_TO_K = 1e3*1.602176634e-12/1.380649e-16


def common_rectangles(sources):
    """Same original-cell rectangles as import_tops_mixtures.py."""
    tt, rr = sources[0][:2]
    if any(s[:2] != (tt, rr) for s in sources):
        raise ValueError('unaligned family source grids')
    result = {}
    for label, temperatures in [('low', tt), ('high', [t for t in tt if t >= .025])]:
        densities = []
        for rho in rr:
            if any((t, rho) in s[3] for s in sources for t in temperatures):
                break
            densities.append(rho)
        if len(densities) < 4:
            raise ValueError('insufficient family rectangle support')
        result[label] = (densities[0], densities[-1])
    return result


def active_top_weight(t_keV, rho, rectangles, gas_logr_range):
    """TOPS weight at a source node inside the current stellar blend domain.

    This does not bound interpolation between thermal/density source nodes,
    or establish that the star visits any particular supported source point.
    Joins match StellarMixtureOpacity: logT 4.4..4.5 and 5.6..5.7.
    """
    lt = math.log10(t_keV*KEV_TO_K)
    if lt <= 4.4:
        return None
    label = 'high' if lt >= 5.7 else 'low'
    rmin, rmax = rectangles[label]
    if not rmin <= rho <= rmax:
        return None
    if 5.6 < lt < 5.7 and not rectangles['high'][0] <= rho <= rectangles['high'][1]:
        return None
    if lt >= 4.5:
        return 1.
    if not gas_logr_range[0] <= math.log10(rho)-3*(lt-6) <= gas_logr_range[1]:
        return None
    u = (lt-4.4)/.1
    return u*u*(3-2*u)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['family_manifest', 'heldout_manifest', 'output']:
        p.add_argument(name, type=Path)
    p.add_argument('--summary-only', action='store_true', help='omit individual rows; retain domains, counts, extrema and provenance')
    p.add_argument('--gas-logr-range', type=float, nargs=2, default=[-8.,6.],
                   help='AESOPUS density-coordinate bounds in the 4.4..4.5 blend (default -8 6)')
    a = p.parse_args()
    family = json.loads(a.family_manifest.read_text())['planes']
    heldouts = json.loads(a.heldout_manifest.read_text())['planes']
    if not a.gas_logr_range[0] < a.gas_logr_range[1]:
        raise ValueError('invalid gas density support')
    cache = {}
    def source(manifest, record):
        path = (manifest.parent/record['file']).resolve()
        key = (path, json.dumps(record, sort_keys=True))
        if key not in cache:
            if sha(manifest.parent/record['request']) != record['request_sha256']:
                raise ValueError('changed source request')
            cache[key] = read(path, record)
        return cache[key]
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
        sources = [source(manifest, record) for manifest, record in records]
        rectangles = common_rectangles([source(a.family_manifest, record) for record in same_z])
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
        active = []
        for row in rows:
            weight_top = active_top_weight(row['T_keV'], row['density_g_cm3'], rectangles, a.gas_logr_range)
            if weight_top is not None:
                active.append({**row, 'TOPS_log_opacity_weight': weight_top,
                    'blended_relative_difference': math.expm1(weight_top*math.log1p(row['relative_difference']))})
        comparisons.append({'X': target['X'], 'Z': target['Z'], 'sources': [lo, hi, target],
            'original_cell_count': len(rows), 'excluded_cell_count': len(tt)*len(rr)-len(rows),
            'worst': max(rows, key=lambda r: abs(r['relative_difference'])),
            'hot_rectangle_worst': max(hot, key=lambda r: abs(r['relative_difference'])) if hot else None,
            'family_density_rectangles': rectangles, 'active_source_node_count': len(active),
            'active_source_nodes_worst': max(active, key=lambda r: abs(r['blended_relative_difference'])) if active else None,
            **({} if a.summary_only else {'rows': rows})})
    report = {'scope': __doc__, 'manifest_sha256': {str(path): sha(path)
              for path in [a.family_manifest, a.heldout_manifest]},
              'active_domain_scope': 'Source nodes within the original family rectangles and current logT joins 4.4..4.5, 5.6..5.7. Only the TOPS contribution changes; AESOPUS is held fixed. Not a stellar-profile error or a bound between source nodes.',
              'gas_logR_support': a.gas_logr_range,
              'audit_script_sha256': sha(Path(__file__)), 'comparisons': comparisons}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps([{k: v for k, v in c.items() if k not in ['rows', 'sources']}
                      for c in comparisons], indent=2))


if __name__ == '__main__':
    main()
