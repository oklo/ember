#!/usr/bin/env python3
"""Merge complete requested grey batches into one explicitly uncut source plane.

This is a source assembly, not a stellar opacity: refraction, family
interpolation and the final temperature joins are separate calculations.
Missing batches and inconsistent overlap values or exclusions are rejected.
"""
import argparse
import json
import math
from pathlib import Path

from audit_tops_electron_dispersion import KEV, KB
from fetch_tops_composition import digest
from import_tops_composition import read
from reduce_tops_group_factors import add_inputs, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('plan', type=Path)
    parser.add_argument('X', type=float)
    parser.add_argument('Z', type=float)
    parser.add_argument('output', type=Path)
    parser.add_argument('report', type=Path)
    parser.add_argument('--base-only', action='store_true', help='Require the complete original density range through 10000 g/cm3')
    args = parser.parse_args()
    if args.output.exists() or args.report.exists():
        raise FileExistsError('use new assembly outputs')
    plan = json.loads(args.plan.read_text())
    inputs = {str(p.resolve()): digest(p) for p in
              (args.plan, Path(__file__), Path('scripts/import_tops_composition.py'),
               Path('scripts/reduce_tops_group_factors.py'),
               Path('scripts/audit_tops_electron_dispersion.py'),
               Path('scripts/fetch_tops_composition.py'))}
    jobs = list(plan['requests'])
    for retained in plan.get('retained_requests', []):
        path = Path(retained['plan'])
        if digest(path) != retained['plan_sha256']:
            raise ValueError('retained source plan changed')
        add_inputs(inputs, {str(path.resolve()): digest(path)})
        jobs.append(json.loads(path.read_text())['requests'][retained['request_index']])
    selected = [r for r in jobs if (r['X'], r['Z']) == (args.X, args.Z) and
                (not args.base_only or max(r['densities_g_cm3']) <= 10000)]
    if not selected:
        raise ValueError('no requested source batches in the stated scope')
    cells, exclusions, records = {}, {}, []
    overlaps = 0
    for job in selected:
        root = Path(job['work'])
        receipt = json.loads((root / 'receipt.json').read_text())
        request = json.loads((root / 'request.json').read_text())
        recipe = json.loads((root / 'recipe.json').read_text())
        if ((receipt['X'], receipt['Z']) != (args.X, args.Z) or
                digest(root / 'request.json') != receipt['request_sha256'] or
                request['datype'] != 'gray' or request['plasnu'] != 'off' or request['lib'] != 'new' or
                recipe['job'] != job):
            raise ValueError('source identity or scientific request differs')
        dimensions = (len(job['temperatures_keV']), len(job['densities_g_cm3']))
        tt, rr, part, excluded = read(root / 'source.txt', receipt, dimensions=dimensions)
        if (tt != job['temperatures_keV'] or
                any(abs(x / y - 1) > 5e-5 for x, y in zip(rr, job['densities_g_cm3'], strict=True))):
            raise ValueError('returned axes differ from source plan')
        for key, value in part.items():
            is_excluded = key in excluded
            if key in cells:
                if exclusions[key] != is_excluded or (not is_excluded and cells[key] != value):
                    raise ValueError('inconsistent overlap value or source exclusion')
                overlaps += 1
            else:
                cells[key], exclusions[key] = value, is_excluded
        for path in root.iterdir():
            if path.is_file():
                add_inputs(inputs, {str(path.resolve()): digest(path)})
        records.append({'work': str(root), 'states': len(part), 'excluded_states': len(excluded)})
    tt, rr = sorted({t for t, r in cells}), sorted({r for t, r in cells})
    if args.base_only:
        manifest_path = Path(plan['baseline_manifest'])
        if digest(manifest_path) != plan['baseline_manifest_sha256']:
            raise ValueError('baseline manifest changed')
        add_inputs(inputs, {str(manifest_path.resolve()): digest(manifest_path)})
        original = next(r for r in json.loads(manifest_path.read_text())['planes']
                        if (r['X'], r['Z']) == (args.X, args.Z))
        source = manifest_path.parent / original['file']
        native_t, native_r, _, _ = read(source, original)
        add_inputs(inputs, {str(source.resolve()): digest(source)})
        if tt != native_t or rr != native_r:
            raise ValueError('incomplete original temperature-density axes')
    if len(tt) < 4 or len(rr) < 4:
        raise ValueError('insufficient table axes')
    if args.base_only and set(cells) != {(t, r) for t in tt for r in rr}:
        raise ValueError('incomplete original source rectangle')
    lines = [f'EMBER_OPACITY_TABLE 2 1 {len(tt)} {len(rr)} TOPS native UNCUT grey; refraction not applied',
             ' '.join(format(math.log10(r), '.17g') for r in rr),
             ' '.join(format(math.log10(t * KEV / KB), '.17g') for t in tt),
             f'{args.X:.17g} {args.Z:.17g}']
    prefixes = []
    for t in tt:
        values = []
        missing = False
        for rho in rr:
            key = t, rho
            if key not in cells or exclusions[key]:
                missing = True
            elif missing:
                raise ValueError('valid source states do not form a density prefix')
            else:
                values.append(cells[key])
        if len(values) < 4:
            raise ValueError('insufficient supported densities')
        prefixes.append(len(values))
        lines.append(str(len(values)) + ' ' + ' '.join(format(math.log10(v), '.17g') for v in values))
    verify(inputs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text('\n'.join(lines) + '\n')
    result = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
              'X': args.X, 'Z': args.Z, 'base_only': args.base_only,
              'source_batches': len(selected), 'source_states': len(cells),
              'source_exclusions': sum(exclusions.values()), 'exact_overlap_states': overlaps,
              'temperatures_keV': tt, 'densities_atomic_g_cm3': rr,
              'density_prefix_sizes': prefixes, 'sources': records,
              'input_sha256': inputs, 'output_sha256': {str(args.output.resolve()): digest(args.output)}}
    args.report.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in
                     ('scope', 'temperatures_keV', 'densities_atomic_g_cm3',
                      'density_prefix_sizes', 'sources', 'input_sha256')}), flush=True)


if __name__ == '__main__':
    main()
