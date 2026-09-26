#!/usr/bin/env python3
"""Reuse the MESA comparison and exact TOPS cutoff-on/off state matches."""
import hashlib
import json
from pathlib import Path
import statistics

from import_tops_composition import read

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'docs/results/mesa_oplib_mean_comparison_v2.json'
OUTPUT = ROOT / 'docs/results/mesa_oplib_cutoff_comparison_v1.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def summarize(rows):
    keys = ['mesa_relative_to_uncut', 'mesa_relative_to_cutoff', 'cutoff_relative_to_uncut']
    return {key:dict(median=statistics.median(row[key] for row in rows),
                     minimum=min(row[key] for row in rows),
                     maximum=max(row[key] for row in rows),
                     worst=max(rows, key=lambda row:abs(row[key]))) for key in keys}


def main():
    assert not OUTPUT.exists()
    assert digest(INPUT) == '90e28efe30a22258d6c36dd96f8d6ef6c6e52c885f8cee4f04961cf02e682e20'
    report = json.loads(INPUT.read_text())
    sources, rows = [], []
    for comparison in report['comparisons']:
        x, z = comparison['X'], comparison['Z']
        path = Path('/tmp/ember-tops-cool-family-v1') / f'seed-x{round(x*1000):03d}-z{round(z*1000):03d}'
        receipt = json.loads((path/'receipt.json').read_text())
        assert (receipt['X'], receipt['Z']) == (x, z)
        assert digest(path/receipt['request']) == receipt['request_sha256']
        request = json.loads((path/receipt['request']).read_text())
        assert request['plasnu'] == 'on' and request['lib'] == 'new' and request['datype'] == 'gray'
        _, _, cells, excluded = read(path/receipt['file'], receipt, dimensions=tuple(receipt['dimensions']))
        sources.append(dict(path=str(path), receipt_sha256=digest(path/'receipt.json'),
                            source_sha256=receipt['sha256'], request_sha256=receipt['request_sha256']))
        matched = 0
        for cell in comparison['rows']:
            key = cell['T_keV'], cell['rho']
            if key not in cells or key in excluded:
                continue
            on, off, mesa = cells[key], cell['tops_uncut'], cell['mesa_four_nodes']
            assert min(on, off, mesa) > 0
            rows.append(dict(X=x, Z=z, T_keV=key[0], rho=key[1], logR=cell['logR'],
                             tops_cutoff=on, tops_uncut=off, mesa=mesa,
                             mesa_relative_to_uncut=mesa/off-1,
                             mesa_relative_to_cutoff=mesa/on-1,
                             cutoff_relative_to_uncut=on/off-1))
            matched += 1
        assert matched > 0
        sources[-1]['exact_matches'] = matched
    result = dict(outcome='completed_exact_state_comparison', accepted_for_stellar_opacity=False,
                  input_sha256={str(p):digest(p) for p in [Path(__file__),INPUT,ROOT/'scripts/import_tops_composition.py']},
                  sources=sources, rows=rows, summaries=summarize(rows),
                  limitations=['All matches have exactly equal printed T, density, X and Z; substituted densities are excluded.',
                               'This reuses the four-point MESA values; their interpolation and 25-versus-21-element mixture differences remain.',
                               'The TOPS cutoff-on option is not the independently calculated refractive transport correction.',
                               'Better agreement with one convention does not by itself identify the MESA table-generation setting.'])
    with OUTPUT.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(report=str(OUTPUT),sha256=digest(OUTPUT),matches=len(rows),summaries=result['summaries']),indent=2))


if __name__ == '__main__':
    main()
