#!/usr/bin/env python3
"""Import supplied pure-H WD boundary values without the MESA gray blend.

Input is MESA's wd_25.tbl (Rohrmann atmospheres), containing Pgas and T
at Rosseland tau=25.1188. The output explicitly approximates small helium
and metal fractions by the pure-H boundary; the interior retains its actual
composition. This importer does not select the table for stellar evolution.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import re

from write_scientific_result import write_result


def read_table(path):
    lines = Path(path).read_text().splitlines()
    assert lines[0].strip() == '#Table Version   5', 'unrecognized source version'
    starts = [i for i, s in enumerate(lines) if s.startswith('#Teff(K)|')]
    assert len(starts) == 2 and 'Pgas@' in lines[starts[0]] and 'T@' in lines[starts[1]]
    groups = []
    for k, start in enumerate(starts):
        gravity = list(map(float, re.findall(r'log g\s*=\s*([0-9.]+)', lines[start])))
        end = starts[k+1] if k+1 < len(starts) else len(lines)
        rows = [list(map(float, s.split())) for s in lines[start+1:end] if s.strip()]
        assert len(gravity) >= 2 and len(rows) >= 2
        assert all(len(row) == len(gravity)+1 for row in rows)
        assert all(math.isfinite(v) and v > 0 for row in rows for v in row)
        assert all(b > a for a, b in zip(gravity, gravity[1:]))
        teff = [row[0] for row in rows]
        assert all(b > a for a, b in zip(teff, teff[1:]))
        groups.append((teff, gravity, [row[1:] for row in rows]))
    assert groups[0][:2] == groups[1][:2], 'pressure/temperature coordinates differ'
    teff, gravity, pressure = groups[0]
    temperature = groups[1][2]
    valid = list(map(int, lines[1].split('VALID RANGE:', 1)[1].split()))
    # This source's descriptive header has six surplus copies of381. MESA
    # reads one limit per actual gravity column; require every entry, including
    # the surplus ones, to declare the same complete temperature coverage.
    assert len(valid) >= len(gravity) and all(v == len(teff) for v in valid), \
        'partial table needs an explicit source mask'
    return teff, gravity, temperature, pressure


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--maximum-helium', type=float, default=.005)
    p.add_argument('--maximum-metallicity', type=float, default=1e-20)
    a = p.parse_args()
    assert 0 < a.maximum_helium <= .005 and 0 < a.maximum_metallicity <= 1e-20
    assert not a.output.exists(), 'do not overwrite an existing source import'
    ts, gs, temperatures, pressures = read_table(a.source)
    # This is a physical pure-H proxy, not additional calculated compositions.
    pattern = [0.20138655436774022, 0.0022627702737948339,
               0.055673552701166383, 0.4606771226572986, 0.28]
    approximation = ('Supplied pure-H Rohrmann WD atmosphere; matching at tau25.1188. '
        'Boundary temperature and gas pressure neglect the explicitly bounded helium '
        'and metal content; density uses the actual interior EOS and composition. '
        'Logarithmic bilinear interpolation retains every supplied node. No gray blend '
        'or extrapolation. Composition, matching-depth and stellar response need review.')
    lines = ['EMBER_HYDROGEN_DOMINATED_ATMOSPHERE 1',
        'source '+json.dumps('Rohrmann pure-H WD boundary supplied with MESA; wd_25.tbl version5'),
        'approximation '+json.dumps(approximation), 'basis baryon_mass', 'tau 25.1188',
        f'maximum_helium {a.maximum_helium:.17g}', 'source_helium 0',
        'metal_pattern '+' '.join(format(x, '.17g') for x in pattern),
        f'metallicity 2 0 {a.maximum_metallicity:.17g}',
        'log_teff '+str(len(ts))+' '+' '.join(format(math.log10(t), '.17g') for t in ts),
        'log_g '+str(len(gs))+' '+' '.join(format(g, '.17g') for g in gs), 'data']
    for _ in range(2):
        for tr, pr in zip(temperatures, pressures):
            lines += [f'{math.log10(t):.17g} {math.log10(pg):.17g}' for t, pg in zip(tr, pr)]
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text('\n'.join(lines)+'\n')
    sha = lambda f: hashlib.sha256(Path(f).read_bytes()).hexdigest()
    write_result(a.output.with_suffix('.provenance.json'), dict(
        selected_for_evolution=False, source=str(a.source.resolve()), source_sha256=sha(a.source),
        importer_sha256=sha(__file__), table_sha256=sha(a.output),
        unique_source_states=len(ts)*len(gs), teff=ts, logg=gs,
        matching_tau=25.1188, source_helium=0., source_metallicity=0.,
        maximum_helium=a.maximum_helium, maximum_metallicity=a.maximum_metallicity,
        approximation=approximation,
        references=['https://arxiv.org/abs/1209.2452', 'https://arxiv.org/abs/1301.0319'],
        limits=['The supplied MESA grid extends to logg5.5; the2012 atmosphere paper describes logg6.5..9.5. The MESA source grid, its preparation code and the MESA2013 description support the imported extent; no lower-gravity values are generated by this importer.',
                'The two metallicity planes contain identical pure-H values under the stated trace-composition approximation. They are not distinct solved atmospheres.',
                'These supplied matching states do not include full atmosphere profiles or independent source-convergence records. They provide an independent published prescription rather than a verification of the Ember atmosphere solver.']))


if __name__ == '__main__':
    main()
