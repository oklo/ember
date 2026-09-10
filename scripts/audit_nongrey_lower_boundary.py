#!/usr/bin/env python3
"""Compare independently validated atmosphere solutions at different columns.

The final optical depths, not TAULAS or the seed trimming parameter, define
the tested lower boundaries. This local comparison does not bound all grid
cells, spectral errors or the physical uncertainty of the source model.
"""
import argparse
import json
import math
from pathlib import Path

from assemble_nongrey_grid import load_continuation, physical_identity
from prepare_nongrey_sources import digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('output', type=Path)
    p.add_argument('models', type=Path, nargs='+')
    p.add_argument('--tolerance', type=float, default=1e-4,
                   help='relative T/Pgas criterion, default 0.01 percent')
    a = p.parse_args()
    if len(a.models) < 2 or not 0 < a.tolerance < .01:
        raise ValueError('need at least two models and a finite positive criterion')
    rows = []; inputs = {}; identity = None; coordinate = None; reference = None
    for work in a.models:
        key, state, spec, record, files = load_continuation(work, identity)
        source = record['continuation_provenance']
        if identity is None:
            identity = physical_identity(spec, source['prepared'])
            coordinate = key
            reference = {k: spec[k] for k in ['depths', 'atmosphere_frequencies', 'convection_derivative_step']}
            opacity = source['opacity_sha256']
        if key != coordinate or any(spec[k] != v for k, v in reference.items()) or source['opacity_sha256'] != opacity:
            raise ValueError('comparison must isolate bottom column at one source state and resolution')
        inputs.update(files)
        rows.append({'work': str(work.resolve()), 'state': state,
                     'seed_bottom_tau': source.get('initial_bottom_tau'),
                     'initializer': source['initializer'], 'validation_sha256': digest(work/'final/validated.json')})
    rows.sort(key=lambda r: r['state']['column_mass_range'][1])
    deepest = rows[-1]['state']
    if deepest['optical_depth_range'][1]/rows[0]['state']['optical_depth_range'][1] < 1.5:
        raise ValueError('lower-boundary depths span less than a factor of 1.5')
    if any(a['state']['column_mass_range'][1] >= b['state']['column_mass_range'][1] for a, b in zip(rows, rows[1:])):
        raise ValueError('repeated bottom column is not an independent depth comparison')
    for row in rows:
        row['relative_difference_from_deepest'] = {k: row['state'][k]/deepest[k]-1 for k in ['T', 'Pgas', 'source_density']}
    maximum = max(abs(r['relative_difference_from_deepest'][k]) for r in rows for k in ['T', 'Pgas'])
    if any(digest(path) != expected for path, expected in inputs.items()):
        raise ValueError('source changed during lower-boundary audit')
    report = {'scope': __doc__, 'coordinates': coordinate, 'models': rows,
              'relative_matching_state_tolerance': a.tolerance,
              'maximum_matching_state_relative_difference': maximum,
              'passed': math.isfinite(maximum) and maximum <= a.tolerance,
              'input_files_sha256': inputs, 'script_sha256': digest(__file__)}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: report[k] for k in ['coordinates', 'maximum_matching_state_relative_difference', 'passed']}))
    if not report['passed']:
        raise SystemExit('lower-boundary comparison failed its matching-state criterion')


if __name__ == '__main__':
    main()
