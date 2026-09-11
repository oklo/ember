#!/usr/bin/env python3
"""Request a separate high-density TOPS coverage diagnostic for a pinned mixture.

Retains the original temperature grid and includes the original upper density
as an overlap check. Substituted source states remain explicitly excluded.
No runtime opacity table is assembled or installed.
"""
import argparse
import json
import math
from pathlib import Path
import secrets
import time
import urllib.parse
import urllib.request

from fetch_tops_composition import Form, Text, digest
from import_tops_composition import read, validate_mixture


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('baseline_manifest', type=Path)
    p.add_argument('work', type=Path)
    p.add_argument('--hydrogen', type=float, required=True)
    p.add_argument('--metallicity', type=float, default=.02)
    p.add_argument('--maximum-density', type=float, default=1e6)
    p.add_argument('--density-points', type=int, default=11)
    a = p.parse_args()
    spec = json.loads(a.baseline_manifest.read_text())
    matches = [r for r in spec['planes'] if r['X'] == a.hydrogen and r['Z'] == a.metallicity]
    if len(matches) != 1:
        raise ValueError('expected one pinned source mixture')
    baseline = matches[0]
    old_source = a.baseline_manifest.parent/baseline['file']
    old_request = a.baseline_manifest.parent/baseline['request']
    if digest(old_request) != baseline['request_sha256']:
        raise ValueError('baseline request changed')
    # A second adjacent density interval can use the first interval's
    # verified source as its overlap reference, without querying that
    # interval again. Original source manifests still require 50 by 71.
    tt, rr, cells, excluded = read(old_source, baseline,
                                   dimensions=tuple(baseline.get('dimensions', [50, 71])))
    if (not math.isfinite(a.maximum_density) or a.maximum_density <= rr[-1]
            or not 2 <= a.density_points <= 100):
        raise ValueError('invalid extended density grid')
    dimensions = (len(tt), a.density_points)
    a.work.mkdir(parents=True, exist_ok=True)
    recipe = {'baseline_manifest': str(a.baseline_manifest.resolve()),
              'baseline_manifest_sha256': digest(a.baseline_manifest),
              'baseline': baseline, 'maximum_density': a.maximum_density,
              'density_points': a.density_points, 'script_sha256': digest(__file__)}
    saved = a.work/'recipe.json'
    if saved.exists() and json.loads(saved.read_text()) != recipe:
        raise ValueError('changed request requires a new work directory')
    saved.write_text(json.dumps(recipe, indent=2)+'\n')
    request_path = a.work/'request.json'
    request = json.loads(old_request.read_text())
    request.update(rgrid='range', rlow=format(rr[-1], '.17g'),
                   rup=format(a.maximum_density, '.17g'), nr=str(a.density_points), rspace='log')
    if request_path.exists():
        pending = json.loads(request_path.read_text())
        if {k: v for k, v in pending.items() if k != 'userid'} != {
                k: v for k, v in request.items() if k != 'userid'}:
            raise ValueError('pending request changed')
        request = pending
    else:
        request['userid'] = 'e'+''.join(secrets.choice('abcdefghijklmnopqrstuvwxyz') for _ in range(2))
        request_path.write_text(json.dumps(request, indent=2)+'\n')
    source_path, receipt_path = a.work/'source.txt', a.work/'receipt.json'
    if receipt_path.exists():
        receipt = json.loads(receipt_path.read_text())
        if receipt['request_sha256'] != digest(request_path):
            raise ValueError('completed request changed')
        result = read(source_path, receipt, dimensions=dimensions)
    else:
        if source_path.exists():
            raise ValueError('incomplete source result retained; use a new work directory')
        def post(endpoint, data):
            req = urllib.request.Request('https://aphysics2.lanl.gov/'+endpoint,
                                         data=urllib.parse.urlencode(data).encode())
            return urllib.request.urlopen(req, timeout=180).read().decode()
        # Preserve the exact request before contacting the calculation service.
        submitted = post('submit', request)
        (a.work/'submit.html').write_text(submitted)
        form = Form(); form.feed(submitted)
        if not form.done or form.data.get('output') != 'tabcol':
            raise ValueError('unexpected TOPS response form')
        for attempt in range(6):
            html = post('results', form.data)
            (a.work/f'result-{attempt}.html').write_text(html)
            parser = Text(); parser.feed(html)
            text = '\n'.join(v.strip() for v in ''.join(parser.parts).replace('\u00a0',' ').splitlines() if v.strip())+'\n'
            try:
                validate_mixture(text, baseline, dimensions=dimensions)
            except (ValueError, IndexError, KeyError) as error:
                print(f'rejected source response: {error}', flush=True)
                time.sleep(5)
            else:
                break
        else:
            raise ValueError('service did not return the requested mixture and grid size')
        source_path.write_text(text)
        receipt = {**baseline, 'file': source_path.name, 'sha256': digest(source_path),
                   'request': request_path.name, 'request_sha256': digest(request_path),
                   'dimensions': list(dimensions)}
        result = read(source_path, receipt, dimensions=dimensions)
    temperatures, densities, new_cells, new_excluded = result
    wanted = [rr[-1]*(a.maximum_density/rr[-1])**(i/(a.density_points-1)) for i in range(a.density_points)]
    if temperatures != tt or any(abs(x/y-1) > 5e-5 for x, y in zip(densities, wanted, strict=True)):
        raise ValueError('returned temperature or density coordinates differ from request')
    overlap = []
    for t in tt:
        point = (t, rr[-1])
        if point not in excluded and point not in new_excluded:
            overlap.append({'temperature_keV': t, 'density_g_cm3': rr[-1],
                            'relative_difference': new_cells[point]/cells[point]-1})
    if not overlap or any(row['relative_difference'] != 0 for row in overlap):
        raise ValueError('retained source overlap changed')
    if not receipt_path.exists():
        receipt_path.write_text(json.dumps(receipt, indent=2)+'\n')
    report = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
              'recipe': recipe, 'request_sha256': digest(request_path),
              'source_sha256': digest(source_path), 'dimensions': list(dimensions),
              'overlap': overlap, 'density_g_cm3': densities,
              'coverage': [{'temperature_keV': t,
                  'original_densities_g_cm3': [r for r in densities if (t, r) not in new_excluded],
                  'substituted_densities_g_cm3': [r for r in densities if (t, r) in new_excluded]} for t in tt]}
    (a.work/'coverage.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'overlap_states_exact': len(overlap), 'substituted_states': len(new_excluded),
                      'coverage': str(a.work/'coverage.json')}))


if __name__ == '__main__':
    main()
