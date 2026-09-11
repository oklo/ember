#!/usr/bin/env python3
"""Resolve a small TOPS density interval and isolate the plasma-cutoff option.

These are separate source diagnostics. Neither calculation is installed in
stellar opacity. Every request, source reply and composition check is retained.
"""
import argparse
import json
from pathlib import Path
import secrets
import time
import urllib.parse
import urllib.request

from fetch_tops_composition import Form, Text, digest
from import_tops_composition import read, validate_mixture


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('reference', type=Path, help='completed source work with receipt/request/source.txt')
    p.add_argument('output', type=Path)
    p.add_argument('--temperature-kev', type=float, required=True)
    p.add_argument('--density-min', type=float, required=True)
    p.add_argument('--density-max', type=float, required=True)
    p.add_argument('--density-points', type=int, default=91)
    p.add_argument('--frequency-dependent', action='store_true')
    a = p.parse_args()
    if (a.output.exists() or not 0 < a.density_min < a.density_max
            or not 2 <= a.density_points <= (6 if a.frequency_dependent else 100)):
        raise ValueError('use a new diagnostic directory and a positive density interval')
    receipt = json.loads((a.reference/'receipt.json').read_text())
    reference_request = json.loads((a.reference/'request.json').read_text())
    if digest(a.reference/'request.json') != receipt['request_sha256']:
        raise ValueError('reference request changed')
    tt, rr, cells, excluded = read(a.reference/'source.txt', receipt,
                                  dimensions=tuple(receipt['dimensions']))
    if a.temperature_kev not in tt:
        raise ValueError('temperature must be an existing source isotherm')
    controls = [r for r in rr if a.density_min <= r <= a.density_max
                and (a.temperature_kev, r) not in excluded]
    if not controls:
        raise ValueError('density interval needs a converged reference point')
    densities = [a.density_min+(a.density_max-a.density_min)*i/(a.density_points-1)
                 for i in range(a.density_points)]
    if any(not any(abs(r/v-1) < 1e-12 for v in densities) for r in controls):
        raise ValueError('linear diagnostic grid must include the reference controls')
    dimensions = (1, len(densities))
    inputs = {str(path.resolve()): digest(path) for path in
              [a.reference/'receipt.json', a.reference/'request.json', a.reference/'source.txt', Path(__file__)]}
    a.output.mkdir(parents=True)
    results = []
    for cutoff in ['on', 'off']:
        work = a.output/('cutoff-'+cutoff); work.mkdir()
        request = {**reference_request, 'userid': 'e'+''.join(secrets.choice('abcdefghijklmnopqrstuvwxyz') for _ in range(2)),
                   'tgrid': 'specific', 'temps': format(a.temperature_kev, '.17g'),
                   'rgrid': 'range', 'rlow': format(a.density_min, '.17g'),
                   'rup': format(a.density_max, '.17g'), 'nr': str(len(densities)), 'rspace': 'linear',
                   'plasnu': cutoff, 'datype': 'cont' if a.frequency_dependent else 'gray'}
        (work/'request.json').write_text(json.dumps(request, indent=2)+'\n')
        def post(endpoint, data):
            req = urllib.request.Request('https://aphysics2.lanl.gov/'+endpoint,
                                         data=urllib.parse.urlencode(data).encode())
            return urllib.request.urlopen(req, timeout=180).read().decode()
        html = post('submit', request)
        (work/'submit.html').write_text(html)
        form = Form(); form.feed(html)
        if not form.done or form.data.get('output') != 'tabcol':
            raise ValueError('unexpected source response form')
        for attempt in range(6):
            html = post('results', form.data); (work/f'results-{attempt}.html').write_text(html)
            parser = Text(); parser.feed(html)
            text = '\n'.join(v.strip() for v in ''.join(parser.parts).replace('\u00a0', ' ').splitlines() if v.strip())+'\n'
            try:
                validate_mixture(text, receipt, dimensions=dimensions)
            except (ValueError, IndexError, KeyError):
                time.sleep(5)
            else:
                break
        else:
            raise ValueError('source did not return the requested mixture or grid')
        source = work/'source.txt'; source.write_text(text)
        record = {**receipt, 'file': source.name, 'sha256': digest(source),
                  'request': 'request.json', 'request_sha256': digest(work/'request.json'),
                  'dimensions': list(dimensions), 'diagnostic_plasma_cutoff': cutoff}
        (work/'receipt.json').write_text(json.dumps(record, indent=2)+'\n')
        t, r, k, e = read(source, record, dimensions=dimensions)
        if (t != [a.temperature_kev] or e
                or any(abs(x/y-1) > 5e-5 for x, y in zip(r, densities, strict=True))):
            raise ValueError('source query coordinates changed or were substituted')
        if cutoff == 'on' and any(k[a.temperature_kev, rho] != cells[a.temperature_kev, rho] for rho in controls):
            raise ValueError('original source control changed')
        results.append({'plasma_cutoff': cutoff, 'density': r, 'kappa': [k[a.temperature_kev, rho] for rho in r],
                        'source_sha256': digest(source), 'request_sha256': digest(work/'request.json')})
    if results[0]['density'] != results[1]['density'] or any(digest(Path(p)) != h for p, h in inputs.items()):
        raise ValueError('source coordinates or inputs changed between controls')
    report = {'scope': __doc__, 'reference_input_sha256': inputs, 'reference_controls': controls,
              'temperature_keV': a.temperature_kev, 'results': results,
              'maximum_cutoff_relative_effect': max(abs(x/y-1) for x, y in
                  zip(results[0]['kappa'], results[1]['kappa'], strict=True)),
              'accepted_for_stellar_opacity': False}
    (a.output/'diagnostic.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'density_points': len(densities), 'maximum_cutoff_relative_effect': report['maximum_cutoff_relative_effect']}))


if __name__ == '__main__':
    main()
