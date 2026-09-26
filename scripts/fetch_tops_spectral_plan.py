#!/usr/bin/env python3
"""Fetch only missing groups in a bounded TOPS spectral plan.

Requests, service handles and responses are retained. A resumed request polls
the same handle. Completed source receipts are verified without another request.
"""
import argparse
import json
import math
from pathlib import Path
import secrets
import time
import urllib.error
import urllib.parse
import urllib.request

from fetch_tops_composition import Form, Text, digest
from import_tops_composition import read, validate_mixture
from audit_tops_spectral_means import source


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('plan', type=Path)
    p.add_argument('--request-index', type=int, nargs='+', help='Run only these independent request groups; retain the complete plan identity.')
    a = p.parse_args()
    spec = json.loads(a.plan.read_text())
    seed = Path(spec['baseline'])
    original = json.loads((seed/'receipt.json').read_text())
    native, _, _, _ = read(seed/'source.txt', original,
                           dimensions=tuple(original['dimensions']))
    if digest(seed/'request.json') != original['request_sha256']:
        raise ValueError('seed request changed')
    base_request = json.loads((seed/'request.json').read_text())
    reused = set()
    for path in spec['reused_sources']:
        root = Path(path)
        old = json.loads((root/'receipt.json').read_text())
        if any(old[k] != original[k] for k in ('X', 'Z', 'metals')):
            raise ValueError('reused mixture differs')
        spectra, _ = source(root)
        reused.update(spectra)
    wanted = {(t, r) for job in spec['requests'] for t in job['temperatures_keV']
              for r in job['densities_g_cm3']}
    if len(wanted) != sum(len(v['temperatures_keV'])*len(v['densities_g_cm3'])
                          for v in spec['requests']) or wanted & reused:
        raise ValueError('plan would duplicate a source calculation')
    cap = spec['maximum_source_bytes_per_request']

    def post(endpoint, data):
        request = urllib.request.Request('https://aphysics2.lanl.gov/'+endpoint,
            data=urllib.parse.urlencode(data).encode())
        with urllib.request.urlopen(request, timeout=180) as response:
            raw = response.read(cap+1)
        if len(raw) > cap:
            raise ValueError('response exceeds the planned data cap')
        return raw.decode()

    indices = list(range(len(spec['requests']))) if a.request_index is None else a.request_index
    if len(set(indices)) != len(indices) or any(i < 0 or i >= len(spec['requests']) for i in indices):
        raise ValueError('invalid or repeated request index')
    for index in indices:
        job = spec['requests'][index]
        tt, rr = job['temperatures_keV'], job['densities_g_cm3']
        if any(t not in native for t in tt) or len(set(rr)) != len(rr) or not 1 <= len(tt)*len(rr) <= 6:
            raise ValueError('non-native temperature, repeated density or more than six spectral states')
        dimensions = len(tt), len(rr)
        root = Path(job['work']); root.mkdir(parents=True, exist_ok=True)
        recipe = {'plan_sha256': digest(a.plan), 'script_sha256': digest(Path(__file__)),
                  'seed_receipt_sha256': digest(seed/'receipt.json'), 'request': job}
        recipe_path = root/'recipe.json'
        if recipe_path.exists() and json.loads(recipe_path.read_text()) != recipe:
            raise ValueError('changed plan requires a separate directory')
        recipe_path.write_text(json.dumps(recipe, indent=2)+'\n')
        request = {**base_request, 'tgrid': 'specific', 'temps': ' '.join(map(str, tt)),
                   'rgrid': 'specific', 'dens': ' '.join(map(str, rr)),
                   'plasnu': 'on', 'datype': 'cont'}
        if 'density_spacing' in job:
            spacing = job['density_spacing']
            if spacing not in ('linear', 'log') or len(rr) < 2 or any(r <= 0 for r in rr):
                raise ValueError('invalid density range')
            expected = [rr[0]+(rr[-1]-rr[0])*i/(len(rr)-1) if spacing == 'linear'
                        else rr[0]*(rr[-1]/rr[0])**(i/(len(rr)-1)) for i in range(len(rr))]
            if any(not math.isclose(x, y, rel_tol=1e-12) for x, y in zip(rr, expected, strict=True)):
                raise ValueError('listed densities differ from the requested range')
            request.update(rgrid='range', rlow=str(rr[0]), rup=str(rr[-1]), nr=str(len(rr)), rspace=spacing)
        request_path = root/'request.json'
        if request_path.exists():
            saved = json.loads(request_path.read_text())
            if any(saved[k] != v for k, v in request.items() if k != 'userid'):
                raise ValueError('pending request changed')
            request = saved
        else:
            request['userid'] = 'e'+''.join(secrets.choice('abcdefghijklmnopqrstuvwxyz') for _ in range(2))
            request_path.write_text(json.dumps(request, indent=2)+'\n')
        receipt_path = root/'receipt.json'
        if not receipt_path.exists():
            submitted = root/'submit.html'
            if not submitted.exists():
                submitted.write_text(post('submit', request))
            form = Form(); form.feed(submitted.read_text())
            if not form.done or form.data.get('output') != 'tabcol':
                raise ValueError('unexpected source handle')
            for attempt in range(6):
                try:
                    response = post('results', form.data)
                except urllib.error.HTTPError as error:
                    (root/f'error-{time.time_ns()}.txt').write_bytes(error.read(cap))
                    # A repeated server error can mean an unsupported spectral
                    # temperature. Preserve it and stop; do not resubmit.
                    raise
                (root/f'results-{time.time_ns()}.html').write_text(response)
                parser = Text(); parser.feed(response)
                text = '\n'.join(v.strip() for v in ''.join(parser.parts).replace('\u00a0', ' ').splitlines() if v.strip())+'\n'
                try:
                    validate_mixture(text, original, dimensions=dimensions)
                except (ValueError, KeyError, IndexError):
                    time.sleep(5)
                else:
                    break
            else:
                raise ValueError('source did not return the requested mixture')
            path = root/'source.txt'; path.write_text(text)
            receipt = {**original, 'file': path.name, 'sha256': digest(path),
                       'request': request_path.name, 'request_sha256': digest(request_path),
                       'dimensions': list(dimensions), 'validation_role': 'spectral_interpolation_control'}
            returned_t, returned_r, _, excluded = read(path, receipt, dimensions=dimensions)
            if returned_t != tt or excluded or any(abs(x/y-1) > 5e-5 for x, y in zip(returned_r, rr, strict=True)):
                raise ValueError('returned coordinates differ or contain substitutions')
            receipt_path.write_text(json.dumps(receipt, indent=2)+'\n')
        spectra, _ = source(root)
        actual_t = sorted({k[0] for k in spectra}); actual_r = sorted({k[1] for k in spectra})
        if actual_t != tt or len(actual_r) != len(rr) or any(abs(x/y-1) > 5e-5 for x, y in zip(actual_r, rr, strict=True)):
            raise ValueError('completed source does not match planned coordinates')
        print(json.dumps({'work': str(root), 'states': len(spectra), 'verified': True}), flush=True)


if __name__ == '__main__':
    main()
