#!/usr/bin/env python3
"""Fetch bounded batches of uncut TOPS group data with explicit source identity."""
import argparse
import json
import math
from pathlib import Path
import secrets
import urllib.parse
import urllib.request

from fetch_tops_composition import Form, Text, digest
from import_tops_composition import read
from tops_groups import source_groups


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('plan', type=Path)
    p.add_argument('--request-index', type=int, nargs='+')
    a = p.parse_args(); spec = json.loads(a.plan.read_text())
    seed = Path(spec['baseline']); original = json.loads((seed/'receipt.json').read_text())
    native, _, _, _ = read(seed/'source.txt', original, dimensions=tuple(original['dimensions']))
    if digest(seed/'request.json') != original['request_sha256']:
        raise ValueError('seed request changed')
    base = json.loads((seed/'request.json').read_text())
    indices = list(range(len(spec['requests']))) if a.request_index is None else a.request_index
    if len(set(indices)) != len(indices) or any(i < 0 or i >= len(spec['requests']) for i in indices):
        raise ValueError('invalid request selection')
    cap = spec['maximum_response_bytes']

    def post(endpoint, values):
        request = urllib.request.Request('https://aphysics2.lanl.gov/'+endpoint,
            data=urllib.parse.urlencode(values).encode())
        with urllib.request.urlopen(request, timeout=180) as response:
            raw = response.read(cap+1)
        if len(raw) > cap:
            raise ValueError('response exceeds data cap')
        return raw.decode()

    for index in indices:
        job = spec['requests'][index]; tt, rr = job['temperatures_keV'], job['densities_g_cm3']
        if (any(t not in native for t in tt) or not 2 <= len(rr) <= 100
                or len(set(tt)) != len(tt) or tt != sorted(tt) or rr != sorted(set(rr))):
            raise ValueError('invalid source coordinates')
        spacing = job.get('density_spacing', 'linear')
        expected = [rr[0]+(rr[-1]-rr[0])*i/(len(rr)-1) if spacing == 'linear'
                    else rr[0]*(rr[-1]/rr[0])**(i/(len(rr)-1)) for i in range(len(rr))]
        if spacing not in ('linear', 'log') or any(not math.isclose(x, y, rel_tol=1e-12) for x, y in zip(rr, expected, strict=True)):
            raise ValueError('density range differs from listed values')
        n = spec['photon_boundaries']
        if not 2 <= n <= 1000 or n*len(tt)*len(rr) > spec['maximum_group_states']:
            raise ValueError('unbounded group request')
        root = Path(job['work']); root.mkdir(parents=True, exist_ok=True)
        recipe = {'plan_sha256': digest(a.plan), 'script_sha256': digest(Path(__file__)),
                  'seed_receipt_sha256': digest(seed/'receipt.json'), 'job': job}
        recipe_path = root/'recipe.json'
        if recipe_path.exists() and json.loads(recipe_path.read_text()) != recipe:
            raise ValueError('changed group recipe needs a separate directory')
        recipe_path.write_text(json.dumps(recipe, indent=2)+'\n')
        q = {**base, 'datype': 'groups', 'plasnu': 'off', 'tgrid': 'specific', 'temps': ' '.join(map(str, tt)),
             'rgrid': 'range', 'rlow': str(rr[0]), 'rup': str(rr[-1]), 'nr': str(len(rr)), 'rspace': spacing,
             'egrid': 'range', 'egplow': str(spec['photon_min_keV']), 'egphigh': str(spec['photon_max_keV']),
             'ngpengs': str(n), 'espace': 'log'}
        request_path = root/'request.json'
        if request_path.exists():
            saved = json.loads(request_path.read_text())
            if any(saved[k] != v for k, v in q.items() if k != 'userid'):
                raise ValueError('pending group request changed')
            q = saved
        else:
            q['userid'] = 'e'+''.join(secrets.choice('abcdefghijklmnopqrstuvwxyz') for _ in range(2))
            request_path.write_text(json.dumps(q, indent=2)+'\n')
        receipt_path = root/'receipt.json'
        if not receipt_path.exists():
            submitted = root/'submit.html'
            if not submitted.exists(): submitted.write_text(post('submit', q))
            form = Form(); form.feed(submitted.read_text())
            if not form.done or form.data.get('output') != 'tabcol':
                raise ValueError('unexpected source handle')
            html = post('results', form.data); (root/'results.html').write_text(html)
            parser = Text(); parser.feed(html)
            text = '\n'.join(v.strip() for v in ''.join(parser.parts).replace('\u00a0',' ').splitlines() if v.strip())+'\n'
            path = root/'source.txt'; path.write_text(text)
            receipt = {**original, 'file': path.name, 'sha256': digest(path), 'request': request_path.name,
                       'request_sha256': digest(request_path), 'dimensions': [len(tt), len(rr)],
                       'validation_role': 'group_transport_control', 'plasma_cutoff': 'off'}
            returned_t, returned_r, _, excluded = read(path, receipt, dimensions=(len(tt), len(rr)))
            if returned_t != tt or excluded or any(abs(x/y-1) > 5e-5 for x, y in zip(returned_r, rr, strict=True)):
                raise ValueError('group source returned different or substituted coordinates')
            receipt_path.write_text(json.dumps(receipt, indent=2)+'\n')
        _, data, _ = source_groups(root)
        if sorted({k[0] for k in data}) != tt or any(abs(x/y-1) > 5e-5 for x, y in zip(sorted({k[1] for k in data}), rr, strict=True)):
            raise ValueError('completed group coordinates differ')
        print(json.dumps({'work': str(root), 'states': len(data), 'groups_per_state': n-1}), flush=True)


if __name__ == '__main__': main()
