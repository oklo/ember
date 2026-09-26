#!/usr/bin/env python3
"""Fetch bounded native uncut grey opacity rectangles for verified mixtures.

Keep source substitutions explicitly excluded. Preserve requests, calculation
handles and responses so a completed batch is never submitted again. No
cutoff-on overlap equality or stellar-input acceptance is implied.
"""
import argparse
import json
import math
from pathlib import Path
import secrets
import urllib.parse
import urllib.request

from fetch_tops_composition import Form, Text, digest
from import_tops_composition import read


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('plan', type=Path)
    p.add_argument('--request-index', type=int, nargs='+')
    a = p.parse_args()
    spec = json.loads(a.plan.read_text())
    manifest_path = Path(spec['baseline_manifest'])
    if digest(manifest_path) != spec['baseline_manifest_sha256']:
        raise ValueError('baseline manifest changed')
    baseline = json.loads(manifest_path.read_text())
    indices = list(range(len(spec['requests']))) if a.request_index is None else a.request_index
    if len(indices) != len(set(indices)) or any(i < 0 or i >= len(spec['requests']) for i in indices):
        raise ValueError('invalid request indices')
    cap = spec['maximum_response_bytes']
    if not 0 < cap <= 20000000:
        raise ValueError('invalid response cap')

    def post(endpoint, values):
        request = urllib.request.Request('https://aphysics2.lanl.gov/' + endpoint,
                                        data=urllib.parse.urlencode(values).encode())
        with urllib.request.urlopen(request, timeout=180) as response:
            raw = response.read(cap + 1)
        if len(raw) > cap:
            raise ValueError('response exceeds data cap')
        return raw.decode()

    for index in indices:
        job = spec['requests'][index]
        matches = [r for r in baseline['planes'] if (r['X'], r['Z']) == (job['X'], job['Z'])]
        if len(matches) != 1:
            raise ValueError('expected one pinned mixture')
        old = matches[0]
        request_path = manifest_path.parent / old['request']
        if digest(request_path) != old['request_sha256']:
            raise ValueError('baseline request changed')
        native_t, _, _, _ = read(manifest_path.parent / old['file'], old,
                                 dimensions=tuple(old.get('dimensions', [50, 71])))
        tt, rr = job['temperatures_keV'], job['densities_g_cm3']
        if (not tt or tt != sorted(set(tt)) or any(t not in native_t for t in tt)
                or not 2 <= len(rr) <= 100 or rr != sorted(set(rr)) or rr[0] <= 0
                or not all(math.isfinite(r) for r in rr)
                or len(tt) * len(rr) > spec['maximum_states_per_request']):
            raise ValueError('invalid request coordinates')
        expected_rho = [rr[0] * (rr[-1] / rr[0]) ** (i / (len(rr) - 1)) for i in range(len(rr))]
        if any(not math.isclose(x, y, rel_tol=1e-12) for x, y in zip(rr, expected_rho, strict=True)):
            raise ValueError('listed densities do not match logarithmic range')
        root = Path(job['work'])
        root.mkdir(parents=True, exist_ok=True)
        recipe = {'plan_sha256': digest(a.plan), 'script_sha256': digest(Path(__file__)),
                  'baseline_manifest_sha256': digest(manifest_path), 'baseline': old, 'job': job}
        recipe_path = root / 'recipe.json'
        if recipe_path.exists() and json.loads(recipe_path.read_text()) != recipe:
            raise ValueError('changed recipe requires another directory')
        if not recipe_path.exists():
            recipe_path.write_text(json.dumps(recipe, indent=2) + '\n')
        q = {**json.loads(request_path.read_text()), 'datype': 'gray', 'plasnu': 'off',
             'tgrid': 'specific', 'temps': ' '.join(map(str, tt)), 'rgrid': 'range',
             'rlow': str(rr[0]), 'rup': str(rr[-1]), 'nr': str(len(rr)), 'rspace': 'log'}
        local_request = root / 'request.json'
        if local_request.exists():
            saved = json.loads(local_request.read_text())
            if {k: v for k, v in saved.items() if k != 'userid'} != {k: v for k, v in q.items() if k != 'userid'}:
                raise ValueError('saved request differs')
            q = saved
        else:
            q['userid'] = 'e' + ''.join(secrets.choice('abcdefghijklmnopqrstuvwxyz') for _ in range(2))
            local_request.write_text(json.dumps(q, indent=2) + '\n')
        receipt_path = root / 'receipt.json'
        source_path = root / 'source.txt'
        if receipt_path.exists():
            receipt = json.loads(receipt_path.read_text())
            if receipt['request_sha256'] != digest(local_request):
                raise ValueError('completed request changed')
        else:
            submit = root / 'submit.html'
            if not submit.exists():
                submit.write_text(post('submit', q))
            form = Form()
            form.feed(submit.read_text())
            if not form.done or form.data.get('output') != 'tabcol':
                raise ValueError('unexpected source handle')
            # Reuse a completed response even after a local parsing failure.
            result_path = root / 'results.html'
            if not result_path.exists():
                result_path.write_text(post('results', form.data))
            parser = Text()
            parser.feed(result_path.read_text())
            text = '\n'.join(v.strip() for v in ''.join(parser.parts).replace('\u00a0', ' ').splitlines() if v.strip()) + '\n'
            if source_path.exists() and source_path.read_text() != text:
                raise ValueError('retained source differs from the response')
            if not source_path.exists():
                source_path.write_text(text)
            receipt = {**old, 'file': source_path.name, 'sha256': digest(source_path),
                       'request': local_request.name, 'request_sha256': digest(local_request),
                       'dimensions': [len(tt), len(rr)], 'plasma_cutoff': 'off',
                       'validation_role': 'native_uncut_grey_candidate'}
        actual_t, actual_r, cells, excluded = read(source_path, receipt, dimensions=(len(tt), len(rr)))
        if actual_t != tt or any(abs(x / y - 1) > 5e-5 for x, y in zip(actual_r, rr, strict=True)):
            raise ValueError('returned coordinates differ from plan')
        if not receipt_path.exists():
            receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
        print(json.dumps({'request_index': index, 'work': str(root), 'states': len(cells),
                          'source_exclusions': len(excluded)}), flush=True)


if __name__ == '__main__':
    main()
