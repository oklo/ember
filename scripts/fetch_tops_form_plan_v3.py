#!/usr/bin/env python3
"""Fetch bounded grey, group or full-spectrum plans using identifiers supplied by the TOPS form.

Retain official forms, scientific requests, calculation handles and results.
Identifiers are allocated serially and checked against a local inventory.
Completed sources are revalidated locally. A failed submitted request reuses
its handle; neither failed scientific checks nor source exclusions are hidden.
"""
import argparse
import fcntl
import os
import re
import uuid
from html.parser import HTMLParser
import http.cookiejar
import json
import math
from pathlib import Path
import time
import urllib.error
import urllib.parse
import urllib.request

from fetch_tops_composition import Form, Text, digest
from import_tops_composition import read
from tops_groups import source_groups
from audit_tops_spectral_means import source as spectral_source


class SubmissionForm(HTMLParser):
    def __init__(self):
        super().__init__()
        self.active = False
        self.fields = {}

    def handle_starttag(self, tag, attrs):
        values = dict(attrs)
        if tag == 'form' and values.get('action') == '/submit':
            self.active = True
        if (self.active and tag.lower() == 'input' and
                values.get('type', '').lower() == 'hidden' and 'name' in values):
            self.fields[values['name']] = values.get('value', '')

    def handle_endtag(self, tag):
        if tag == 'form':
            self.active = False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kind', choices=('grey', 'groups', 'spectrum'))
    parser.add_argument('plan', type=Path)
    parser.add_argument('--request-index', type=int, nargs='+')
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    indices = list(range(len(plan['requests']))) if args.request_index is None else args.request_index
    if len(indices) != len(set(indices)) or any(i < 0 or i >= len(plan['requests']) for i in indices):
        raise ValueError('invalid request selection')
    cap = plan['maximum_response_bytes']
    if not 0 < cap <= 20000000:
        raise ValueError('invalid response cap')
    if args.kind in ('grey', 'spectrum'):
        manifest_path = Path(plan['baseline_manifest'])
        if digest(manifest_path) != plan['baseline_manifest_sha256']:
            raise ValueError('baseline manifest changed')
        manifest = json.loads(manifest_path.read_text())
    else:
        seed = Path(plan['baseline'])

    for index in indices:
        job = plan['requests'][index]
        if args.kind in ('grey', 'spectrum'):
            matches = [r for r in manifest['planes'] if (r['X'], r['Z']) == (job['X'], job['Z'])]
            if len(matches) != 1:
                raise ValueError('expected one pinned mixture')
            old = matches[0]
            request_path = manifest_path.parent / old['request']
            source_path = manifest_path.parent / old['file']
        else:
            old = json.loads((seed / 'receipt.json').read_text())
            request_path, source_path = seed / old['request'], seed / old['file']
        if digest(request_path) != old['request_sha256']:
            raise ValueError('baseline request changed')
        native_t, _, _, _ = read(source_path, old, dimensions=tuple(old.get('dimensions', [50, 71])))
        tt, rr = job['temperatures_keV'], job['densities_g_cm3']
        if (not tt or tt != sorted(set(tt)) or any(t not in native_t for t in tt) or
                not 2 <= len(rr) <= 100 or rr != sorted(set(rr)) or
                any(not math.isfinite(r) or r <= 0 for r in rr)):
            raise ValueError('invalid source coordinates')
        spacing = 'log' if args.kind == 'grey' else job.get('density_spacing', 'linear')
        if spacing not in ('linear', 'log'):
            raise ValueError('invalid density spacing')
        expected = [rr[0] * (rr[-1] / rr[0])**(i / (len(rr) - 1)) if spacing == 'log'
                    else rr[0] + (rr[-1] - rr[0]) * i / (len(rr) - 1) for i in range(len(rr))]
        if any(not math.isclose(x, y, rel_tol=1e-12) for x, y in zip(rr, expected, strict=True)):
            raise ValueError('density range differs from listed values')
        count = len(tt) * len(rr)
        if args.kind == 'grey' and count > plan['maximum_states_per_request']:
            raise ValueError('unbounded grey request')
        if args.kind == 'spectrum' and count > 6:
            raise ValueError('full-spectrum request exceeds six-state limit')
        q = {**json.loads(request_path.read_text()), 'datype': {'grey':'gray', 'groups':'groups', 'spectrum':'cont'}[args.kind],
             'plasnu': 'on' if args.kind == 'spectrum' else 'off', 'tgrid': 'specific', 'temps': ' '.join(map(str, tt)),
             'rgrid': 'range', 'rlow': str(rr[0]), 'rup': str(rr[-1]), 'nr': str(len(rr)), 'rspace': spacing}
        if args.kind == 'groups':
            n = plan['photon_boundaries']
            if not 2 <= n <= 1000 or count > 72 or n * count > plan['maximum_group_states']:
                raise ValueError('group request exceeds verified batch size')
            q.update(egrid='range', egplow=str(plan['photon_min_keV']),
                     egphigh=str(plan['photon_max_keV']), ngpengs=str(n), espace='log')
        root = Path(job['work'])
        root.mkdir(parents=True, exist_ok=True)
        recipe = {'kind': args.kind, 'plan_sha256': digest(args.plan), 'script_sha256': digest(Path(__file__)),
                  'baseline': old, 'job': job, 'request_timeout_seconds': 75,
                  'input_sha256': {str(p.resolve()): digest(p) for p in
                                   (request_path, source_path, Path('scripts/fetch_tops_composition.py'),
                                    Path('scripts/import_tops_composition.py'), Path('scripts/tops_groups.py'),
                                    Path('scripts/audit_tops_spectral_means.py'))}}
        recipe_path = root / 'recipe.json'
        if recipe_path.exists() and json.loads(recipe_path.read_text()) != recipe:
            raise ValueError('changed recipe requires another directory')
        if not recipe_path.exists():
            recipe_path.write_text(json.dumps(recipe, indent=2) + '\n')
        jar = http.cookiejar.CookieJar()
        opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(jar))
        calls_path = root / 'http_calls.json'
        calls = json.loads(calls_path.read_text()) if calls_path.exists() else []

        def fetch(stage, url, values=None):
            start = time.monotonic()
            req = urllib.request.Request(url, data=None if values is None else urllib.parse.urlencode(values).encode(), headers={'Cache-Control': 'no-cache'})
            try:
                with opener.open(req, timeout=75) as response:
                    raw, status = response.read(cap + 1), response.status
            except urllib.error.HTTPError as error:
                raw = error.read(cap + 1)
                (root / f'{stage}-error-{len(calls):03d}.html').write_bytes(raw)
                calls.append({'stage': stage, 'status': error.code, 'elapsed_seconds': time.monotonic() - start})
                calls_path.write_text(json.dumps(calls, indent=2) + '\n')
                raise
            if len(raw) > cap:
                raise ValueError('response exceeds data cap')
            calls.append({'stage': stage, 'status': status, 'elapsed_seconds': time.monotonic() - start,
                          'cookie_count': len(jar)})
            calls_path.write_text(json.dumps(calls, indent=2) + '\n')
            return raw.decode()

        local_request = root / 'request.json'
        session_keys = {'userid', 'pageid', 'idnumber', 'numhits', 'unique'}
        if local_request.exists():
            saved = json.loads(local_request.read_text())
            if {k: v for k, v in saved.items() if k not in session_keys} != {k: v for k, v in q.items() if k not in session_keys}:
                raise ValueError('saved scientific request differs')
            q = saved
        else:
            claims = Path('/tmp/ember-tops-form-id-claims-v1')
            if not (claims / 'inventory.json').is_file():
                raise ValueError('initialize the local identifier inventory before fetching')
            # Simultaneous GETs returned the same identifier in a saved control.
            # Serialize allocation only; submitted calculations can run in parallel.
            with (claims / 'allocation.lock').open('a') as lock:
                fcntl.flock(lock, fcntl.LOCK_EX)
                claimed = False
                for attempt in range(4):
                    form_path = root / f'form-{attempt:02d}.html'
                    if form_path.exists():
                        raise ValueError('incomplete form allocation needs review before another attempt')
                    url = 'https://aphysics2.lanl.gov/opacity/lanl/?ember_request=' + uuid.uuid4().hex
                    form_path.write_text(fetch('form', url))
                    form = SubmissionForm()
                    form.feed(form_path.read_text())
                    if set(form.fields) != session_keys or not re.fullmatch('[a-z]{3}', form.fields['userid']):
                        raise ValueError('unexpected official session fields')
                    claim_path = claims / (form.fields['userid'] + '.json')
                    try:
                        with claim_path.open('x') as stream:
                            json.dump({'owner': 'primary', 'pid': os.getpid(),
                                       'source_root': str(root.resolve()),
                                       'recipe_sha256': digest(recipe_path),
                                       'official_form_sha256': digest(form_path)}, stream, indent=2)
                            stream.write('\n')
                    except FileExistsError:
                        continue
                    claimed = True
                    (root / 'identifier_claim.json').write_bytes(claim_path.read_bytes())
                    q.update(form.fields)
                    local_request.write_text(json.dumps(q, indent=2) + '\n')
                    break
                if not claimed:
                    raise ValueError('official form repeatedly supplied an already claimed identifier')
        receipt_path, source_path = root / 'receipt.json', root / 'source.txt'
        if receipt_path.exists():
            receipt = json.loads(receipt_path.read_text())
            if receipt['request_sha256'] != digest(local_request):
                raise ValueError('completed request changed')
        else:
            submitted = root / 'submit.html'
            if not submitted.exists():
                submitted.write_text(fetch('submit', 'https://aphysics2.lanl.gov/submit', q))
            handle = Form()
            handle.feed(submitted.read_text())
            if not handle.done or handle.data.get('output') != 'tabcol':
                raise ValueError('unexpected source handle')
            result_path = root / 'results.html'
            if not result_path.exists():
                result_path.write_text(fetch('results', 'https://aphysics2.lanl.gov/results', handle.data))
            result = Text()
            result.feed(result_path.read_text())
            source = '\n'.join(v.strip() for v in ''.join(result.parts).replace('\u00a0', ' ').splitlines() if v.strip()) + '\n'
            if source_path.exists() and source_path.read_text() != source:
                raise ValueError('retained source differs from response')
            if not source_path.exists():
                source_path.write_text(source)
            receipt = {**old, 'file': source_path.name, 'sha256': digest(source_path),
                       'request': local_request.name, 'request_sha256': digest(local_request),
                       'dimensions': [len(tt), len(rr)],
                       'plasma_cutoff': 'on' if args.kind == 'spectrum' else 'off',
                       'validation_role': 'official_form_spectrum_control' if args.kind == 'spectrum' else 'official_form_uncut_' + args.kind}
        actual_t, actual_r, cells, excluded = read(source_path, receipt, dimensions=(len(tt), len(rr)))
        if actual_t != tt or any(abs(x / y - 1) > 5e-5 for x, y in zip(actual_r, rr, strict=True)):
            raise ValueError('returned coordinates differ from plan')
        if args.kind in ('groups', 'spectrum') and excluded:
            raise ValueError('substituted states cannot supply requested groups')
        if not receipt_path.exists():
            receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
        if args.kind == 'groups':
            boundaries, groups, _ = source_groups(root)
            if len(boundaries) != n or len(groups) != count:
                raise ValueError('incomplete source groups')
        if args.kind == 'spectrum':
            spectra, _ = spectral_source(root)
            if set(spectra) != set(cells) or len(spectra) != count:
                raise ValueError('incomplete source spectra')
        print(json.dumps({'request_index': index, 'kind': args.kind, 'work': str(root),
                          'states': len(cells), 'source_exclusions': len(excluded)}), flush=True)


if __name__ == '__main__':
    main()
