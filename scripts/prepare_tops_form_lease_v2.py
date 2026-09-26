#!/usr/bin/env python3
"""Prepare one source request with an official identifier safe to reuse locally.

An older identifier is reusable only after its uncut ATOMIC response is saved,
checksum-verified, and scientifically validated. Its allocation record is
archived before replacement. Active and unverifiable requests remain reserved.
The frozen source fetcher independently verifies the prepared scientific fields.
Malformed or incomplete previous group responses remain reserved; the allocator
skips their identifiers instead of aborting unrelated new requests.
"""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import time
import urllib.request
import uuid

from fetch_tops_composition import digest
from fetch_tops_form_plan_v2 import SubmissionForm
from import_tops_composition import read
from tops_groups import source_groups


def completed_source(claim, userid):
    if claim.get('owner') != 'primary' or 'source_root' not in claim:
        return None
    root = Path(claim['source_root'])
    receipt_path, request_path = root / 'receipt.json', root / 'request.json'
    if not receipt_path.is_file() or not request_path.is_file():
        return None
    receipt, request = json.loads(receipt_path.read_text()), json.loads(request_path.read_text())
    if (request.get('userid') != userid or request.get('lib') != 'new' or
            request.get('plasnu') != 'off' or digest(request_path) != receipt['request_sha256']):
        return None
    source_path = root / receipt['file']
    read(source_path, receipt, dimensions=tuple(receipt['dimensions']))
    if request['datype'] == 'groups':
        source_groups(root)
    elif request['datype'] != 'gray':
        return None
    return {str(p.resolve()): digest(p) for p in (receipt_path, request_path, source_path)}


def checked_completed_source(claim, userid):
    """Keep malformed, missing or scientifically incomplete prior sources reserved."""
    try:
        return completed_source(claim, userid), None
    except (OSError, ValueError, KeyError, IndexError) as error:
        return None, type(error).__name__ + ': ' + str(error)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kind', choices=('grey', 'groups'))
    parser.add_argument('plan', type=Path)
    parser.add_argument('request_index', type=int)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    if not 0 <= args.request_index < len(plan['requests']):
        raise ValueError('invalid request index')
    job = plan['requests'][args.request_index]
    root = Path(job['work'])
    root.mkdir(parents=True, exist_ok=True)
    request_path = root / 'request.json'
    if request_path.exists():
        print(json.dumps({'request_index': args.request_index, 'already_prepared': True}), flush=True)
        return
    if args.kind == 'grey':
        manifest_path = Path(plan['baseline_manifest'])
        if digest(manifest_path) != plan['baseline_manifest_sha256']:
            raise ValueError('baseline manifest changed')
        baseline = next(p for p in json.loads(manifest_path.read_text())['planes']
                        if (p['X'], p['Z']) == (job['X'], job['Z']))
        seed_request = manifest_path.parent / baseline['request']
    else:
        seed = Path(plan['baseline'])
        baseline = json.loads((seed / 'receipt.json').read_text())
        seed_request = seed / baseline['request']
    if digest(seed_request) != baseline['request_sha256']:
        raise ValueError('baseline request changed')
    tt, rr = job['temperatures_keV'], job['densities_g_cm3']
    q = {**json.loads(seed_request.read_text()), 'datype': 'gray' if args.kind == 'grey' else 'groups',
         'plasnu': 'off', 'tgrid': 'specific', 'temps': ' '.join(map(str, tt)),
         'rgrid': 'range', 'rlow': str(rr[0]), 'rup': str(rr[-1]), 'nr': str(len(rr)),
         'rspace': 'log' if args.kind == 'grey' else job.get('density_spacing', 'linear')}
    if args.kind == 'groups':
        q.update(egrid='range', egplow=str(plan['photon_min_keV']), egphigh=str(plan['photon_max_keV']),
                 ngpengs=str(plan['photon_boundaries']), espace='log')
    claims = Path('/tmp/ember-tops-form-id-claims-v1')
    if not (claims / 'inventory.json').is_file():
        raise ValueError('identifier inventory is missing')
    inputs = {str(p.resolve()): digest(p) for p in
              (args.plan, seed_request, Path(__file__), Path('scripts/fetch_tops_form_plan_v2.py'),
               Path('scripts/import_tops_composition.py'), Path('scripts/tops_groups.py'))}
    attempts = []
    with (claims / 'allocation.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        for attempt in range(32):
            form_path = root / f'lease-form-{attempt:02d}.html'
            if form_path.exists():
                raise ValueError('incomplete lease allocation requires review')
            url = 'https://aphysics2.lanl.gov/opacity/lanl/?ember_request=' + uuid.uuid4().hex
            start = time.monotonic()
            with urllib.request.urlopen(urllib.request.Request(url, headers={'Cache-Control': 'no-cache'}), timeout=75) as response:
                raw = response.read(1000001)
            if len(raw) > 1000000:
                raise ValueError('official form exceeds the data cap')
            form_path.write_bytes(raw)
            form = SubmissionForm()
            form.feed(raw.decode())
            if set(form.fields) != {'userid', 'pageid', 'idnumber', 'numhits', 'unique'} or not re.fullmatch('[a-z]{3}', form.fields['userid']):
                raise ValueError('unexpected official session fields')
            userid = form.fields['userid']
            claim_path = claims / (userid + '.json')
            previous = None
            attempts.append({'userid': userid, 'elapsed_seconds': time.monotonic() - start,
                             'form_file': form_path.name, 'form_sha256': digest(form_path)})
            if claim_path.exists():
                old_bytes = claim_path.read_bytes()
                previous, rejection = checked_completed_source(json.loads(old_bytes), userid)
                if rejection is not None:
                    attempts[-1]['unverifiable_prior_source'] = rejection
                if previous is None:
                    attempts[-1]['reserved_or_unverifiable'] = True
                    continue
                history = claims / 'history'
                history.mkdir(exist_ok=True)
                archived = history / (userid + '-' + hashlib.sha256(old_bytes).hexdigest() + '.json')
                if archived.exists():
                    if archived.read_bytes() != old_bytes:
                        raise ValueError('allocation history differs')
                else:
                    archived.write_bytes(old_bytes)
                inputs.update(previous)
                inputs[str(archived.resolve())] = digest(archived)
            q.update(form.fields)
            allocation = {'owner': 'primary', 'pid': os.getpid(), 'source_root': str(root.resolve()),
                          'plan_sha256': digest(args.plan), 'request_index': args.request_index,
                          'official_form_sha256': digest(form_path),
                          'previous_completed_source_sha256': previous,
                          'allocator_sha256': digest(Path(__file__))}
            tmp = claim_path.with_suffix('.json.next')
            tmp.write_text(json.dumps(allocation, indent=2) + '\n')
            os.replace(tmp, claim_path)
            (root / 'identifier_lease.json').write_bytes(claim_path.read_bytes())
            request_path.write_text(json.dumps(q, indent=2) + '\n')
            receipt = {'scope': __doc__, 'request_index': args.request_index, 'userid': userid,
                       'reused_completed_identifier': previous is not None,
                       'request_sha256': digest(request_path), 'attempts': attempts, 'input_sha256': inputs}
            (root / 'lease_allocation.json').write_text(json.dumps(receipt, indent=2) + '\n')
            print(json.dumps({k: v for k, v in receipt.items() if k not in ('scope', 'attempts', 'input_sha256')}), flush=True)
            return
    raise ValueError('no reusable official identifier found within the bounded attempts')


if __name__ == '__main__':
    main()
