"""Identify and verify an optional tabulated bulk EOS used by an atmosphere."""
import hashlib
import math
import re
from pathlib import Path


def bulk_identity(prepared):
    record = prepared.get('bulk_eos')
    if record is None:
        if prepared.get('bulk_eos_experiment', {}).get('enabled'):
            raise ValueError('bulk EOS experiment has no reviewed source identity')
        return None
    family = record.get('family', {})
    if (record.get('format') != 1 or not family.get('model')
            or not family.get('source_sha256') or not family.get('recipe_sha256')):
        raise ValueError('incomplete bulk EOS source identity')
    hashes = [*family['source_sha256'].values(), *family['recipe_sha256'].values(),
              record.get('table_sha256', '')]
    if any(not isinstance(v, str) or not re.fullmatch('[0-9a-f]{64}', v) for v in hashes):
        raise ValueError('invalid bulk EOS source checksum')
    return family


def validate_bulk_eos(directory, spec, prepared, log):
    family = bulk_identity(prepared)
    marker = 'EMBER NONIDEAL BULK EOS:'
    if family is None:
        if marker in log:
            raise ValueError('atmosphere used an unrecorded bulk EOS')
        return None
    record = prepared['bulk_eos']
    path = Path(directory)/record['table']
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != record['table_sha256']:
        raise ValueError('bulk EOS table checksum mismatch')
    lines = raw.decode().splitlines()
    if not lines or lines[0] not in ('EMBER_ATMOSPHERE_BULK_EOS 1',
                                       'EMBER_ATMOSPHERE_BULK_EOS 2'):
        raise ValueError('unsupported atmosphere bulk EOS table')
    nt, np, x = lines[1].split()
    nt, np, x = int(nt), int(np), float(x)
    version = int(lines[0].rsplit(maxsplit=1)[1])
    end = 4+nt*np
    if nt < 2 or np < 2 or len(lines) < end or (version == 1 and len(lines) != end):
        raise ValueError('invalid bulk EOS table dimensions')
    if version == 2:
        if len(lines) <= end+1 or lines[end] != 'cells':
            raise ValueError('missing bulk EOS cell mask')
        mask = ' '.join(lines[end+1:]).split()
        if len(mask) != (nt-1)*(np-1) or any(v not in ('0', '1') for v in mask):
            raise ValueError('invalid bulk EOS cell mask')
    for line, n in zip(lines[2:4], [nt, np]):
        a = list(map(float, line.split()))
        if len(a) != n or any(not math.isfinite(v) for v in a) or any(v >= w for v, w in zip(a, a[1:])):
            raise ValueError('invalid bulk EOS table axis')
    for line in lines[4:end]:
        a = list(map(float, line.split()))
        if len(a) != 12 or any(not math.isfinite(v) for v in a):
            raise ValueError('invalid bulk EOS table entry')
    if (not 0 < x < 1 or spec['hydrogen'] != [x] or spec['helium3'] != [0.]
            or sum(spec['metals']) > 1e-8):
        raise ValueError('bulk EOS table composition mismatch')
    loaded = [line.rsplit(maxsplit=1)[-1] for line in log.splitlines() if marker in line]
    values = [float(v.replace('D', 'E')) for v in loaded]
    if not values or any(not math.isfinite(v) or abs(v-x) > 1e-8 for v in values):
        raise ValueError('source log does not confirm the declared bulk EOS')
    return dict(table_sha256=record['table_sha256'], hydrogen=x, family=family)
