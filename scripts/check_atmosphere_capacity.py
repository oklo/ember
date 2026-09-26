#!/usr/bin/env python3
"""Replay an existing atmosphere with a larger source array capacity.

This checks the actual input decks and full saved numerical outputs. It does
not authorize combining executable identities in an atmosphere table.
"""
import argparse
import gzip
import json
import math
from pathlib import Path
import re

from assemble_nongrey_grid import load_continuation
from generate_nongrey_grid import execute, input_fingerprint, sequence, temperatures
from import_nongrey_grid import source_state
from prepare_nongrey_sources import digest


def contents(path):
    return path.read_bytes() if path.exists() else gzip.decompress(
        path.with_name(path.name + '.gz').read_bytes())


def phase_timings(path):
    rows = []
    previous = 0.
    for line in contents(path).decode().splitlines():
        iteration, stage, cumulative, seconds, label = line.split(maxsplit=4)
        cumulative, seconds = float(cumulative), float(seconds)
        if (not all(math.isfinite(v) and v >= 0 for v in [cumulative, seconds])
                or cumulative < previous or abs(cumulative-previous-seconds) > .02001):
            raise ValueError('invalid printed CPU timing')
        rows.append({'iteration': int(iteration), 'stage': int(stage), 'label': label,
                     'cumulative_cpu_seconds': cumulative, 'stage_cpu_seconds': seconds})
        previous = cumulative
    if not rows:raise ValueError('missing phase timings')
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('prepared', type=Path)
    p.add_argument('model', type=Path, help='completed continuation directory')
    p.add_argument('work', type=Path)
    p.add_argument('report', type=Path)
    p.add_argument('--completed-replay', action='store_true',
                   help='check existing completed outputs; never launch or rewrite their inputs')
    a = p.parse_args()
    prepared = json.loads(a.prepared.read_text())
    variant = prepared['array_capacity']
    base = variant['canonical_prepared']
    original = Path(base['tlusty']).parent
    enlarged = Path(prepared['tlusty']).parent
    if variant['parameter'] not in ['MDEPTH', 'MTABT']:
        raise ValueError('unrecognized capacity parameter')
    if not variant['value'] > variant['original']:
        raise ValueError('capacity did not increase')
    for source in [base, prepared]:
        if digest(source['tlusty']) != source['executables']['tlusty']:
            raise ValueError('executable changed')
    if digest(original/'BASICS.FOR') != variant['original_parameter_file_sha256']:
        raise ValueError('original parameter file changed')
    pattern = r'(\b' + variant['parameter'] + r'\s*=\s*)' + str(variant['original']) + r'\b'
    expected, count = re.subn(pattern, lambda m: m[1] + str(variant['value']),
                              (original/'BASICS.FOR').read_text())
    if count != 1:
        raise ValueError('unrecognized original capacity')
    for name, change in variant.get('additional_capacities', {}).items():
        if name != 'MFREQ' or not 2 <= change['value'] < change['original']:
            raise ValueError('unrecognized additional capacity change')
        expected, count = re.subn(r'(\b'+name+r'\s*=\s*)'+str(change['original'])+r'\b',
                                  lambda m:m[1]+str(change['value']),expected)
        if count != 1:raise ValueError('unrecognized original frequency capacity')
    if (enlarged/'BASICS.FOR').read_text() != expected:
        raise ValueError('source differs beyond the recorded array capacities')
    for name, expected_hash in variant['unchanged_source_sha256'].items():
        if digest(original/name) != expected_hash or digest(enlarged/name) != expected_hash:
            raise ValueError('equation source changed: ' + name)
    key, _, spec, record, dependencies = load_continuation(a.model)
    if 'MFREQ' in variant.get('additional_capacities', {}):
        if spec['atmosphere_frequencies'] > variant['additional_capacities']['MFREQ']['value']:
            raise ValueError('replay exceeds the frequency capacity')
    if record['continuation_provenance']['prepared'] != base:
        raise ValueError('model uses another original executable or source')
    old = a.model/'final'
    if not a.completed_replay:
        a.work.mkdir(parents=True, exist_ok=False)
        for name in ['fort.5', 'fort.15', 'tas', 'ember-masses.dat',
                     'opacity.sha256', 'physics.json', 'fort.8']:
            (a.work/name).write_bytes(contents(old/name))
        for name in ['fort.2', 'fort.55', 'ember-condensates.cfg',
                     'condensate-abundances.dat', 'condensates.sha256']:
            if (old/name).exists() and (old/name).stat().st_size:
                (a.work/name).write_bytes((old/name).read_bytes())
        for name in ['data', 'opacity.bin']:
            (a.work/name).symlink_to((old/name).resolve(), target_is_directory=name == 'data')
    saved = json.loads((old/'completed.json').read_text())
    if input_fingerprint(base['executables']['tlusty'], a.work) != saved['input_sha256']:
        raise ValueError('copied inputs differ from the completed model')
    if a.completed_replay:
        complete = json.loads((a.work/'completed.json').read_text())
        if (complete['input_sha256'] != input_fingerprint(prepared['executables']['tlusty'], a.work)
                or not {'fort.7', 'fort.9', 'run.log'}.issubset(complete['outputs'])
                or any(digest(a.work/name) != checksum for name, checksum in complete['outputs'].items())):
            raise ValueError('existing replay outputs or input fingerprint changed')
    else:
        execute(prepared['tlusty'], a.work, ['fort.7', 'fort.9'])
    state = source_state((a.work/'run.log').read_text(), (a.work/'fort.9').read_text(),
                         key[2], key[3], {'temperature_K': temperatures(spec),
                         'density_g_cm3': sequence(spec['log_density'])}, spec['tau'])
    comparisons = {}
    for name in ['fort.7', 'fort.9', 'fort.10', 'fort.11', 'fort.13',
                 'fort.14', 'fort.17', 'fort.18']:
        comparisons[name] = {'original_sha256': digest(old/name),
                             'replay_sha256': digest(a.work/name),
                             'exact': contents(old/name) == contents(a.work/name)}
    diagnostics_equal = state == record['diagnostics']
    timings = {'original': phase_timings(old/'fort.69'), 'replay': phase_timings(a.work/'fort.69')}
    signatures = [[(r['iteration'], r['stage'], r['label']) for r in rows] for rows in timings.values()]
    phases_equal = signatures[0] == signatures[1]
    passed = all(row['exact'] for row in comparisons.values()) and diagnostics_equal and phases_equal
    dependencies.update({str(a.prepared.resolve()): digest(a.prepared),
                         str(Path(__file__).resolve()): digest(Path(__file__)),
                         str(Path(prepared['tlusty']).resolve()): digest(prepared['tlusty'])})
    report = {'passed': passed, 'scope': __doc__, 'input_sha256': dependencies,
              'parameter': variant['parameter'], 'original_capacity': variant['original'],
              'new_capacity': variant['value'], 'coordinates': key,
              'additional_capacities': variant.get('additional_capacities', {}),
              'depths': spec['depths'], 'frequencies': spec['atmosphere_frequencies'],
              'original_executable_sha256': base['executables']['tlusty'],
              'new_executable_sha256': prepared['executables']['tlusty'],
              'outputs': comparisons, 'diagnostics_exact': diagnostics_equal,
              'timed_phase_sequence_exact': phases_equal, 'phase_timings': timings,
              'timing_note': 'fort.69 contains CPU timings, not physical outputs. Compare its iteration/stage labels exactly; elapsed CPU values are recorded separately and may differ. These phases exclude final diagnostic output processing.',
              'timing_files_sha256': {str(v.resolve()): digest(v) for v in [old/'fort.69', a.work/'fort.69']},
              'diagnostics': state, 'work': str(a.work.resolve()),
              'seconds': json.loads((a.work/'completed.json').read_text())['seconds']}
    a.report.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': passed, 'report': str(a.report),
                      'diagnostics_exact': diagnostics_equal}), flush=True)
    if not passed:
        raise ValueError('capacity replay differs; inspect the complete report')


if __name__ == '__main__':
    main()
