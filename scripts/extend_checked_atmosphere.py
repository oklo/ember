#!/usr/bin/env python3
"""Fill missing nodes in an already checked atmosphere table.

The existing table is pinned by a completed stellar-run receipt and a primary
source review, and its manifest by a separate retained input record. Existing rows and axes remain
byte for byte unchanged. New profiles undergo full source revalidation.
Old source recipes remain linked through the parent manifest; this operation
does not claim to replay old spectral binaries that are no longer retained.
"""
import argparse
import itertools
import json
import math
from pathlib import Path

from assemble_nongrey_grid import complete_cells, load_continuation
from atmosphere_capacity_identity import CapacityEquivalence
from audit_atmosphere_retention import read
from prepare_nongrey_sources import digest


def pinned(record, path):
    """A path alias is allowed only when the record names the same file."""
    path = path.resolve(strict=True)
    matches = [checksum for name, checksum in record.get('input_sha256', {}).items()
               if Path(name).resolve() == path]
    if not matches or any(checksum != digest(path) for checksum in matches):
        raise ValueError('input is not pinned by the supplied record: '+str(path))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['previous_table', 'previous_manifest', 'stellar_receipt', 'source_review',
                 'manifest_pin', 'continuation_plan', 'output']:
        p.add_argument(name, type=Path)
    a = p.parse_args()
    output_manifest = a.output.with_suffix('.manifest.json')
    if a.output.exists() or output_manifest.exists():
        raise FileExistsError('preserve the existing candidate')
    parent = json.loads(a.previous_manifest.read_text())
    stellar = json.loads(a.stellar_receipt.read_text())
    # The CN driver's acceptance flag deliberately stays false because its
    # complete physical trajectory is unfinished. It is not a table-status flag.
    if (stellar.get('native_exit_code') != 0 or stellar.get('failure') is not None
            or not stellar.get('history')
            or stellar.get('outcome') not in ['requested_checkpoint_stop', 'completed',
                                             'requested_interval_reached',
                                             'surface_hydrogen_review_trigger']):
        raise ValueError('a clean completed stellar-run receipt is required')
    source_review = json.loads(a.source_review.read_text())
    if (source_review.get('outcome') != 'passed_raw_sources_and_retention'
            or not any(Path(name).suffix == '.dat' and checksum == digest(a.previous_table)
                       for name, checksum in source_review.get('sha256', {}).items())):
        raise ValueError('the parent table lacks the required primary source review')
    pinned(stellar, a.previous_table)
    pinned(json.loads(a.manifest_pin.read_text()), a.previous_manifest)
    if digest(a.previous_table) != parent['table_sha256']:
        raise ValueError('parent table and manifest disagree')
    _, raw_axes, before = read(a.previous_table)
    axes = parent['axes']
    if len(axes) != 4:
        raise ValueError('unexpected parent axes')
    for i, (raw, physical) in enumerate(zip(raw_axes, axes, strict=True)):
        expected = [math.log10(v) if i == 2 else v for v in physical]
        if raw != expected:
            raise ValueError('table and source axes differ')
    if sum(row[0] for row in before.values()) != parent['accepted_states']:
        raise ValueError('parent state count differs')
    paths = json.loads(a.continuation_plan.read_text())['continuations']
    if not paths or len(set(paths)) != len(paths):
        raise ValueError('empty or duplicate continuation list')
    dependencies = {str(path.resolve()): digest(path) for path in
                    [a.previous_table, a.previous_manifest, a.stellar_receipt, a.source_review,
                     a.manifest_pin, a.continuation_plan, Path(__file__)]}
    capacity = CapacityEquivalence(parent.get('capacity_replays', []))
    dependencies.update(capacity.inputs)
    lines = a.previous_table.read_text().splitlines()
    start = lines.index('data')+1
    keys = list(itertools.product(*raw_axes))
    index = {key: i for i, key in enumerate(keys)}
    added = []
    for work in paths:
        key, state, spec, record, source_inputs = load_continuation(
            Path(work), parent['physics'], capacity)
        if any(value not in axis for value, axis in zip(key, axes, strict=True)):
            raise ValueError('this operation fills existing axes only')
        raw_key = (key[0], key[1], math.log10(key[2]), key[3])
        location = start+index[raw_key]
        if lines[location] != '0':
            raise ValueError('cannot replace an existing source state')
        lines[location] = f'1 {math.log10(state["T"]):.17g} {math.log10(state["Pgas"]):.17g}'
        dependencies.update(source_inputs)
        added.append(dict(coordinates=key, state=state, work=str(Path(work).resolve()),
                          source_specification=spec,
                          validation_sha256=digest(Path(work)/'final/validated.json')))
    if any(digest(path) != checksum for path, checksum in dependencies.items()):
        raise ValueError('an input changed during assembly')
    physical_keys = list(itertools.product(*axes))
    present = {key for key, row in zip(physical_keys, lines[start:], strict=True) if row != '0'}
    with a.output.open('x') as f:
        f.write('\n'.join(lines)+'\n')
    _, after_axes, after = read(a.output)
    if after_axes != raw_axes or any(after[key] != value for key, value in before.items() if value[0]):
        raise ValueError('existing table values changed')
    expected_added = {(k['coordinates'][0], k['coordinates'][1],
                       math.log10(k['coordinates'][2]), k['coordinates'][3]) for k in added}
    actual_added = {key for key in before if not before[key][0] and after[key][0]}
    if actual_added != expected_added or len(present) != parent['accepted_states']+len(added):
        raise ValueError('unexpected change in missing states')
    report = dict(scope=__doc__, assembly_method='fill_missing_nodes_in_checked_table',
                  base_table=str(a.previous_table.resolve()),
                  parent_source_manifest=str(a.previous_manifest.resolve()),
                  inherited_source_states=parent['accepted_states'], accepted_states=len(present),
                  missing_states=len(before)-len(present), axes=axes,
                  axis_order=parent['axis_order'], physics=parent['physics'],
                  complete_cells=complete_cells(axes, present), extension=added,
                  capacity_replays=parent.get('capacity_replays', []),
                  input_files_sha256=dependencies, table_sha256=digest(a.output),
                  retained_rows_exact=True, new_states=len(added),
                  selected_for_stellar_evolution=False)
    with output_manifest.open('x') as f:
        json.dump(report, f, indent=2, allow_nan=False); f.write('\n')
    print(json.dumps({k: report[k] for k in ['inherited_source_states', 'accepted_states',
                                          'new_states', 'table_sha256']}))


if __name__ == '__main__':
    main()
