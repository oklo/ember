#!/usr/bin/env python3
"""Append higher-hydrogen atmosphere planes to a checked source table.

Existing source rows and missing-state flags are retained literally. Other
axes and physical prescriptions stay fixed. New source profiles are fully
revalidated; empty corners remain unsupported. This does not select the table.
"""
import argparse
import itertools
import json
import math
from pathlib import Path

from assemble_nongrey_grid import complete_cells, load_continuation
from atmosphere_capacity_identity import CapacityEquivalence
from audit_atmosphere_retention import read
from extend_checked_atmosphere import pinned
from prepare_nongrey_sources import digest


def append_hydrogen_rows(lines, raw_axes, states):
    """Append validated (physical coordinates, matching state) pairs."""
    if not states:
        raise ValueError('no new source states')
    additions = {}
    for key, state in states:
        if len(key) != 4 or not all(math.isfinite(x) for x in key):
            raise ValueError('invalid source coordinates')
        if not raw_axes[0][-1] < key[0] < 1 or key[2] <= 0:
            raise ValueError('only higher hydrogen coordinates may be added')
        raw_key = (key[0], key[1], math.log10(key[2]), key[3])
        if any(raw_key[i] not in raw_axes[i] for i in [1, 2, 3]):
            raise ValueError('other source axes must remain unchanged')
        if raw_key in additions:
            raise ValueError('duplicate source coordinate')
        if not all(math.isfinite(state[k]) and state[k] > 0 for k in ['T', 'Pgas']):
            raise ValueError('invalid matching state')
        additions[raw_key] = f'1 {math.log10(state["T"]):.17g} {math.log10(state["Pgas"]):.17g}'
    extra = sorted({key[0] for key in additions})
    axes = [raw_axes[0]+extra, *raw_axes[1:]]
    header_end = lines.index('data')+1
    if len(lines)-header_end != math.prod(map(len, raw_axes)):
        raise ValueError('parent table size differs')
    output = lines.copy()
    hydrogen_line = next(i for i, line in enumerate(lines[:header_end])
                         if line.startswith('hydrogen '))
    fields = lines[hydrogen_line].split()
    # Preserve existing coordinate text, including its round-trip precision.
    output[hydrogen_line] = ' '.join(['hydrogen', str(len(axes[0])), *fields[2:],
                                     *(format(x, '.17g') for x in extra)])
    output.extend(additions.get(key, '0')
                  for key in itertools.product(extra, *raw_axes[1:]))
    return output, axes, set(additions)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['previous_table', 'previous_manifest', 'stellar_receipt',
                 'source_review', 'manifest_pin', 'continuation_plan', 'output']:
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    manifest_path = args.output.with_suffix('.manifest.json')
    if args.output.exists() or manifest_path.exists():
        raise FileExistsError('preserve existing candidates')
    parent = json.loads(args.previous_manifest.read_text())
    stellar = json.loads(args.stellar_receipt.read_text())
    if (stellar.get('native_exit_code') != 0 or stellar.get('failure') is not None
            or not stellar.get('history') or stellar.get('outcome') not in
            ['requested_checkpoint_stop', 'completed', 'requested_interval_reached',
             'surface_hydrogen_review_trigger']):
        raise ValueError('a clean completed stellar receipt is required')
    review = json.loads(args.source_review.read_text())
    if (review.get('outcome') != 'passed_raw_sources_and_retention'
            or not any(Path(p).suffix == '.dat' and v == digest(args.previous_table)
                       for p, v in review.get('sha256', {}).items())):
        raise ValueError('parent source review missing or mismatched')
    pinned(stellar, args.previous_table)
    pinned(json.loads(args.manifest_pin.read_text()), args.previous_manifest)
    if digest(args.previous_table) != parent['table_sha256']:
        raise ValueError('parent table and manifest disagree')
    header, raw_axes, before = read(args.previous_table)
    axes = parent['axes']
    if len(axes) != 4 or any(raw != [math.log10(x) if i == 2 else x for x in axis]
                            for i, (raw, axis) in enumerate(zip(raw_axes, axes, strict=True))):
        raise ValueError('parent physical and stored axes disagree')
    if sum(row[0] for row in before.values()) != parent['accepted_states']:
        raise ValueError('parent state count differs')
    paths = json.loads(args.continuation_plan.read_text())['continuations']
    if not paths or len(paths) != len(set(paths)):
        raise ValueError('empty or duplicate continuation list')
    dependencies = {str(p.resolve()): digest(p) for p in
                    [args.previous_table, args.previous_manifest, args.stellar_receipt,
                     args.source_review, args.manifest_pin, args.continuation_plan,
                     Path(__file__), Path(__file__).with_name('extend_checked_atmosphere.py')]}
    capacity = CapacityEquivalence(parent.get('capacity_replays', []))
    dependencies.update(capacity.inputs)
    added = []
    for work in paths:
        key, state, spec, record, source_inputs = load_continuation(
            Path(work), parent['physics'], capacity)
        dependencies.update(source_inputs)
        added.append(dict(coordinates=key, state=state, work=str(Path(work).resolve()),
                          source_specification=spec,
                          validation_sha256=digest(Path(work)/'final/validated.json')))
    lines = args.previous_table.read_text().splitlines()
    output, new_raw_axes, expected_new = append_hydrogen_rows(
        lines, raw_axes, [(x['coordinates'], x['state']) for x in added])
    new_axes = [new_raw_axes[0], *axes[1:]]
    if any(digest(p) != v for p, v in dependencies.items()):
        raise ValueError('an input changed during assembly')
    with args.output.open('x') as stream:
        stream.write('\n'.join(output)+'\n')
    new_header, actual_axes, after = read(args.output)
    if (new_header != header or actual_axes != new_raw_axes
            or any(after[key] != value for key, value in before.items())):
        raise ValueError('existing metadata, source values or missing flags changed')
    actual_new = {key for key, value in after.items() if value[0] and key not in before}
    if actual_new != expected_new:
        raise ValueError('new source coordinates differ')
    start = output.index('data')+1
    present = {key for key, row in zip(itertools.product(*new_axes), output[start:], strict=True)
               if row != '0'}
    report = dict(scope=__doc__, assembly_method='append_higher_hydrogen_planes',
                  base_table=str(args.previous_table.resolve()),
                  parent_source_manifest=str(args.previous_manifest.resolve()),
                  inherited_source_states=parent['accepted_states'], accepted_states=len(present),
                  missing_states=len(after)-len(present), axes=new_axes,
                  axis_order=parent['axis_order'], physics=parent['physics'],
                  complete_cells=complete_cells(new_axes, present), extension=added,
                  capacity_replays=parent.get('capacity_replays', []),
                  input_files_sha256=dependencies, table_sha256=digest(args.output),
                  retained_rows_exact=True, retained_missing_masks_exact=True,
                  new_states=len(added), selected_for_stellar_evolution=False)
    with manifest_path.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps({k: report[k] for k in
                      ['inherited_source_states', 'accepted_states', 'new_states', 'table_sha256']}))


if __name__ == '__main__':
    main()
