#!/usr/bin/env python3
"""Check that an atmosphere temperature extension retains every source and mask.

This compares the stored matching states, including missing-state flags. It
does not compare derivatives at the old upper temperature boundary, where
the adjacent interpolation interval changes when the grid is extended.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import shlex

from prepare_nongrey_sources import digest


def read(path):
    lines = path.read_text().splitlines()
    if lines[0] != 'EMBER_COMPOSITION_ATMOSPHERE 2':
        raise ValueError('expected the masked atmosphere format')
    header, axes = {}, []
    for line in lines[1:]:
        fields = shlex.split(line)
        if fields == ['data']:
            break
        key = fields[0]
        if key in header:
            raise ValueError('repeated atmosphere metadata')
        header[key] = fields[1:]
    if set(header) != {'source', 'approximation', 'basis', 'tau', 'metals',
                       'hydrogen', 'helium3', 'log_teff', 'log_g'}:
        raise ValueError('unexpected atmosphere metadata')
    for key in ['hydrogen', 'helium3', 'log_teff', 'log_g']:
        fields = header.pop(key)
        axis = list(map(float, fields[1:]))
        if (len(axis) != int(fields[0]) or len(axis) < 2
                or not all(math.isfinite(v) for v in axis)
                or any(a >= b for a, b in zip(axis, axis[1:]))):
            raise ValueError('invalid atmosphere axis')
        axes.append(axis)
    rows = lines[lines.index('data')+1:]
    states = {}
    for key, line in zip(itertools.product(*axes), rows, strict=True):
        fields = line.split()
        if fields == ['0']:
            states[key] = (0,)
        elif len(fields) == 3 and fields[0] == '1':
            values = tuple(map(float, fields[1:]))
            if not all(math.isfinite(v) for v in values):
                raise ValueError('invalid atmosphere matching state')
            states[key] = (1, *values)
        else:
            raise ValueError('invalid atmosphere mask or matching state')
    return header, axes, states


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['previous', 'candidate', 'output']:
        p.add_argument(name, type=Path)
    a = p.parse_args()
    inputs = {str(f.resolve()): digest(f) for f in [a.previous, a.candidate, Path(__file__)]}
    old_header, old_axes, old = read(a.previous)
    new_header, new_axes, new = read(a.candidate)
    if old_header != new_header:
        raise ValueError('atmosphere source prescription changed')
    if (any(old_axes[i] != new_axes[i] for i in [0, 1, 3])
            or new_axes[2][:len(old_axes[2])] != old_axes[2]
            or len(new_axes[2]) <= len(old_axes[2])):
        raise ValueError('expected only an upper-temperature extension')
    if any(new.get(key) != value for key, value in old.items()):
        raise ValueError('an existing source state or missing-state mask changed')
    if any(digest(f) != checksum for f, checksum in inputs.items()):
        raise ValueError('an input changed during the comparison')
    report = {'scope': __doc__, 'passed': True, 'input_sha256': inputs,
              'retained_source_states': sum(v[0] for v in old.values()),
              'retained_missing_states': sum(v[0] == 0 for v in old.values()),
              'candidate_source_states': sum(v[0] for v in new.values()),
              'source_values_and_masks_exact': True}
    a.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ['scope', 'input_sha256']}))


if __name__ == '__main__':
    main()
