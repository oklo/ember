#!/usr/bin/env python3
"""Check hotter opacity rows while retaining every previous numeric cell.

This checks finite source tables, unchanged composition, wavelengths and density
coordinates, and exact retention of all earlier temperatures, electron densities
and absorption coefficients. It does not establish accuracy between source rows.
"""
import argparse
import json
from pathlib import Path

from nongrey_opacity import read_table
from prepare_nongrey_sources import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--pair', nargs=2, action='append', type=Path, required=True,
                        metavar=('OLD', 'NEW'))
    args = parser.parse_args()
    inputs = {str(p.resolve()): digest(p) for pair in args.pair for p in pair}
    records = []
    for old_path, new_path in args.pair:
        old, new = read_table(old_path), read_table(new_path)
        nf, nr, nt = old['shape']
        new_nt = new['shape'][2]
        if new['shape'][:2] != (nf, nr) or new_nt <= nt:
            raise ValueError('expected only additional temperature rows')
        for key in ['abundance', 'ifmol', 'tmolim', 'flags', 'log_density', 'frequency']:
            if new[key] != old[key]:
                raise ValueError('changed opacity input: '+key)
        if new['log_temperature'][:nt] != old['log_temperature']:
            raise ValueError('previous temperature coordinates changed')
        if new['log_electron_density'][:nt*nr] != old['log_electron_density']:
            raise ValueError('previous electron densities changed')
        for i in range(nf*nr):
            if new['log_opacity'][i*new_nt:i*new_nt+nt] != old['log_opacity'][i*nt:(i+1)*nt]:
                raise ValueError('previous absorption coefficients changed')
        records.append({'old': str(old_path.resolve()), 'new': str(new_path.resolve()),
                        'retained_temperature_rows': nt, 'new_temperature_rows': new_nt-nt,
                        'retained_absorption_values': nf*nr*nt,
                        'new_absorption_values': nf*nr*(new_nt-nt),
                        'retained_cells_exact': True, 'all_source_values_finite': True})
    if any(digest(path) != expected for path, expected in inputs.items()):
        raise ValueError('opacity inputs changed during audit')
    report = {'scope': __doc__, 'passed': True, 'pairs': records,
              'input_files_sha256': inputs, 'script_sha256': digest(__file__)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'passed': True, 'pairs': len(records)}))


if __name__ == '__main__':
    main()
