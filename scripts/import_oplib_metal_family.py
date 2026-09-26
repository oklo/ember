#!/usr/bin/env python3
"""Import original OPLIB composition planes within a complete source rectangle.

Missing native values are rejected. No temperature, density or composition
extrapolation is performed. These are comparison/response sources, not a
replacement selected for the retained interior opacities.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import tarfile

import numpy as np

from audit_mesa_oplib_native import read_native
from write_scientific_result import write_result


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['archive', 'work', 'report']:
        p.add_argument(name, type=Path)
    p.add_argument('--metals', nargs='+', type=float, default=[.004, .006, .01, .02, .03])
    a = p.parse_args()
    expected = '042eb369c8f15c808758bee7938ed138f686a7d31bb48475bbb719be24cec997'
    if a.work.exists() or a.report.exists() or digest(a.archive) != expected:
        raise ValueError('existing output or different source archive')
    metals = a.metals
    if (len(metals) < 2 or sorted(set(metals)) != metals
            or not all(np.isfinite(z) and 0 < z <= .05 for z in metals)):
        raise ValueError('ordered positive source metallicities through Z=.05 required')
    hydrogen = [i/20 for i in range(20)]
    inputs = {str(q.resolve()): digest(q) for q in
              [a.archive, Path(__file__), Path(__file__).with_name('audit_mesa_oplib_native.py')]}
    a.work.mkdir(parents=True)
    source = a.work/'source'
    source.mkdir()
    records, tables = [], {}
    pattern = re.compile(r'kap_input_data/oplib/OPLIB_GS98_LOG_1194_25E/oplib_gs98_z([0-9.]+)_x([0-9.]+)\.data')
    with tarfile.open(a.archive) as archive:
        for member in archive:
            match = pattern.fullmatch(member.name)
            if not match:
                continue
            z, x = map(float, match.groups())
            if z not in metals or x not in hydrogen:
                continue
            if not member.isfile() or member.size > 1000000 or (z, x) in tables:
                raise ValueError('invalid source member')
            raw = archive.extractfile(member).read()
            path = source/Path(member.name).name
            path.write_bytes(raw)
            entry = dict(path=str(path.resolve()), sha256=digest(path), archive_member=member.name)
            table = read_native(entry)
            if table['X'] != x or table['Z'] != z:
                raise ValueError('source composition label differs')
            records.append(entry)
            tables[z, x] = table
    if set(tables) != {(z, x) for z in metals for x in hydrogen}:
        raise ValueError('missing native composition plane')
    example = tables[metals[0], hydrogen[0]]
    count = example['logT'].index(7.392)+1
    ts, rs = example['logT'][:count], example['logR']
    outputs = []
    for z in metals:
        label = f'{z:.3f}'
        if float(label) != z:
            label = format(z, '.17g')
        path = a.work/f'oplib_z{label}.dat'
        lines = [f'{len(hydrogen)} {len(ts)} {len(rs)} Original OPLIB GS98; composition-response source only',
                 ' '.join(format(v, '.17g') for v in rs),
                 ' '.join(format(v, '.17g') for v in ts)]
        for x in hydrogen:
            table = tables[z, x]
            if (table['logT'][:count] != ts or table['logR'] != rs
                    or np.any(table['mask'][:count])):
                raise ValueError('inconsistent axes or missing source value')
            lines.append(f'{x:.17g} {z:.17g}')
            lines.extend(' '.join(format(v, '.17g') for v in row)
                         for row in table['logk'][:count])
        path.write_text('\n'.join(lines)+'\n')
        # Every stored number must round-trip to the unaltered source value.
        written = path.read_text().splitlines()
        for i, x in enumerate(hydrogen):
            begin = 3+i*(count+1)
            if list(map(float, written[begin].split())) != [x, z]:
                raise ValueError('output composition differs')
            values = np.loadtxt(written[begin+1:begin+1+count])
            if not np.array_equal(values, tables[z, x]['logk'][:count]):
                raise ValueError('output source values changed')
        outputs.append(dict(Z=z, path=str(path), sha256=digest(path)))
    manifest = a.work/'oplib_metal_family.dat'
    manifest.write_text(f'EMBER_OPACITY_MIXTURE 1 {len(metals)} logR OPLIB_GS98_composition_response\n'+
                        ''.join(f'{r["Z"]:.17g} "{Path(r["path"]).name}"\n' for r in outputs))
    for path, expected in inputs.items():
        if digest(path) != expected:
            raise ValueError('source changed during import')
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  outcome='imported_unchanged_native_source_values', selected_for_evolution=False,
                  input_sha256=inputs, source_planes=records,
                  hydrogen=hydrogen, metals=metals, logT_range=[ts[0], ts[-1]],
                  logR_range=[rs[0], rs[-1]], output_files=outputs,
                  manifest=str(manifest), manifest_sha256=digest(manifest),
                  source_values=len(hydrogen)*len(metals)*len(ts)*len(rs),
                  missing_values_used=0,
                  limitations=['Original 25-element OPLIB mixture and plasma conventions differ from the retained TOPS material.',
                               'This source cannot provide densities beyond logR=1.5; dense warm-envelope coverage is a separate problem.',
                               'An explicit response approximation and finite stellar checks are required before selecting changes.'])
    write_result(a.report, report)
    print(report['source_values'], 'original source values imported')


if __name__ == '__main__':
    main()
