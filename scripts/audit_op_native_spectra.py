#!/usr/bin/env python3
"""Check decoded OP spectra against their independently stored native means."""
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from op_native_spectra import read_mesh, read_spectra, native_rosseland

ROOT = Path('/tmp/ember-op-selected-spectra-v1')
OUTPUT = Path('docs/results/op_native_spectral_reader_v1.json')


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    if OUTPUT.exists():raise FileExistsError(OUTPUT)
    manifest = json.loads((ROOT/'extraction_manifest.json').read_text())
    inputs = {str(Path(p).resolve()): digest(p) for p in
              [__file__, 'scripts/op_native_spectra.py', ROOT/'extraction_manifest.json']}
    meshes = {}; records = []; failures = []
    for entry in manifest['files']:
        p = Path(entry['path'])
        if digest(p) != entry['sha256']:raise ValueError('extracted source changed')
        inputs[str(p)] = entry['sha256']
        if p.name.endswith('.mesh'):
            meshes[int(p.name[1:3])] = read_mesh(p)
    base_dv, base_u = meshes[1]
    for z, (dv, u) in meshes.items():
        if dv != base_dv or not np.array_equal(u, base_u):
            raise ValueError('element frequency grids differ')
    for entry in sorted(manifest['files'], key=lambda r: r['path']):
        p = Path(entry['path'])
        if not re.fullmatch(r'm\d\d\.\d\d\d', p.name):continue
        r = read_spectra(p); z = r['atomic_number']; dv, u = meshes[z]
        checks = []
        for row in r['states']:
            actual = float(native_rosseland(row['cross_section_atomic'], dv))
            error = actual / row['rosseland_atomic'] - 1
            charge = float(np.dot(z - 1 - row['ion_indices'], row['ion_fractions']))
            normalization = float(row['ion_fractions'].sum())
            check = dict(electron_index=row['electron_index'], reported_rosseland_atomic=row['rosseland_atomic'],
                         decoded_rosseland_atomic=actual, relative_mean_difference=error,
                         ion_fraction_sum=normalization, reported_electrons_per_atom=row['electrons_per_atom'],
                         ion_sum_electrons_per_atom=charge, charge_absolute_difference=charge-row['electrons_per_atom'],
                         packed_points=row['packed_points'])
            # Allow the documented packing error plus decimal source rounding.
            check['mean_passed'] = abs(error) <= r['packing_tolerance'] + .00011
            check['ion_passed'] = abs(normalization-1) <= .0003 and abs(charge-row['electrons_per_atom']) <= .0003 * z
            if not check['mean_passed'] or not check['ion_passed']:
                failures.append(dict(file=p.name, **check))
            checks.append(check)
        records.append(dict(file=str(p), atomic_number=z, temperature_index=r['temperature_index'],
                            packing_tolerance=r['packing_tolerance'], checks=checks))
    for p,h in inputs.items():
        if digest(p) != h:raise ValueError('source changed during read')
    r = dict(outcome='passed' if not failures else 'failed', input_sha256=inputs,
             files=len(records), states=sum(len(r['checks']) for r in records),
             shared_frequency_mesh=dict(points=len(base_u), u_min=float(base_u[0]), u_max=float(base_u[-1]), dv=base_dv),
             maximum_mean_difference=max(abs(c['relative_mean_difference']) for r in records for c in r['checks']),
             records=records, failures=failures, accepted_for_stellar_opacity=False,
             limitations=['Native finite-grid means only; no frequency-tail, thermodynamic-interpolation or mixture-screening acceptance.',
                          'Source ion fractions and means have limited printed precision; tolerances explicitly include it.',
                          'This does not establish agreement with TOPS or completeness for the tracked element mixture.'])
    OUTPUT.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:r[k] for k in ('outcome','files','states','maximum_mean_difference')}),flush=True)
    if failures:raise SystemExit(1)


if __name__ == '__main__':main()
