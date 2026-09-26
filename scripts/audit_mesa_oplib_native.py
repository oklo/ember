#!/usr/bin/env python3
"""Check original MESA OPLIB values without filling their missing entries."""
import hashlib
import json
import math
from pathlib import Path
import statistics

import numpy as np

from compare_mesa_oplib_means_v2 import interpolate, stencil

ROOT = Path(__file__).resolve().parents[1]
WORK = Path('/tmp/ember-mesa-native-input-fetch-v1')
COMPARISON = ROOT / 'docs/results/mesa_oplib_mean_comparison_v2.json'
OUTPUT = ROOT / 'docs/results/mesa_oplib_native_comparison_v1.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_native(entry):
    assert digest(entry['path']) == entry['sha256']
    lines = Path(entry['path']).read_text().splitlines()
    assert '25 Elements' in lines[0] and 'Grevesse & Sauval (1998)' in lines[0]
    metadata = list(map(float, lines[2].split()))
    assert len(metadata) == 10 and metadata[:2] == [1, 37]
    assert metadata[4:] == [39, -8, 1.5, 74, 3.764, 9.065]
    rr = list(map(float, lines[5].split()))
    a = np.loadtxt(lines[6:])
    assert a.shape == (74, 40) and np.isfinite(a).all()
    assert rr == [-8 + .25*i for i in range(39)] and np.all(np.diff(a[:,0]) > 0)
    mask = a[:,1:] >= 9.99
    assert np.all(np.abs(a[:,1:][~mask]) < 8)
    missing = [dict(logT=float(a[i,0]),logR=rr[j],source_value=float(a[i,j+1]))
               for i,j in np.argwhere(mask)]
    return dict(X=metadata[2],Z=metadata[3],logT=a[:,0].tolist(),logR=rr,
                logk=a[:,1:].tolist(),mask=mask,missing=missing,source=entry)


def evaluate(table,t,r,order):
    ii,_ = stencil(table['logT'],t,order)
    jj,_ = stencil(table['logR'],r,order)
    if any(table['mask'][i,j] for i in ii for j in jj):
        raise ValueError('interpolation touches an explicitly missing native value')
    return interpolate(table,t,r,order)


def summary(rows,key):
    return dict(median=statistics.median(row[key] for row in rows),
                minimum=min(row[key] for row in rows),maximum=max(row[key] for row in rows),
                worst=max(rows,key=lambda row:abs(row[key])))


def main():
    assert not OUTPUT.exists()
    manifest = json.loads((WORK/'manifest.json').read_text())
    assert len(manifest['selected']) == 19
    assert digest(COMPARISON) == '90e28efe30a22258d6c36dd96f8d6ef6c6e52c885f8cee4f04961cf02e682e20'
    tables = [read_native(entry) for entry in manifest['selected']]
    checks = []
    for table in tables:
        for cell in table['missing']:
            try:
                evaluate(table,cell['logT'],cell['logR'],2)
            except ValueError:
                pass
            else:
                raise AssertionError('missing native value was accepted')
        checks.append(dict(X=table['X'],Z=table['Z'],missing_values_rejected=len(table['missing'])))
    compared = json.loads(COMPARISON.read_text())
    rows = []
    for comparison in compared['comparisons']:
        table = next(t for t in tables if (t['X'],t['Z']) == (comparison['X'],comparison['Z']))
        for cell in comparison['rows']:
            t,r = cell['logT'],cell['logR']
            linear = 10**evaluate(table,t,r,2)
            cubic = 10**evaluate(table,t,r,4)
            # These source temperatures coincide with native OPLIB isotherms
            # to the printed 0.001-dex precision. Evaluate on that printed
            # isotherm as a sensitivity check, not as a corrected temperature.
            nearest = min(table['logT'],key=lambda v:abs(v-t))
            assert abs(nearest-t) <= .0005
            native_row = 10**evaluate(table,nearest,r,4)
            rows.append(dict(X=table['X'],Z=table['Z'],T_keV=cell['T_keV'],rho=cell['rho'],
                             logT=t,logR=r,nearest_native_logT=nearest,
                             tops_uncut=cell['tops_uncut'],processed_mesa=cell['mesa_four_nodes'],
                             native_linear=linear,native_cubic=cubic,native_row_cubic=native_row,
                             native_relative_to_tops=cubic/cell['tops_uncut']-1,
                             processed_relative_to_native=cell['mesa_four_nodes']/cubic-1,
                             native_interpolation_change=cubic/linear-1,
                             printed_temperature_change=native_row/cubic-1))
    profile_path = ROOT/'docs/reports/2026-09-11/evolution_latest_profile.csv'
    profile = np.genfromtxt(profile_path,delimiter=',',names=True)
    # Report names rather than guessing profile field meanings. Domain
    # assessment is separate from the fixed-source comparisons above.
    table_info = [{k:v for k,v in table.items() if k not in ['logk','mask']} for table in tables]
    result = dict(outcome='completed_native_table_comparison',accepted_for_stellar_opacity=False,
                  input_sha256={str(p):digest(p) for p in [Path(__file__),ROOT/'scripts/compare_mesa_oplib_means_v2.py',
                                 COMPARISON,WORK/'manifest.json',profile_path]},
                  tables=table_info,missing_value_rejection_checks=checks,rows=rows,
                  summaries={key:summary(rows,key) for key in ['native_relative_to_tops','processed_relative_to_native',
                                                              'native_interpolation_change','printed_temperature_change']},
                  profile_columns=list(profile.dtype.names),
                  limitations=['The original archive flags missing values at logT >= 8.463; no such values enter these cooler comparisons.',
                               'The original tables are mixture means on a logR grid, not the underlying element spectra or their EOS tables.',
                               'The 25-versus-21-element mixture difference remains. Interpolation sensitivity is not a physical error bound.',
                               'The native temperature is printed to 0.001 dex; its underlying exact value is not assumed recovered.',
                               'No missing value or exterior domain is extrapolated; the MESA generation convention is not inferred from numerical agreement.'])
    with OUTPUT.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(dict(report=str(OUTPUT),sha256=digest(OUTPUT),row_count=len(rows),
                          missing_values=sum(len(t['missing']) for t in tables),summaries=result['summaries']),indent=2))


if __name__ == '__main__':
    main()
