#!/usr/bin/env python3
"""Separate density-interpolation sensitivity from MESA/TOPS differences.

Reverse the comparison onto MESA's exact logR nodes. TOPS means are then
interpolated only across contiguous, verified density samples; no gap larger
than 0.4001 dex or substituted source state is bridged. This is a source
comparison, not a physical error bound or a selected stellar opacity.
"""
import hashlib
import json
import math
from pathlib import Path
import statistics

import numpy as np

from audit_mesa_oplib_native import read_native, evaluate
from compare_mesa_oplib_means_v2 import stencil, KELVIN_PER_KEV
from import_tops_composition import read as read_tops

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT/'docs/results/mesa_oplib_density_comparison_v1.json'
BASE = ROOT/'docs/results/mesa_oplib_mean_comparison_v2.json'
MANIFEST = Path('/tmp/ember-mesa-native-input-fetch-v1/manifest.json')


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def sample(grid, values, q, order):
    ii, ww = stencil(grid, q, order)
    if max(np.diff(np.asarray(grid)[ii])) > .4001:
        raise ValueError('disconnected saved density support')
    return math.fsum(w*values[i] for i, w in zip(ii, ww)), ii


def stats(rows, key):
    return dict(count=len(rows), median=statistics.median(r[key] for r in rows),
                minimum=min(r[key] for r in rows), maximum=max(r[key] for r in rows),
                worst=max(rows, key=lambda r: abs(r[key])))


def main():
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    inputs = {str(p): digest(p) for p in [Path(__file__), BASE, MANIFEST,
        ROOT/'scripts/audit_mesa_oplib_native.py',
        ROOT/'scripts/compare_mesa_oplib_means_v2.py',
        ROOT/'scripts/import_tops_composition.py']}
    base = json.loads(BASE.read_text())
    tables = [read_native(v) for v in json.loads(MANIFEST.read_text())['selected']]
    rows, skipped, original = [], [], []
    for c in base['comparisons']:
        table = next(t for t in tables if (t['X'],t['Z']) == (c['X'],c['Z']))
        entry = table['source']; inputs[entry['path']] = entry['sha256']
        states = {}
        for source in c['sources']:
            work = Path(source['work']); p = work/'receipt.json'
            assert digest(p) == source['receipt_sha256']
            receipt = json.loads(p.read_text()); inputs[str(p)] = digest(p)
            raw = work/receipt['file']; inputs[str(raw)] = digest(raw)
            _, _, cells, excluded = read_tops(raw, receipt, dimensions=tuple(receipt['dimensions']))
            for key, value in cells.items():
                if key not in excluded:
                    if key in states and states[key] != value:
                        raise ValueError('conflicting duplicate source state')
                    states[key] = value
        for tkev in sorted(set(t for t,r in states)):
            t = math.log10(tkev*KELVIN_PER_KEV)
            pairs = sorted((math.log10(r), math.log10(k)) for (tt,r),k in states.items() if tt == tkev)
            densities, opacities = zip(*pairs)
            for r in table['logR']:
                logrho = r+3*t-18
                try:
                    lin, _ = sample(densities, opacities, logrho, 2)
                    cub, ii = sample(densities, opacities, logrho, 4)
                except ValueError as error:
                    skipped.append(dict(X=c['X'],Z=c['Z'],T_keV=tkev,logR=r,reason=str(error)))
                    continue
                # Exactly on a density node: only the temperature direction
                # contributes interpolation to the MESA mean.
                mesa = 10**evaluate(table,t,r,4)
                row = dict(X=c['X'],Z=c['Z'],T_keV=tkev,logT=t,logR=r,
                           rho=10**logrho,mesa=mesa,tops_linear=10**lin,tops_cubic=10**cub,
                           tops_source_densities=[10**densities[i] for i in ii],
                           mesa_relative_to_tops=mesa/10**cub-1,
                           tops_interpolation_change=10**(cub-lin)-1)
                rows.append(row)
        for cell in c['rows']:
            t,r = cell['logT'],cell['logR']
            ti,tw = stencil(table['logT'],t,4)
            values = [math.fsum(w*table['logk'][i][j] for i,w in zip(ti,tw))
                      for j in range(len(table['logR']))]
            ii4,ww4 = stencil(table['logR'],r,4)
            # An independently enlarged density stencil measures sensitivity;
            # it is not automatically more accurate near an ionization feature.
            lower = int(np.searchsorted(table['logR'],r,side='right'))-1
            first = max(0,min(len(table['logR'])-6,lower-2))
            ii6 = list(range(first,first+6))
            ww6 = [math.prod((r-table['logR'][j])/(table['logR'][i]-table['logR'][j])
                            for j in ii6 if j != i) for i in ii6]
            k4 = 10**math.fsum(w*values[i] for i,w in zip(ii4,ww4))
            k6 = 10**math.fsum(w*values[i] for i,w in zip(ii6,ww6))
            original.append(dict(X=c['X'],Z=c['Z'],T_keV=cell['T_keV'],rho=cell['rho'],
                                 mesa_four=k4,mesa_six=k6,tops=cell['tops_uncut'],
                                 density_stencil_change=k6/k4-1,
                                 mesa_six_relative_to_tops=k6/cell['tops_uncut']-1))
    # This extra interpolator is checked by a known cubic and an explicit gap.
    grid=[0.,.1,.3,.6,.8]; values=[1+.2*x-.3*x*x+.7*x**3 for x in grid]
    value,_=sample(grid,values,.43,4)
    cubic_error=abs(value-(1+.2*.43-.3*.43**2+.7*.43**3))
    assert cubic_error < 1e-14
    try:
        sample([0.,1.,2.,3.],[1.,2.,3.,4.],1.5,4)
    except ValueError:
        pass
    else:
        raise AssertionError('disconnected density interval accepted')
    assert rows and original
    result=dict(outcome='completed_density_source_comparison',accepted_for_stellar_opacity=False,
                input_sha256=inputs,rows=rows,skipped_nodes=skipped,original_source_rows=original,
                controls=dict(cubic_absolute_error=cubic_error,disconnected_interval_rejected=True),
                summaries={k:stats(v,k) for v,keys in [(rows,['mesa_relative_to_tops','tops_interpolation_change']),
                            (original,['density_stencil_change','mesa_six_relative_to_tops'])] for k in keys},
                limitations=['The source mixtures still have 25 versus 21 elements.',
                    'Exact MESA density nodes remove local density interpolation only on the MESA side.',
                    'The TOPS density interpolation and enlarged MESA stencil measure sensitivity, not absolute accuracy.',
                    'Temperature rounding and upstream mixture calculations remain separate uncertainties.',
                    'No opacity is extrapolated or selected for stellar evolution.'])
    for p,h in inputs.items():
        assert digest(p) == h
    with OUTPUT.open('x') as f:
        json.dump(result,f,indent=2,allow_nan=False); f.write('\n')
    print(json.dumps(dict(report=str(OUTPUT),sha256=digest(OUTPUT),
                         rows=len(rows),original_rows=len(original),summaries=result['summaries'])))


if __name__ == '__main__':
    main()
