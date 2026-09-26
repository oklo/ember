#!/usr/bin/env python3
"""Export existing, checked collision interpolants for the C++ evaluator.

No collision integrals are recalculated. The text file contains spline
coefficients and quadrature rules, not frozen stellar transport coefficients.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PPoly
from scipy.special import roots_jacobi, roots_legendre
from electron_pair_table import ElectronPairTable
from ion_collision_table import ExtendedIonTable


def export(output):
    root = Path(__file__).resolve().parents[1]
    pair_path = root/'docs/results/electron_pair_table_sources_v2.json'
    acceptance_path = root/'docs/results/electron_pair_interpolation_v2.json'
    acceptance = json.loads(acceptance_path.read_text())
    assert acceptance['outcome'] == 'passed_independent_controls'
    pair = ElectronPairTable(pair_path)
    assert pair.report_sha256 == acceptance['input_sha256'][str(pair_path)]
    ext_report = root/'docs/results/ion_collision_table_extension_v1.json'
    ext = json.loads(ext_report.read_text())
    table_path = Path(ext['table'])
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    assert digest(table_path) == ext['table_sha256']
    ions = ExtendedIonTable(table_path)
    extension = PPoly.from_bernstein_basis(ions.extension)
    inputs = {str(p): digest(p) for p in [pair_path, acceptance_path, ext_report,
              table_path, Path(__file__), root/'scripts/electron_pair_table.py',
              root/'scripts/ion_collision_table.py']}
    base = json.loads(table_path.read_text())
    inputs[base['base_table']] = base['base_table_sha256']
    sources = json.loads(pair_path.read_text())['completed_sources']
    inputs.update({s['path']: s['sha256'] for s in sources})
    for p, h in inputs.items():
        assert digest(p) == h
    with output.open('x') as f:
        f.write('EMBER_COLLISION_TRANSPORT_V1\n')
        def line(a):
            f.write(' '.join(format(float(v), '.17g') for v in np.asarray(a).ravel())+'\n')
        tx, ty, _ = pair.splines[0].tck
        line([len(tx), len(ty)])
        line(tx); line(ty)
        for spline in pair.splines:
            assert np.array_equal(tx, spline.tck[0]) and np.array_equal(ty, spline.tck[1])
            line(spline.get_coeffs())
        for polynomial in [ions.base, extension]:
            line([len(polynomial.x), polynomial.c.shape[0]])
            line(polynomial.x)
            line(polynomial.c)
        # All integrals have smooth endpoints on these variables. The Jacobi
        # rule preserves the checked x^(3/2) polynomial normalization.
        for count, kind in [(384, 'jacobi'), (512, 'legendre')]:
            x, w = roots_jacobi(count, 0, 1.5) if kind == 'jacobi' else roots_legendre(count)
            line([count]); line(x); line(w)
    report = dict(outcome='exported_checked_interpolants', input_sha256=inputs,
                  output=str(output), output_sha256=digest(output),
                  output_bytes=output.stat().st_size, new_collision_integrals=0,
                  pair_eta_range=[float(pair.eta[0]), float(pair.eta[-1])],
                  pair_b_range=[float(pair.b[0]), float(pair.b[-1])],
                  ion_log10_strength_range=[ions.minimum, ions.maximum])
    output.with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='input_sha256'}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    export(parser.parse_args().output)
