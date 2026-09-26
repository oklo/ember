#!/usr/bin/env python3
"""Measure absorption/scattering sensitivity in existing TOPS spectral controls.

These collisionless comparisons isolate conventions. They do not adopt a
dielectric function, approximate all collective scattering corrections, or
change the selected stellar opacity. The absorption-only result excludes
scattering; it is not a physical replacement for a scattering calculation.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from audit_tops_spectral_means import source, digest


def means(data, temperature, cutoff, order):
    u = data[:, 0]/temperature
    last = int(np.searchsorted(u, 700))
    first = int(np.searchsorted(u, cutoff, side='right')-1)
    if not 0 <= first < last < len(u) or u[0] > .002:
        raise ValueError('incomplete frequency support')
    i = np.arange(first, last)
    lo, hi = u[i].copy(), u[i+1]
    lo[0] = cutoff
    nodes, weights = leggauss(order)
    x = (hi-lo)[:, None]*nodes/2+(hi+lo)[:, None]/2
    t = np.log(x/u[i, None])/np.log(u[i+1]/u[i])[:, None]
    logk = np.log(data[:, 1:])
    k = np.exp(logk[i, None, :]+t[:, :, None]*(logk[i+1]-logk[i])[:, None, :])
    total, absorption, scatter = (k[:, :, j] for j in range(3))
    n = np.sqrt(np.maximum(0., 1-(cutoff/x)**2))
    wr = x**4*np.exp(-x)/(-np.expm1(-x))**2
    inverse = {
        'n2_over_reported_total': n**2/total,
        'n3_over_reported_total': n**3/total,
        'n3_over_summed_components': n**3/(absorption+scatter),
        'absorption_divided_by_n_scattering_unchanged': n**3/(absorption+n*scatter),
        'absorption_divided_by_n_scattering_times_n': n**3/(absorption+n*n*scatter),
        'absorption_only_divided_by_n': n**3/absorption,
    }
    out = {key: (4*math.pi**4/15)/np.sum((hi-lo)[:, None]/2*weights*wr*v)
           for key, v in inverse.items()}
    # Independent algebraic check: account for n^2 in diffusion and 1/n in absorption.
    reconstructed = n**2/(absorption/n+scatter)
    if not np.allclose(reconstructed, inverse['absorption_divided_by_n_scattering_unchanged'], rtol=1e-14, atol=0):
        raise ValueError('inconsistent absorption convention')
    if not all(math.isfinite(v) and v > 0 for v in out.values()):
        raise ValueError('invalid mean')
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['work', 'classical_check', 'output']:
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    check = json.loads(args.classical_check.read_text())
    inputs = {str(p.resolve()): digest(p) for p in [Path(__file__), args.classical_check,
              Path(__file__).with_name('audit_tops_spectral_means.py')]}
    for path, checksum in check['input_sha256'].items():
        if digest(Path(path)) != checksum:
            raise ValueError('classical source-control input changed')
        inputs[path] = checksum
    spectra, _ = source(args.work/'cutoff-on')
    records = []
    for record in check['records']:
        key = next(k for k in spectra if k[1] == record['density'] and
                   abs(k[0]*1.602176634e-9/1.380649e-16/record['temperature_K']-1) < 1e-12)
        result = [means(spectra[key], key[0], record['classical_cutoff_u'], n) for n in [8, 16]]
        error = max(abs(result[0][k]/result[1][k]-1) for k in result[1])
        if error > 1e-7:
            raise ValueError('quadrature comparison failed')
        baseline = record['full_normalization_classical_refractive_rosseland']
        if abs(result[1]['n3_over_reported_total']/baseline-1) > 1e-12:
            raise ValueError('n-cubed control does not reproduce independent integration')
        kc = record['conductive_opacity']
        combined = {k: 1/(1/v+1/kc) for k, v in result[1].items()}
        records.append({'temperature_keV': key[0], 'density': key[1],
                        'classical_cutoff_u': record['classical_cutoff_u'],
                        'radiative_opacities': result[1], 'combined_opacities': combined,
                        'combined_relative_to_n3_total': {k: v/combined['n3_over_reported_total']-1 for k, v in combined.items()},
                        'quadrature_relative_difference': error})
    report = {'scope': __doc__, 'accepted_for_stellar_opacity': False, 'input_sha256': inputs,
              'references': ['https://doi.org/10.1103/PhysRevE.91.053102',
                             'https://articles.adsabs.harvard.edu/pdf/1996ApJ...466L.115I'],
              'assumptions': 'Classical collisionless n from the existing source-electron controls. Scattering variants vary only its explicit factor; the density-dependent electron structure factor is not recalculated.',
              'records': records}
    if any(digest(Path(p)) != h for p, h in inputs.items()):
        raise ValueError('diagnostic input changed')
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps([r['combined_relative_to_n3_total'] for r in records]))


if __name__ == '__main__':
    main()
