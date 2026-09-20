#!/usr/bin/env python3
"""Compare TOPS grey means with direct integration of the returned spectra.

This is an averaging diagnostic, not a replacement opacity table. The cutoff
is inferred from the reported Rosseland mean; the Planck mean is then an
independent check of that inference. Rosseland integrals use total opacity;
Planck integrals use absorption. Both use log-linear interpolation between
positive spectral samples.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq

from import_tops_composition import read


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source(root):
    receipt = json.loads((root/'receipt.json').read_text())
    if digest(root/'request.json') != receipt['request_sha256']:
        raise ValueError('changed source request')
    text = (root/'source.txt').read_text()
    tt, rr, cells, excluded = read(root/'source.txt', receipt,
                                  dimensions=tuple(receipt['dimensions']))
    if excluded:
        raise ValueError('substituted density is not a spectral control')
    means = {}
    for part in text.split('Density     Ross opa    Planck opa  No. Free    Av Sq Free  T=  ')[1:]:
        lines = part.splitlines(); t = float(lines[0])
        for line in lines[1:]:
            v = line.split()
            if len(v) != 5 or not all(re.fullmatch(r'[0-9.E+-]+', x) for x in v):
                break
            rho, ross, planck, free, square = map(float, v)
            means[t, rho] = dict(rosseland=ross, planck=planck, free_electrons_per_ion=free)
    spectra = {}
    for part in text.split('T (keV), Density (gm/cc) =')[1:]:
        lines = part.splitlines(); key = tuple(map(float, lines[0].split()))
        if lines[1].strip() != 'Photon energy(keV), Total, Absorp, Scatt (cm**2/gm)':
            raise ValueError('unknown spectrum columns')
        rows = []
        for line in lines[2:]:
            v = line.split()
            if len(v) != 4 or not all(re.fullmatch(r'[0-9.E+-]+', x) for x in v):
                break
            rows.append(list(map(float, v)))
        data = np.array(rows)
        if (len(rows) < 1000 or not np.isfinite(data).all() or np.any(data <= 0)
                or np.any(np.diff(data[:, 0]) <= 0)):
            raise ValueError('incomplete, unordered or nonpositive spectrum')
        if key in spectra:
            raise ValueError('duplicate spectrum')
        spectra[key] = data
    if set(spectra) != set(cells) or set(means) != set(cells):
        raise ValueError('spectra and reported means differ in coordinates')
    return spectra, means


class Integrals:
    """Integrate interval by interval; moving a cutoff never jumps a whole bin."""
    def __init__(self, data, temperature, order):
        self.u = data[:, 0]/temperature
        self.log_kappa = np.log(data[:, 1])
        self.log_absorption = np.log(data[:, 2])
        self.nodes, self.weights = leggauss(order)
        # Weight above u=700 is negligible even for the dense controls here.
        self.last = int(np.searchsorted(self.u, 700.))
        if self.last >= len(self.u) or self.u[0] > .002:
            raise ValueError('spectrum has insufficient frequency coverage')
        self.intervals = self.integrate(self.u[:self.last], self.u[1:self.last+1],
                                        np.arange(self.last))
        self.tail = np.vstack([np.cumsum(self.intervals[::-1], axis=0)[::-1], np.zeros(4)])

    def integrate(self, left, right, index, refractive_cutoff=None):
        left, right, index = np.atleast_1d(left), np.atleast_1d(right), np.atleast_1d(index)
        u = .5*((right-left)[:, None]*self.nodes+(right+left)[:, None])
        fraction = np.log(u/self.u[index, None])/np.log(self.u[index+1]/self.u[index])[:, None]
        k = np.exp(self.log_kappa[index, None]+fraction*(self.log_kappa[index+1]-self.log_kappa[index])[:, None])
        absorption = np.exp(self.log_absorption[index, None]+fraction*(self.log_absorption[index+1]-self.log_absorption[index])[:, None])
        z = np.exp(-u); denominator = -np.expm1(-u)
        wr = u**4*z/denominator**2
        wp = u**3*z/denominator
        inverse = wr/k
        if refractive_cutoff is not None:
            inverse *= np.maximum(0., 1-(refractive_cutoff/u)**2)**1.5
        f = np.stack([wr, inverse, wp, wp*absorption], axis=-1)
        return .5*(right-left)[:, None]*np.einsum('ijk,j->ik', f, self.weights)

    def above(self, cutoff):
        i = int(np.searchsorted(self.u, cutoff, side='right')-1)
        if not 0 <= i < self.last:
            raise ValueError('cutoff outside integrated source support')
        return self.integrate(cutoff, self.u[i+1], i)[0]+self.tail[i+1]

    def refractive(self, cutoff):
        i = int(np.searchsorted(self.u, cutoff, side='right')-1)
        indices = np.arange(i, self.last)
        left = self.u[indices].copy(); left[0] = cutoff
        return self.integrate(left, self.u[indices+1], indices, cutoff).sum(axis=0)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('work', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--cutoff-bracket', type=float, nargs=2, required=True)
    a = p.parse_args()
    inputs = {str(path.resolve()): digest(path) for root in [a.work/'cutoff-on', a.work/'cutoff-off']
              for path in [root/'receipt.json', root/'request.json', root/'source.txt']}
    inputs[str(Path(__file__).resolve())] = digest(Path(__file__))
    on, off = [source(a.work/('cutoff-'+label)) for label in ['on', 'off']]
    if set(on[0]) != set(off[0]):
        raise ValueError('on/off coordinates differ')
    records = []
    for key, data in on[0].items():
        if not np.array_equal(data, off[0][key]):
            raise ValueError('plasma-cutoff option changed the spectrum')
        calculations = []
        for order in [8, 16]:
            integral = Integrals(data, key[0], order)
            whole = integral.tail[0]
            def residual(cutoff):
                v = integral.above(cutoff)
                return v[0]/v[1]-on[1][key]['rosseland']
            cutoff = brentq(residual, *a.cutoff_bracket, xtol=1e-11)
            above = integral.above(cutoff)
            calculations.append(dict(
                quadrature_order=order, inferred_cutoff_u=cutoff,
                uncut_rosseland=(4*math.pi**4/15)/whole[1],
                uncut_planck=whole[3]/(math.pi**4/15),
                conditional_rosseland=above[0]/above[1],
                conditional_planck=above[3]/above[2],
                rosseland_weight_fraction_above_cutoff=above[0]/(4*math.pi**4/15),
                planck_weight_fraction_above_cutoff=above[2]/(math.pi**4/15),
                full_normalization_step_cutoff_rosseland=(4*math.pi**4/15)/above[1],
                diagnostic_collisionless_refractive_rosseland=(4*math.pi**4/15)/integral.refractive(cutoff)[1]))
        final = calculations[-1]
        errors = dict(
            uncut_rosseland=final['uncut_rosseland']/off[1][key]['rosseland']-1,
            uncut_planck=final['uncut_planck']/off[1][key]['planck']-1,
            independently_predicted_cutoff_planck=final['conditional_planck']/on[1][key]['planck']-1)
        records.append(dict(temperature_keV=key[0], density=key[1], spectral_points=len(data),
                            spectra_identical=True, reported_on=on[1][key], reported_off=off[1][key],
                            calculations=calculations, relative_differences=errors,
                            quadrature_relative_difference=max(abs(final[k]/calculations[0][k]-1)
                                for k in final if k != 'quadrature_order')))
    report = dict(scope=__doc__, input_sha256=inputs, states=records,
                  accepted_for_stellar_opacity=False,
                  interpretation='A full-normalization step cutoff is a diagnostic of the averaging convention, not a complete refractive-index transport prescription.',
                  refractive_diagnostic='Also integrates n(u)^3=[1-(u_cut/u)^2]^(3/2) above the inferred cutoff. This collisionless free-electron expression is a separate sensitivity diagnostic, not a selected dielectric model.',
                  frequency_limits='Integrals use the provided spectrum from its first frequency through the first node above u=700. No missing spectral opacity is invented.',
                  reference='Colgan et al. 2016, equation 1 and section 2.8; https://www.osti.gov/servlets/purl/1335616')
    if any(digest(Path(path)) != value for path, value in inputs.items()):
        raise ValueError('an input changed during the audit')
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps([r['relative_differences'] for r in records]))


if __name__ == '__main__':
    main()
