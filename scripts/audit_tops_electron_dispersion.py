#!/usr/bin/env python3
"""Measure finite-temperature electron dispersion in existing TOPS controls.

Use Braaten & Segel (1993), equations 1, 4, 11 and 18, for collisionless free
electrons. Keep the atomic absorption and the scattering prescription fixed.
These comparisons do not include collision damping or bound-electron
dispersion and do not accept a production opacity table.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expit

from audit_tops_spectral_means import source, digest

ME = 9.1093837015e-28
C = 2.99792458e10
HBAR = 1.054571817e-27
KB = 1.380649e-16
KEV = 1.602176634e-9
RW = 4*math.pi**4/15


def electron_moments(temperature, ne):
    theta = KB*temperature/(ME*C*C)
    if not 0 < theta < .02 or not ne > 0:
        raise ValueError('electron control outside the pair-negligible temperature domain')
    target = ne*math.pi**2*(HBAR/(ME*C))**3
    qf = (3*target)**(1/3)
    ef = qf*qf/(math.sqrt(1+qf*qf)+1)

    def integral(eta, kind, precision=1e-10):
        top_energy = theta*max(eta+60, 60)
        top = math.sqrt(top_energy*(top_energy+2))
        points = []
        for delta in [-20, -5, 0, 5, 20]:
            energy = theta*(eta+delta)
            if 0 < energy < top_energy:
                points.append(math.sqrt(energy*(energy+2)))

        def f(q):
            gamma = math.sqrt(1+q*q)
            kinetic = q*q/(gamma+1)
            occupation = expit(eta-kinetic/theta)
            if kind == 'density':
                return q*q*occupation
            v2 = q*q/(gamma*gamma)
            factor = 1-v2/3 if kind == 'plasma' else 5*v2/3-v2*v2
            return q*q/gamma*factor*occupation

        return quad(f, 0, top, points=points, epsabs=max(target*1e-13, 1e-30),
                    epsrel=precision, limit=200)[0]

    eta = brentq(lambda x: integral(x, 'density')/target-1,
                 -50, ef/theta+50, xtol=1e-10, rtol=1e-14)
    density = integral(eta, 'density', 1e-12)
    plasma = integral(eta, 'plasma', 1e-12)
    omega1 = integral(eta, 'omega1', 1e-12)
    inversion_error = abs(density/target-1)
    ratio = plasma/density
    vstar2 = omega1/plasma
    if inversion_error > 1e-9 or not 0 < ratio <= 1 or not 0 <= vstar2 < .5:
        raise ValueError('electron integral or series domain check failed')
    return {'theta': theta, 'kinetic_chemical_potential_over_kT': eta,
            'fermi_momentum_over_mc': qf, 'plasma_frequency_squared_ratio': ratio,
            'vstar_squared': vstar2, 'density_inversion_relative_error': inversion_error,
            'zero_temperature_frequency_squared_ratio': 1/math.sqrt(1+qf*qf)}


def transverse_factor(z):
    # Equation 18 after setting n=kc/omega. This convergent series avoids
    # cancellation in its logarithmic form at n=0 and at small electron speed.
    coefficients = np.array([3/((2*j+1)*(2*j+3)) for j in range(32)])
    factor = np.polynomial.polynomial.polyval(z, coefficients)
    derivative = np.polynomial.polynomial.polyval(z, np.arange(1, 32)*coefficients[1:])
    return factor, derivative


def index_squared(x, cutoff, vstar2):
    r2 = (cutoff/x)**2
    y = (1-r2)/(1+r2*vstar2/5)
    for _ in range(6):
        f, df = transverse_factor(vstar2*y)
        y -= (y+r2*f-1)/(1+r2*vstar2*df)
    f, _ = transverse_factor(vstar2*y)
    if np.max(np.abs(y+r2*f-1)) > 1e-12 or np.any(y < 0) or np.any(y > 1):
        raise ValueError('transverse dispersion root failed')
    return y


def mean(data, temperature_keV, cutoff, vstar2, order):
    u = data[:, 0]/temperature_keV
    first, last = np.searchsorted(u, [cutoff, 700], side='right')-1
    if not 0 <= first < last < len(u)-1 or u[0] > .002:
        raise ValueError('incomplete spectral coverage')
    i = np.arange(first, last+1)
    lo, hi = u[i].copy(), u[i+1].copy()
    lo[0], hi[-1] = cutoff, 700.
    nodes, weights = leggauss(order)
    x = (hi-lo)[:, None]*nodes/2+(hi+lo)[:, None]/2
    t = np.log(x/u[i, None])/np.log(u[i+1]/u[i])[:, None]
    lk = np.log(data[:, 2:4])
    k = np.exp(lk[i, None, :]+t[:, :, None]*(lk[i+1]-lk[i])[:, None, :])
    n = np.sqrt(index_squared(x, cutoff, vstar2))
    wr = x**4*np.exp(-x)/(-np.expm1(-x))**2
    # Keep scattering unchanged and absorption divided by n in all variants.
    inv = n**3/(k[:, :, 0]+n*k[:, :, 1])
    return RW/np.sum((hi-lo)[:, None]/2*weights*wr*inv)


def controls():
    cold = electron_moments(1000, 3e27)
    cold_error = abs(cold['plasma_frequency_squared_ratio']/cold['zero_temperature_frequency_squared_ratio']-1)
    dilute = electron_moments(1e6, 1e20)
    dilute_error = abs(dilute['plasma_frequency_squared_ratio']-(1-2.5*dilute['theta']))
    x = np.geomspace(1+1e-8, 100, 128)
    nonrel_error = float(np.max(np.abs(index_squared(x, 1., 0.)-(1-x**-2))))
    z = np.linspace(.01, .45, 100)
    exact = 3/(2*z)*(1-(1-z)*np.arctanh(np.sqrt(z))/np.sqrt(z))
    series_error = float(np.max(np.abs(transverse_factor(z)[0]/exact-1)))
    if cold_error > 1e-8 or dilute_error > 2e-6 or nonrel_error > 1e-14 or series_error > 1e-10:
        raise ValueError('independent electron/dispersion limit checks failed')
    return {'zero_temperature_relative_error': cold_error,
            'dilute_first_order_absolute_error': dilute_error,
            'nonrelativistic_index_absolute_error': nonrel_error,
            'series_vs_logarithmic_relative_error': series_error}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['work', 'classical_check', 'component_check', 'output']:
        p.add_argument(name, type=Path)
    a = p.parse_args()
    checks = controls()
    classical = json.loads(a.classical_check.read_text())
    components = json.loads(a.component_check.read_text())
    inputs = {str(v.resolve()): digest(v) for v in [Path(__file__), a.classical_check,
              a.component_check, Path(__file__).with_name('audit_tops_spectral_means.py')]}
    for report in [classical, components]:
        for path, checksum in report['input_sha256'].items():
            if digest(Path(path)) != checksum:
                raise ValueError('control input changed')
            inputs[path] = checksum
    spectra, _ = source(a.work/'cutoff-on')
    rows = []
    for row in classical['records']:
        key = next(k for k in spectra if k[1] == row['density'] and
                   abs(k[0]*KEV/KB/row['temperature_K']-1) < 1e-12)
        electron = electron_moments(row['temperature_K'], row['estimated_electron_density'])
        u0 = row['classical_cutoff_u']
        up = u0*math.sqrt(electron['plasma_frequency_squared_ratio'])
        recipes = {'classical': (u0, 0), 'finite_temperature_frequency': (up, 0),
                   'transverse_dispersion': (up, electron['vstar_squared'])}
        values, errors = {}, {}
        for name, (cutoff, velocity) in recipes.items():
            pair = [mean(spectra[key], key[0], cutoff, velocity, order) for order in [8, 16]]
            values[name], errors[name] = pair[1], abs(pair[0]/pair[1]-1)
        control = next(r for r in components['records'] if r['temperature_keV'] == key[0] and r['density'] == key[1])
        baseline = control['radiative_opacities']['absorption_divided_by_n_scattering_unchanged']
        if abs(values['classical']/baseline-1) > 1e-12 or max(errors.values()) > 1e-7:
            raise ValueError('spectral quadrature or independent baseline check failed')
        kc = row['conductive_opacity']
        combined = {name: 1/(1/value+1/kc) for name, value in values.items()}
        rows.append({'temperature_keV': key[0], 'density': key[1], 'electron': electron,
                     'classical_cutoff_u': u0, 'finite_temperature_cutoff_u': up,
                     'radiative_opacities': values, 'combined_opacities': combined,
                     'relative_to_classical': {name: {'radiative': values[name]/values['classical']-1,
                     'combined': combined[name]/combined['classical']-1} for name in recipes},
                     'quadrature_relative_errors': errors})
    if any(digest(Path(path)) != checksum for path, checksum in inputs.items()):
        raise ValueError('source inputs changed during calculation')
    result = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
              'reference': 'https://arxiv.org/abs/hep-ph/9302213',
              'input_sha256': inputs, 'limit_checks': checks, 'records': rows}
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps([{'T_keV': r['temperature_keV'], 'density': r['density'],
                      'changes': r['relative_to_classical']} for r in rows]), flush=True)


if __name__ == '__main__':
    main()
