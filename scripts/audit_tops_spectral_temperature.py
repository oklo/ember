#!/usr/bin/env python3
"""Withhold a native temperature and compare spectral and grey interpolation.

Atomic absorption and scattering are interpolated at fixed photon energy.
Electron dispersion and the full Rosseland weight are evaluated at the query
temperature. The withheld spectrum is used only as an independent reference.
This is an interpolation experiment, not an installed opacity prescription.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from audit_tops_spectral_means import source, digest
from audit_tops_electron_dispersion import electron_moments, index_squared, KEV, KB, HBAR, ME, RW


def temperature_interpolate(temperatures, values, target, method):
    """The runtime's arithmetic-slope monotone Hermite, vectorized in columns."""
    x = np.log(np.asarray(temperatures))
    q = math.log(target)
    if not x[0] <= q <= x[-1]:
        raise ValueError('temperature extrapolation')
    k = min(max(int(np.searchsorted(x, q, side='right'))-1, 0), len(x)-2)
    h = x[k+1]-x[k]; t = (q-x[k])/h
    if method == 'linear':
        return (1-t)*values[k]+t*values[k+1]
    if method != 'hermite':
        raise ValueError('unknown interpolation method')

    def secant(i):
        return (values[i+1]-values[i])/(x[i+1]-x[i])

    def limited(a, b):
        s = .5*(a+b)
        small = np.where(np.abs(a) <= np.abs(b), a, b)
        return np.where(a*b <= 0., 0., np.where(np.abs(s) > 3*np.abs(small), 3*small, s))

    dc = secant(k)
    s0 = limited(secant(k-1), dc) if k else dc
    s1 = limited(dc, secant(k+1)) if k+2 < len(x) else dc
    t2, t3 = t*t, t*t*t
    return ((2*t3-3*t2+1)*values[k]+(t3-2*t2+t)*h*s0
            +(-2*t3+3*t2)*values[k+1]+(t3-t2)*h*s1)


def mean_from_interpolated_spectra(spectra, records, tt, rho, target, method, order):
    arrays = [spectra[t, rho] for t in tt]
    ne = math.exp(float(temperature_interpolate(tt,
        np.log([records[t, rho]['free_electron_density_cm3'] for t in tt]), target, method)))
    electron = electron_moments(target*KEV/KB, ne)
    cutoff = HBAR*math.sqrt(4*math.pi*4.80320471257e-10**2*ne/ME)/(target*KEV)
    cutoff *= math.sqrt(electron['plasma_frequency_squared_ratio'])
    low, high = cutoff*target, 700*target
    if max(v[0, 0] for v in arrays) > low or min(v[-1, 0] for v in arrays) < high:
        raise ValueError('a source lacks common photon-energy coverage')
    # The union uses the interpolation sources only, not the withheld mesh.
    energies = np.unique(np.concatenate([v[:, 0] for v in arrays]))
    energies = np.concatenate(([low], energies[(energies > low) & (energies < high)], [high]))
    lo, hi = energies[:-1]/target, energies[1:]/target
    nodes, weights = leggauss(order)
    u = ((hi-lo)[:, None]*nodes/2+(hi+lo)[:, None]/2).ravel()
    photon_energy = u*target
    log_energy = np.log(photon_energy)
    log_components = np.stack([np.column_stack([
        np.interp(log_energy, np.log(v[:, 0]), np.log(v[:, c])) for c in (2, 3)]) for v in arrays])
    components = np.exp(temperature_interpolate(tt, log_components, target, method))
    n = np.sqrt(index_squared(u, cutoff, electron['vstar_squared']))
    wr = u**4*np.exp(-u)/(-np.expm1(-u))**2
    f = wr*n**3/(components[:, 0]+n*components[:, 1])
    integral = float(np.sum((hi-lo)[:, None]/2*weights*f.reshape(len(lo), order)))
    return {'rosseland': RW/integral, 'electron_density': ne, 'cutoff_u': cutoff,
            'frequency_intervals': len(lo)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('report', 'plan', 'scratch', 'output'):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    if a.output.exists() or a.scratch.exists():
        raise FileExistsError('use new experiment paths')
    a.scratch.mkdir(parents=True)
    report = json.loads(a.report.read_text()); plan = json.loads(a.plan.read_text())
    inputs = dict(report['input_sha256'])
    inputs.update({str(path.resolve()): digest(path) for path in (a.report, a.plan, Path(__file__))})
    for name, h in inputs.items():
        if digest(Path(name)) != h:
            raise ValueError('source dependency changed: '+name)
    records = {(r['temperature_keV'], r['density_atomic_g_cm3']): r for r in report['records']}
    target = plan['temperature_check_keV']
    tt = [t for t in plan['node_temperatures_keV'] if t != target]
    spectra = {}
    for directory in sorted({r['source'] for r in records.values()}):
        data, _ = source(Path(directory))
        for key, value in data.items():
            if key in spectra:
                raise ValueError('duplicate source spectrum')
            spectra[key] = value
    # A log power law has an affine log-temperature dependence and must be exact.
    trial = np.log(np.asarray(tt))[:, None]*np.array([-.5, 0., 2., 3.5])[None, :]+1.2
    analytic_error = float(np.max(np.abs(temperature_interpolate(tt, trial, target, 'hermite')
                                        -(1.2+math.log(target)*np.array([-.5, 0., 2., 3.5])))))
    if analytic_error > 1e-13:
        raise ValueError('power-law interpolation control failed')
    output = []
    for rho in plan['node_densities_g_cm3']:
        reference = records[target, rho]['rosseland_atomic_cm2_g']
        grey = math.exp(float(temperature_interpolate(tt, np.log([
            records[t, rho]['rosseland_atomic_cm2_g'] for t in tt]), target, 'hermite')))
        old = next(r for r in report['comparisons'] if r['kind'] == 'temperature'
                   and r['density_atomic_g_cm3'] == rho)
        runtime_difference = abs(grey/old['runtime_opacity']-1)
        if runtime_difference > 1e-12:
            raise ValueError('vector interpolation differs from the compiled runtime')
        methods = {}
        for method in ('linear', 'hermite'):
            pair = [mean_from_interpolated_spectra(spectra, records, tt, rho, target, method, order)
                    for order in (8, 16)]
            numerical = abs(pair[0]['rosseland']/pair[1]['rosseland']-1)
            if numerical > 1e-6:
                raise ValueError('spectral quadrature did not converge')
            methods[method] = {**pair[1], 'quadrature_relative_change': numerical,
                              'relative_error': pair[1]['rosseland']/reference-1,
                              'electron_density_relative_error': pair[1]['electron_density']/records[target, rho]['free_electron_density_cm3']-1}
        output.append({'density_atomic_g_cm3': rho, 'temperature_keV': target,
                       'reference_rosseland': reference, 'grey_hermite_relative_error': grey/reference-1,
                       'runtime_control_relative_difference': runtime_difference, 'methods': methods})
        (a.scratch/'partial.json').write_text(json.dumps(output, indent=2)+'\n')
        print(json.dumps({'density': rho, 'errors': {name: v['relative_error'] for name, v in methods.items()}}), flush=True)
    if any(digest(Path(name)) != h for name, h in inputs.items()):
        raise ValueError('input changed during comparison')
    result = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
              'power_law_control_absolute_error': analytic_error, 'records': output,
              'maximum_relative_errors': {method: max(abs(r['methods'][method]['relative_error']) for r in output)
                                          for method in ('linear', 'hermite')},
              'radiative_interpolation_relative_criterion': .005,
              'assumptions': ['The withheld temperature supplies reference values only.',
                  'Atomic absorption and scattering are interpolated at fixed photon energy, avoiding artificial movement of bound edges with temperature.',
                  'Interpolated free-electron density enters the independently evaluated electron dispersion at the target temperature.',
                  'Passing a coarsened native-temperature comparison would remain a sampled check, not a global physical error bound.'],
              'input_sha256': inputs}
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result['maximum_relative_errors']), flush=True)


if __name__ == '__main__':
    main()
