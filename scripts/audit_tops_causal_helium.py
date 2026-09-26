"""Test a causal, sum-conserving low-frequency continuation of pure-He spectra.

The join-frequency variation is a physical-model sensitivity, not numerical
convergence. The continued absorption replaces the low-frequency source;
it never multiplies an already corrected spectrum by a second damping factor.
These calculations are diagnostics and do not select an opacity table.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad

from audit_tops_spectral_means import source, digest
from audit_tops_conductivity_weight import integrate
from causal_conductivity import ContinuedConductivity, real_refractive_index

C, HBAR, KEV, KB = 2.99792458e10, 1.054571817e-27, 1.602176634e-9, 1.380649e-16
MU, ME, E, SIGMA_SB = 1.66053906660e-24, 9.1093837015e-28, 4.80320471257e-10, 5.670374419e-5
HE_MASS = 4.002602
ROSS_NORM = 4*math.pi**4/15


def analytic_controls():
    u = np.geomspace(1e-7, 1e5, 16001)
    D, plasma, collision = 1.7, 2., .3
    opacity = plasma**2*collision/(D*(u*u+collision**2))
    model = ContinuedConductivity(u, opacity, D, plasma, .5, tail_power=2)
    points = np.array([1e-5, .01, .1, .49, .5, .51, 1., 2., 3., 10., 100.])
    real, imaginary = model.dielectric(points)
    expected_real = 1-plasma**2/(points**2+collision**2)
    expected_imaginary = plasma**2*collision/(points*(points**2+collision**2))
    errors = {'real_scaled_absolute': float(np.max(abs(real-expected_real)/(1+abs(expected_real)))),
              'imaginary_relative': float(np.max(abs(imaginary/expected_imaginary-1))),
              'collision_relative': abs(model.collision_u/collision-1),
              'sum_relative': abs(model.sum_relative_error)}
    if max(errors.values()) > 2e-5:
        raise ValueError('analytic Drude control failed')
    # Independent numerical integration of the subtracted transform also tests
    # the p=3 tail, which differs from the Drude p=2 asymptote.
    sample = np.geomspace(1e-4, 50., 6001)
    other = ContinuedConductivity(sample, 1/(1+sample**2), 1., 2., .8)
    quadrature = []
    for q in (.23, 1.17, 3.1):
        kq = float(other.absorption(q))
        def f(x):
            if abs(x-q) < 1e-8:
                h = 1e-6*q
                derivative = (other.absorption(q+h)-other.absorption(q-h))/(2*h)
                return float(derivative/(2*q))
            kx = other.amplitude/(x*x+other.collision_u**2) if x == 0 else other.absorption(x)
            return float((kx-kq)/(x*x-q*q))
        bounds = sorted({0., other.join, q, 1., 5., 20., float(other.x[-1])})
        direct = sum(quad(f, a, b, epsabs=2e-6, epsrel=2e-6, limit=500)[0]
                     for a, b in zip(bounds, bounds[1:]))
        b = other.x[-1]
        direct += quad(lambda t: (other.k[-1]*t**3-kq)/(b*(1-(q*t/b)**2)),
                       0, 1, epsabs=1e-11)[0]
        expected = 1+2*other.D/math.pi*direct
        actual = float(other.dielectric(np.array([q]))[0][0])
        quadrature.append(abs(actual-expected)/(1+abs(expected)))
    if max(quadrature) > 1e-5:
        raise ValueError('independent principal-value quadrature failed')
    errors['independent_transform_scaled_absolute'] = max(quadrature)
    return errors


def rosseland(model, data, temperature, intervals):
    edges = np.unique(np.r_[0., np.geomspace(1e-8, 100, intervals+1),
                            model.join, model.collision_u, model.plasma_u])
    nodes, weights = leggauss(8)
    width = np.diff(edges)
    u = ((edges[1:]+edges[:-1])[:, None]+width[:, None]*nodes)/2
    measure = width[:, None]*weights/2
    epsilon1, epsilon2 = model.dielectric(u)
    n = real_refractive_index(epsilon1, epsilon2)
    k = model.absorption(u)
    source_u = data[:, 0]/temperature
    scatter = np.exp(np.interp(np.log(u), np.log(source_u), np.log(data[:, 3])))
    weight = u**4*np.exp(-u)/(-np.expm1(-u))**2
    transport = measure*weight*n**3/(k+n*scatter)
    inverse = float(np.sum(transport))
    no_low_scatter = np.where(u < source_u[0], 0., scatter)
    low_control = float(np.sum(measure*weight*n**3/(k+n*no_low_scatter)))
    original = np.exp(np.interp(np.log(u), np.log(source_u), np.log(data[:, 2])))
    classical = np.sqrt(np.maximum(0., 1-(model.plasma_u/u)**2))
    original_inverse = float(np.sum(measure*weight*classical**3/(original+classical*scatter)))
    if not inverse > 0 or not original_inverse > 0:
        raise ValueError('nonpositive Rosseland integral')
    return {'rosseland': ROSS_NORM/inverse,
            'classical_component_rosseland': ROSS_NORM/original_inverse,
            'inverse_rosseland_fraction_negative_epsilon1': float(np.sum(transport[epsilon1 <= 0]))/inverse,
            'lowest_frequency_scattering_relative_effect': inverse/low_control-1}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('work', 'conduction_probe', 'scratch', 'output'):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    if a.output.exists() or a.scratch.exists():
        raise FileExistsError('use fresh diagnostic output and scratch paths')
    a.scratch.mkdir(parents=True)
    controls = analytic_controls()
    print(json.dumps({'analytic_controls': controls}), flush=True)
    paths = [a.work/('cutoff-'+label)/name for label in ('on', 'off')
             for name in ('receipt.json', 'request.json', 'source.txt')]
    paths += [Path(__file__), a.conduction_probe, Path('data/conduction/condtab21wd_metals.dat')]
    paths.append(Path(__file__).with_name('conduction_composition_probe.cpp'))
    paths += [Path(__file__).with_name(n) for n in ('causal_conductivity.py',
              'audit_tops_spectral_means.py', 'audit_tops_conductivity_weight.py')]
    inputs = {str(v.resolve()): digest(v) for v in paths}
    for label in ('on', 'off'):
        receipt = json.loads((a.work/('cutoff-'+label)/'receipt.json').read_text())
        if receipt['X'] != 0 or receipt['Z'] != 0 or receipt['metals'] != {}:
            raise ValueError('this diagnostic requires pure helium')
    (spectra, means), (off, _) = [source(a.work/('cutoff-'+v)) for v in ('on', 'off')]
    if set(spectra) != set(off) or any(not np.array_equal(v, off[k]) for k, v in spectra.items()):
        raise ValueError('source plasma option changed a spectrum')
    rows = []
    keys = sorted(spectra)
    for temperature, rho in keys:
        T, rho_b = temperature*KEV/KB, rho*4/HE_MASS
        rows.append(f'1 1 {rho_b:.17g} {T:.17g} 1 0 0 1\n')
    request = a.scratch/'conduction.txt'; request.write_text(''.join(rows))
    response = subprocess.run([str(a.conduction_probe.resolve()),
                              str(Path('data/conduction/condtab21wd_metals.dat').resolve()),
                              str(request.resolve()), '0'], capture_output=True, text=True)
    (a.scratch/'conduction.out').write_text(response.stdout)
    (a.scratch/'conduction.err').write_text(response.stderr)
    response.check_returncode()
    lines = response.stdout.splitlines()
    if len(lines) != len(keys)+1 or not lines[-1].startswith('checksum '):
        raise ValueError('unexpected conduction response')
    records = []
    for (temperature, rho), line in zip(keys, lines[:-1], strict=True):
        if means[temperature, rho]['free_electrons_per_ion'] != 2:
            raise ValueError('helium control is not fully ionized')
        data = spectra[temperature, rho]
        T, scale = temperature*KEV/KB, temperature*KEV/HBAR
        ne = 2*rho/(HE_MASS*MU)
        wp2 = 4*math.pi*E**2*ne/ME
        plasma = math.sqrt(wp2)/scale
        D = C*rho/scale
        kcond_b = float(line.split()[0])
        if not kcond_b > 0:raise ValueError('nonpositive conduction opacity')
        mass_scale, rho_b = HE_MASS/4, rho*4/HE_MASS
        kcond_atomic = kcond_b/mass_scale
        thermal_conductivity = 16*SIGMA_SB*T**3/(3*rho_b*kcond_b)
        lorenz = math.pi**2*KB**2/(3*E**2)
        sigma_dc = thermal_conductivity/(lorenz*T)
        collision_WF = wp2/(4*math.pi*sigma_dc)/scale
        measured_area, _, _ = integrate(data[:, 0]/temperature, data[:, 2])
        variants = []
        for fraction in (.5, .75, 1.):
            calculations = []
            for mesh, refinement in ((256, 1), (512, 1), (512, 2)):
                model = ContinuedConductivity(data[:, 0]/temperature, data[:, 2],
                                              D, plasma, fraction*plasma, refinement=refinement)
                value = rosseland(model, data, temperature, mesh)
                calculations.append({**value, 'integration_intervals': mesh,
                                     'spectral_refinement': refinement})
            final = calculations[-1]
            combined = lambda k: 1/(1/k+1/kcond_atomic)
            base = final['classical_component_rosseland']
            variants.append({'join_over_plasma': fraction, 'join_u': model.join,
                'collision_u': model.collision_u, 'collision_over_Wiedemann_Franz': model.collision_u/collision_WF,
                'conductivity_sum_relative_error': model.sum_relative_error,
                'high_frequency_tail_weight_fraction': model.tail_area/model.target_area,
                'rosseland_atomic': final['rosseland'],
                'classical_component_rosseland_atomic': base,
                'combined_relative_change_from_classical': combined(final['rosseland'])/combined(base)-1,
                'combined_opacity_baryonic': combined(final['rosseland'])*mass_scale,
                'frequency_integration_relative_change': calculations[1]['rosseland']/calculations[0]['rosseland']-1,
                'spectral_transform_refinement_relative_change': final['rosseland']/calculations[1]['rosseland']-1,
                'lowest_frequency_scattering_relative_effect': final['lowest_frequency_scattering_relative_effect'],
                'inverse_rosseland_fraction_negative_epsilon1': final['inverse_rosseland_fraction_negative_epsilon1'],
                'calculations': calculations})
            print(json.dumps({'T_keV': temperature, 'rho_atomic': rho, **variants[-1]}), flush=True)
        records.append({'temperature_keV': temperature, 'density_atomic': rho,
                        'density_baryonic': rho_b, 'free_electron_density': ne,
                        'plasma_u': plasma, 'conductive_opacity_baryonic': kcond_b,
                        'thermal_conductivity_cgs': thermal_conductivity,
                        'Wiedemann_Franz_collision_u': collision_WF,
                        'measured_conductivity_sum_fraction': measured_area/(math.pi*plasma**2/(2*D)),
                        'variants': variants})
        (a.scratch/'partial.json').write_text(json.dumps(records, indent=2, allow_nan=False)+'\n')
    numerical = max(abs(v[k]) for r in records for v in r['variants'] for k in
                    ('frequency_integration_relative_change', 'spectral_transform_refinement_relative_change'))
    for name, h in inputs.items():
        if digest(Path(name)) != h:raise ValueError('source changed during the calculation')
    report = {'scope': __doc__, 'accepted_for_stellar_opacity': False,
              'analytic_controls': controls, 'records': records,
              'maximum_numerical_relative_change': numerical,
              'numerical_relative_criterion': 1e-4, 'numerical_checks_passed': numerical <= 1e-4,
              'assumptions': ['Nonrelativistic conductivity sum, fully ionized pure helium.',
                  'Positive index-one absorption is proportional to optical conductivity; its low-frequency part is replaced, not damped twice.',
                  'A continuous Drude low-frequency part, retained high-frequency absorption, and a power-three far tail define one causal response.',
                  'Join-frequency variation is physical approximation sensitivity, not an error bound.',
                  'Rosseland integration uses n_real^3/(kappa_abs+n_real*kappa_scatter) with full thermal normalization; deeply damped modes require a separate transport-validity assessment.',
                  'Wiedemann-Franz uses the degenerate Lorenz number as a diagnostic only; finite-temperature and electron-electron corrections are not inferred from thermal conductivity alone.',
                  'Atomic-mass opacities/densities are converted explicitly for comparison with baryonic conduction.'],
              'references': ['https://doi.org/10.1016/j.hedp.2017.02.008',
                             'https://doi.org/10.1016/j.hedp.2010.01.004'],
              'input_sha256': inputs}
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'report': str(a.output), 'maximum_numerical_relative_change': numerical,
                      'numerical_checks_passed': numerical <= 1e-4}), flush=True)


if __name__ == '__main__':
    main()
