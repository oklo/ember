#!/usr/bin/env python3
"""Continuous mixture corrections for the archived OP atomic spectra.

This is a conditional opacity diagnostic, not a selected stellar input.
It evaluates the long-wavelength, static ring approximation to electron
scattering and the screened Born/Elwert free-free subtraction used by
OPserver, retaining its numerical unit constants. Electron Fermi integrals,
the scattering-angle integral and the free-free thermal integral are evaluated
numerically instead of copying the supplied fits and coarse frequency mesh.

Ion populations and electron counts are per the SAME reference particle;
they need not sum to unity. Returned cross sections are per that reference,
in a0^2, WITHOUT stimulated emission, matching the native OP files.
"""
from functools import lru_cache
import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expit, roots_laguerre, roots_legendre

THOMSON_ATOMIC = 2.37567e-8


def _state(temperature, electron_density, ions, electrons):
    ions = np.asarray(ions, dtype=float)
    if (ions.shape != (28,) or not np.isfinite(ions).all() or
            np.any(ions < 0) or not math.isfinite(temperature) or
            not math.isfinite(electron_density) or
            not math.isfinite(electrons) or
            min(temperature, electron_density, electrons) <= 0):
        raise ValueError('invalid thermodynamic state or positive-ion populations')
    charges = np.arange(1., 29.)
    first = float(ions @ charges)
    second = float(ions @ (charges**2))
    if first <= 0:
        raise ValueError('no positive ions for the mixture correction')
    return ions, first, second


def fermi_integrals(eta):
    """Normalized I_(1/2), I_(-1/2), using p=sqrt(energy/kT)."""
    if not -60 <= eta <= 100:
        raise ValueError('electron degeneracy outside diagnostic range')
    end = math.sqrt(max(eta, 0.) + 65.)
    ph = 4/math.sqrt(math.pi)*quad(
        lambda p:p*p*expit(eta-p*p), 0., end,
        epsabs=1e-32, epsrel=3e-12, limit=150)[0]
    mh = 2/math.sqrt(math.pi)*quad(
        lambda p:expit(eta-p*p), 0., end,
        epsabs=1e-32, epsrel=3e-12, limit=150)[0]
    return ph, mh


@lru_cache(maxsize=512)
def electron_response(temperature, electron_density):
    """OP unit convention, but directly integrated Fermi response."""
    aune = 1.48185e-25*electron_density
    aut = 3.16668e-6*temperature
    target = 7.8748*aune/aut**1.5
    eta = brentq(lambda x:fermi_integrals(x)[0]/target-1., -60., 100.,
                 xtol=2e-12)
    ph, mh = fermi_integrals(eta)
    return dict(eta=eta, normalized_density=target, compressibility_ratio=mh/ph)


@lru_cache(maxsize=12)
def angular_mesh(order):
    nodes, weights = roots_legendre(order)
    return 1-nodes, weights*(1+nodes*nodes)*(1-nodes)


def transport_shape(d, order=64):
    """Integral (1+mu^2)(1-mu) D/(D+1-mu) dmu, without cancellation."""
    d = np.asarray(d, dtype=float)
    if not np.isfinite(d).all() or np.any(d < 0):
        raise ValueError('invalid scattering screening coordinate')
    y, weight = angular_mesh(order)
    # For small D the closed expression is stable and resolves the logarithm
    # at the forward-scattering endpoint better than fixed angular nodes.
    flat = d.reshape(-1)
    answer = np.empty_like(flat)
    small = flat <= .25
    a = flat[small]
    nz = a > 0
    values = np.zeros_like(a)
    b = a[nz]
    values[nz] = b*(8/3+2*b*(1+b)-b*(b*b+2*b+2)*np.log1p(2/b))
    answer[small] = values
    b = flat[~small]
    answer[~small] = np.sum(weight*b[:, None]/(b[:, None]+y), axis=1)
    return answer.reshape(d.shape)


def scattering_factor(temperature, electron_density, ions, u, order=64,
                      relativistic=True):
    ions = np.asarray(ions, dtype=float)
    ions, first, second = _state(temperature, electron_density, ions, 1.)
    u = np.asarray(u, dtype=float)
    if not np.isfinite(u).all() or np.any(u <= 0):
        raise ValueError('invalid photon frequency')
    response = electron_response(temperature, electron_density)
    r = response['compressibility_ratio']
    x = r + second/first
    aune = 1.48185e-25*electron_density
    aut = 3.16668e-6*temperature
    d = (1.1799e5*x*aune/aut**3)/(u*u)
    # Exchange at zero momentum gives 1-R. The ring contribution is the
    # transport-angle average from Boercker (1987), eqs. 3 and 11.
    factor = r - (3/8)*r*r/x*transport_shape(d, order)
    if relativistic:
        v = aut*u
        factor *= 1+v*(-1.0650e-4+v*(1.4746e-8-v*2.0084e-12))
    if not np.isfinite(factor).all() or np.any(factor <= 0):
        raise ValueError('nonpositive scattering factor outside model support')
    return factor


def screened_momentum_difference(alpha_squared, e, photon):
    """Integral of screened minus Coulomb Born momentum transfer.

    Energies are in Rydbergs. q_minus=photon/q_plus avoids cancellation
    between two square roots for small photon energy.
    """
    plus = np.sqrt(e+photon)+np.sqrt(e)
    minus = photon/plus
    a = alpha_squared/(plus*plus)
    b = alpha_squared/(minus*minus)
    return np.log1p(a)-np.log1p(b)+a/(1+a)-b/(1+b)


@lru_cache(maxsize=16)
def thermal_mesh(order):
    return roots_laguerre(order)


def freefree_subtraction(temperature, electron_density, ions, electrons, u,
                         order=64, minimum_ion_fraction=0., atom_count=None):
    """Screened Born/Elwert subtraction, continuous at every frequency.

    Retains the OP recipe's ionic Debye length and classical thermal weight.
    No low-screening skip, endpoint overwrite, or negative-opacity clipping.
    A nonzero ion threshold is only for comparison with the supplied recipe.
    """
    ions, first, second = _state(temperature, electron_density, ions, electrons)
    u = np.asarray(u, dtype=float)
    if u.ndim != 1 or not np.isfinite(u).all() or np.any(u <= 0):
        raise ValueError('invalid photon frequency vector')
    if minimum_ion_fraction:
        if atom_count is None or atom_count <= 0:
            raise ValueError('ion threshold needs the reference atom count')
        use = ions/atom_count > minimum_ion_fraction
    else:
        use = ions > 0
    rydt = temperature/157894.
    aune = 1.48185e-25*electron_density
    alpha_squared = 5.8804e-19*electron_density*second/(electrons*temperature)
    x, weights = thermal_mesh(order)
    e = (x*rydt)[None, :]
    photon = u[:, None]*rydt
    lower = np.sqrt(e)
    upper = np.sqrt(e+photon)
    momentum = screened_momentum_difference(alpha_squared, e, photon)
    result = np.zeros_like(u)
    for i in np.flatnonzero(use):
        z = i+1
        elwert = (upper/lower)*(-np.expm1(-6.283185*z/upper))/(-np.expm1(-6.283185*z/lower))
        result += ions[i]*z*z*((momentum*elwert) @ weights)
    result *= (1.7337*aune/math.sqrt(rydt))/(u*rydt)**3
    if not np.isfinite(result).all() or np.any(result > 0):
        raise ValueError('invalid screened free-free subtraction')
    return result


def corrections(temperature, electron_density, ions, electrons, u, order=64):
    """Components in native OP units; positivity of the sum is caller-owned."""
    factor = scattering_factor(temperature, electron_density, ions, u)
    se = -np.expm1(-np.asarray(u))
    original_scattering = electrons*THOMSON_ATOMIC/se
    scattering = original_scattering*factor
    freefree = freefree_subtraction(temperature, electron_density, ions,
                                   electrons, u, order)
    return dict(scattering=scattering,
                scattering_change=scattering-original_scattering,
                freefree_change=freefree,
                electron_response=electron_response(temperature, electron_density))
