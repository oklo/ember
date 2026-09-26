"""Finite-temperature electron momentum drag for a stated comparison model.

Nonrelativistic Fermi electrons scatter elastically from stationary point ions
through V(r) = -Z e^2 exp(-r/lambda)/r, in the first Born approximation.
The electron perturbation is a displaced Fermi distribution, with no heat
flow. This defines a momentum moment, not a complete kinetic solution or
an electrical/thermal conductivity at arbitrary degeneracy.

At zero temperature the result reduces to Potekhin et al. (1999), equations
7--8, with nonrelativistic kinematics, S(q)=1 and Yukawa screening:
https://www.ioffe.ru/astro/Stars/Paper/pbhy99.pdf .
The finite-temperature drift moment is derived in docs/ELECTRON_ION_DRAG.md.
No prescription for screening, ion correlations or stellar ionization is
selected by this module. Do not use ion collision fits for these electrons.
"""
import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expit

KB = 1.380649e-16
ME = 9.1093837015e-28
HBAR = 1.054571817e-27
E2 = 2.3070775523417355e-19
CLIGHT = 2.99792458e10


def coulomb_bracket(b):
    """ln(1+b)-b/(1+b), retaining its small-b limit without cancellation."""
    b = np.asarray(b, dtype=float)
    if not np.all(np.isfinite(b)) or np.any(b < 0):
        raise ValueError('finite nonnegative squared momentum required')
    small = b < 1e-3
    # Alternating expansion through b^8; omitted relative term < 2e-21.
    series = b*b*(.5+b*(-2/3+b*(3/4+b*(-4/5+b*(5/6+b*(-6/7+b*7/8))))))
    return np.where(small, series, np.log1p(b)-b/(1+b))


def cross_section(momentum, charge, screening_length):
    """Born momentum-transfer cross section in cm^2, including p=0."""
    p = np.asarray(momentum, dtype=float)
    if (not np.all(np.isfinite(p)) or np.any(p < 0) or
            not math.isfinite(charge) or charge <= 0 or
            not math.isfinite(screening_length) or screening_length <= 0):
        raise ValueError('nonnegative momentum, positive ion charge and length required')
    b = (2*p*screening_length/HBAR)**2
    small = b < 1e-3
    series = .5+b*(-2/3+b*(3/4+b*(-4/5+b*(5/6+b*(-6/7+b*7/8)))))
    ratio = np.divide(coulomb_bracket(b), b*b, out=np.zeros_like(b), where=b != 0)
    ratio = np.where(small, series, ratio)
    return 32*math.pi*ME**2*(charge*E2)**2*screening_length**4/HBAR**4*ratio


def fermi_half(eta):
    """Unnormalized integral int_0^infinity sqrt(x)/(exp(x-eta)+1) dx."""
    if not math.isfinite(eta):
        raise ValueError('finite chemical potential required')
    if eta > 20:
        # Subtract the filled Fermi sphere before integrating its edge.
        def edge(y):
            return (math.sqrt(eta+y)-math.sqrt(max(eta-y, 0)))*expit(-y)
        points = [eta] if eta < 50 else None
        correction = quad(edge, 0, 50, points=points, epsabs=1e-13,
                          epsrel=2e-12, limit=160)[0]
        return 2/3*eta**1.5+correction
    # x=t^2 removes the square-root endpoint.
    upper = math.sqrt(max(eta, 0)+50)
    return quad(lambda t: 2*t*t*expit(eta-t*t), 0, upper,
                epsabs=1e-100, epsrel=2e-12, limit=160)[0]


def chemical_potential(number_density, temperature):
    """Solve the nonrelativistic ideal Fermi number equation; return mu/kT."""
    if (not math.isfinite(number_density) or number_density <= 0 or
            not math.isfinite(temperature) or temperature <= 0):
        raise ValueError('positive electron density and temperature required')
    scale = (2*ME*KB*temperature)**1.5/(2*math.pi**2*HBAR**3)
    target = number_density/scale
    classical = math.log(target/(math.sqrt(math.pi)/2))
    degenerate = (1.5*target)**(2/3)
    lo, hi = min(-40., classical-5), max(40., degenerate+10)
    return brentq(lambda eta: math.log(fermi_half(eta)/target), lo, hi,
                  xtol=2e-12, rtol=2e-14)


def drift_integral(eta, b_thermal, *, tail=50., epsrel=2e-11):
    """int_0^infinity L(b_thermal*x) f(x)[1-f(x)] dx."""
    if (not math.isfinite(eta) or not math.isfinite(b_thermal) or
            b_thermal <= 0 or not math.isfinite(tail) or tail < 30):
        raise ValueError('invalid chemical potential, screening or tail')
    # Center the integral on the Fermi edge; this also resolves eta >> 1.
    lo, hi = max(-eta, -tail), max(-eta, 0)+tail
    def integrand(y):
        return float(coulomb_bracket(b_thermal*(eta+y)))*expit(y)*expit(-y)
    points = [0.] if lo < 0 < hi else None
    return quad(integrand, lo, hi, points=points, epsabs=1e-100,
                epsrel=epsrel, limit=180)[0]


def pair_resistance(*, electron_density, ion_density, temperature, charge,
                    screening_length, eta=None):
    """Return displaced-Fermi electron-ion drag K_ei in g/(cm^3 s).

    eta is optionally a previously solved mu/kT at the same ne,T; when supplied
    its number equation is checked. This allows one electron solve per mixture.
    Relativity, moving ions, ion correlations and electron heat perturbations
    are absent and must be assessed by callers before physical selection.
    """
    if (not math.isfinite(ion_density) or ion_density <= 0 or
            not math.isfinite(charge) or charge <= 0 or
            not math.isfinite(screening_length) or screening_length <= 0):
        raise ValueError('positive ion density, charge and screening required')
    if eta is None:
        eta = chemical_potential(electron_density, temperature)
    elif (not math.isfinite(electron_density) or electron_density <= 0 or
          not math.isfinite(temperature) or temperature <= 0 or
          abs((2*ME*KB*temperature)**1.5/(2*math.pi**2*HBAR**3)
              *fermi_half(eta)/electron_density-1) > 1e-9):
        raise ValueError('chemical potential does not match the electron state')
    bt = 8*ME*KB*temperature*screening_length**2/HBAR**2
    moment = drift_integral(eta, bt)
    coefficient = 2*ME**2*(charge*E2)**2/(3*math.pi*HBAR**3)*moment
    pf = HBAR*(3*math.pi**2*electron_density)**(1/3)
    speed = math.sqrt((pf/ME)**2+2*KB*temperature/ME)
    return dict(resistance_g_cm3_s=ion_density*coefficient,
                force_per_ion_per_drift_g_s=coefficient,
                eta_nonrelativistic=eta, dimensionless_moment=moment,
                fermi_momentum_over_mec=pf/(ME*CLIGHT),
                characteristic_Born_parameter=charge*E2/(HBAR*speed),
                static_potential_strength=2*ME*charge*E2*screening_length/HBAR**2)
