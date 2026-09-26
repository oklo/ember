"""Eliminate electron energy dependence for a common elastic scattering kernel.

Ion collision frequencies must share the same energy dependence, differing
only by n_i Z_i^2. This holds for the prescribed common-screening Born model;
it need not hold for charge-dependent quantum cross sections. Ion recoil,
ion heat variables and electron-electron collisions are outside this form.
"""
from dataclasses import dataclass
import numpy as np
from electron_ion_heat import physical_prefactor


@dataclass(frozen=True)
class ElasticMixtureForm:
    collision_matrix: np.ndarray
    scatterer_variance_matrix: np.ndarray
    electron_relative_velocity_map: np.ndarray
    mean_ion_velocity_map: np.ndarray
    total_collision_prefactor: float
    relaxation: str


def elastic_mixture_form(*, velocity_map, ion_density, ion_charges,
                         collision_moments, full_energy_response,
                         relaxation='full_energy'):
    """Return the electron contribution to T times entropy production.

velocity_map maps independent velocity variables to all ion velocities and
the electron velocity (last row). A final variable r_e is appended; its heat
flux is n_e k_B T r_e. The returned matrix acts on [independent velocities,r_e].
Mass/current constraints belong in the caller's velocity map. The matrix
can be added to the ion-only collision form before solving coupled transport.
Inputs I and L are the dimensionless results from electron_ion_heat.py.
"""
    n, z, b, I, L = [np.asarray(v, dtype=float) for v in
                     (ion_density, ion_charges, velocity_map,
                      collision_moments, full_energy_response)]
    if (n.ndim != 1 or not len(n) or z.shape != n.shape or
            b.ndim != 2 or b.shape[0] != len(n)+1 or b.shape[1] < 1 or
            not all(np.all(np.isfinite(v)) for v in [n,z,b,I,L]) or
            np.any(n <= 0) or np.any(z <= 0)):
        raise ValueError('positive ion populations and a finite velocity map required')
    if relaxation not in ['full_energy', 'one_heat_variable']:
        raise ValueError('explicit supported electron relaxation required')
    for a in [I,L]:
        if a.shape != (2,2) or not np.array_equal(a,a.T) or np.any(np.diag(a) <= 0):
            raise ValueError('symmetric positive electron moments required')
        scale = np.sqrt(np.diag(a))
        try:
            np.linalg.cholesky(a/scale[:,None]/scale[None,:])
        except np.linalg.LinAlgError as error:
            raise ValueError('electron response is not positive definite') from error
    weights = n*z*z
    if not np.all(np.isfinite(weights)) or not np.isfinite(np.sum(weights)):
        raise ValueError('ion scattering weights overflow')
    alpha = weights/np.sum(weights)
    mean = alpha@b[:-1]
    relative = b[-1]-mean
    nc = b.shape[1]
    mapping = np.zeros((2,nc+1))
    mapping[0,:nc] = relative
    mapping[1,-1] = 1.
    prefactor = physical_prefactor(float(np.sum(weights)),1.)
    # A sum of positive outer products retains trace-species precision.
    variance = np.zeros((nc+1,nc+1))
    for a, row in zip(alpha,b[:-1],strict=True):
        difference = row-mean
        variance[:nc,:nc] += a*np.outer(difference,difference)
    variance *= prefactor*I[0,0]
    response = I if relaxation == 'one_heat_variable' else np.linalg.inv(L)
    matrix = variance + prefactor*mapping.T@response@mapping
    if not np.all(np.isfinite(matrix)):
        raise ValueError('electron mixture matrix is not finite')
    return ElasticMixtureForm(matrix,variance,relative,mean,prefactor,relaxation)
