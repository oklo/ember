"""Reduce joint electron-ion and electron-electron energy modes for a mixture.

The ion kernel has prescribed common Born screening. Its unrelaxed
scatterer-velocity variance must remain separate from energy relaxation.
Electron-electron scattering has an exact momentum null mode. This module
does not supply ionization, correlations, recoil, or stellar evolution.
"""
from dataclasses import dataclass
import numpy as np


def symmetric(matrix):
    a = np.asarray(matrix, dtype=float)
    if a.ndim != 2 or a.shape[0] != a.shape[1] or not np.isfinite(a).all():
        raise ValueError('finite square collision matrix required')
    scale = max(float(np.max(abs(a))), np.finfo(float).tiny)
    if np.max(abs(a-a.T)) > 2e-12*scale:
        raise ValueError('collision matrix is not symmetric')
    return (a+a.T)/2


def positive_solve(matrix, rhs):
    a = symmetric(matrix)
    if np.any(np.diag(a) <= 0):
        raise ValueError('positive collision diagonal required')
    scale = np.sqrt(np.diag(a))
    normalized = a/scale[:,None]/scale[None,:]
    chol = np.linalg.cholesky(normalized)
    rhs = np.asarray(rhs,dtype=float)
    return np.linalg.solve(chol.T,np.linalg.solve(chol,rhs/scale[:,None]))/scale[:,None]


@dataclass(frozen=True)
class ReducedElectronCollision:
    collision_matrix: np.ndarray
    retained_collision: np.ndarray
    retained_response: np.ndarray
    higher_mode_map: np.ndarray
    electron_electron_prefactor_ratio: float
    scatterer_fractions: np.ndarray


def reduce_electron_collision(*, electron_ion_matrix, electron_electron_matrix,
                              ion_density, ion_charges, electron_electron_scale=1.):
    ei, ee = symmetric(electron_ion_matrix), symmetric(electron_electron_matrix)
    n, z = np.asarray(ion_density,dtype=float), np.asarray(ion_charges,dtype=float)
    if (ei.shape != ee.shape or len(ei) < 2 or n.ndim != 1 or z.shape != n.shape or
            np.any(n <= 0) or np.any(z <= 0) or not np.isfinite(n).all() or
            not np.isfinite(z).all() or not np.isfinite(electron_electron_scale) or
            electron_electron_scale < 0):
        raise ValueError('positive ion populations and matching energy modes required')
    # A momentum penalty here would spuriously oppose common electron drift.
    if np.any(ee[0] != 0) or np.any(ee[:,0] != 0):
        raise ValueError('electron-electron momentum null mode must be exact')
    weights = n*z*z
    total = float(np.sum(weights))
    ne = float(n@z)
    if not np.isfinite(total) or not np.isfinite(ne) or min(total,ne) <= 0:
        raise ValueError('invalid scattering population')
    ratio = electron_electron_scale*ne/total
    joint = ei+ratio*ee
    if len(joint) > 2:
        higher_map = -positive_solve(joint[2:,2:],joint[2:,:2])
        retained = symmetric(joint[:2,:2]+joint[:2,2:]@higher_map)
    else:
        higher_map = np.empty((0,2))
        retained = joint.copy()
    response = symmetric(positive_solve(retained,np.eye(2)))
    return ReducedElectronCollision(joint,retained,response,higher_map,ratio,weights/total)
