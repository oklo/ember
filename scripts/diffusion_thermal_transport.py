"""Material exchange enthalpies and matched diffusion/heat conventions.

This transforms specified fluxes, forces and kinetic coefficients together;
it does not supply physical kinetic coefficients. The neutral material EOS
already includes electrons. Radiation is a separate transport channel.
Independent mass exchanges replace a common reference species after imposing
zero net baryonic mass flux and electric current. See docs/DIFFUSION_THERMAL.md.
"""
from dataclasses import dataclass
import numpy as np


def exchange_enthalpy(*, temperature, gradient_phi_lnT,
                      gradient_phi_lnrho, material_delta):
    """Partial enthalpy differences in erg/g for represented mass exchanges.

    Phi is material specific Helmholtz free energy divided by T, in erg/g/K.
    Inputs are the log-T and log-rho derivatives of its composition gradient.
    material_delta = P_material,lnT / P_material,lnrho at fixed composition.
    Supply only finite, selected composition directions, with no abundance floor.
    """
    gt, gr = [np.asarray(v, dtype=float) for v in
              (gradient_phi_lnT, gradient_phi_lnrho)]
    if (gt.ndim != 1 or not len(gt) or gr.shape != gt.shape or
            not np.all(np.isfinite(gt)) or not np.all(np.isfinite(gr)) or
            not np.isfinite(temperature) or temperature <= 0 or
            not np.isfinite(material_delta)):
        raise ValueError('finite material derivatives and positive temperature required')
    result = temperature * (material_delta * gr - gt)
    if not np.all(np.isfinite(result)):
        raise ValueError('nonfinite exchange enthalpy')
    return result


@dataclass(frozen=True)
class HeatFluxConvention:
    """Transform full material energy flux into reduced heat flux and back.

    The final flux coordinate is Q/e0, where e0 is a positive specific-energy
    scale in erg/g; the other coordinates are independent mass fluxes j_i.
    The final force coordinate is e0*grad(1/T); the others are -grad(q_i/T).
    All flux coordinates therefore have the same units, as do all forces.
    Q is the material energy flux relative to barycentric motion. Its reduced
    value is Q' = Q - sum(h_i*j_i). The conjugate chemical forces change by
    +h_i*grad(1/T). The normalization avoids huge disparate matrix units.

    No collision model, diffusion velocity or stellar heat flux is inferred.
    A supplied transport matrix need not be symmetric or positive for the
    algebraic transformation; those are separate physical requirements.
    """
    exchange_enthalpy_erg_g: np.ndarray
    energy_scale_erg_g: float

    def __post_init__(self):
        h = np.array(self.exchange_enthalpy_erg_g, dtype=float, copy=True)
        if (h.ndim != 1 or not len(h) or not np.all(np.isfinite(h)) or
                not np.isfinite(self.energy_scale_erg_g) or self.energy_scale_erg_g <= 0):
            raise ValueError('finite exchange enthalpies and positive energy scale required')
        h.setflags(write=False)
        object.__setattr__(self, 'exchange_enthalpy_erg_g', h)

    @property
    def scaled_enthalpy(self):
        return self.exchange_enthalpy_erg_g / self.energy_scale_erg_g

    def _vector(self, value):
        a = np.array(value, dtype=float, copy=True)
        if a.shape != (len(self.exchange_enthalpy_erg_g) + 1,) or not np.all(np.isfinite(a)):
            raise ValueError('finite species-plus-energy vector required')
        return a

    def reduced_forces(self, full_forces):
        f = self._vector(full_forces)
        f[:-1] += self.scaled_enthalpy * f[-1]
        return f

    def full_forces(self, reduced_forces):
        f = self._vector(reduced_forces)
        f[:-1] -= self.scaled_enthalpy * f[-1]
        return f

    def reduced_fluxes(self, full_fluxes):
        j = self._vector(full_fluxes)
        j[-1] -= self.scaled_enthalpy @ j[:-1]
        return j

    def full_fluxes(self, reduced_fluxes):
        j = self._vector(reduced_fluxes)
        j[-1] += self.scaled_enthalpy @ j[:-1]
        return j

    def _matrix(self, matrix, inverse):
        n = len(self.exchange_enthalpy_erg_g) + 1
        l = np.asarray(matrix, dtype=float)
        if l.shape != (n, n) or not np.all(np.isfinite(l)):
            raise ValueError('finite species-plus-energy transport matrix required')
        a = np.eye(n)
        a[-1, :-1] = (1 if inverse else -1) * self.scaled_enthalpy
        return a @ l @ a.T

    def reduced_matrix(self, full_matrix):
        return self._matrix(full_matrix, False)

    def full_matrix(self, reduced_matrix):
        return self._matrix(reduced_matrix, True)
