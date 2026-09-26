"""Offline Burgers equation solver; no collision prescription or evolution hook.

The caller supplies collision coefficients and ion pressure gradients in cgs.
Gravity is an input, never inferred from an ideal-gas total pressure. Electrons
are last; their momentum equation is replaced by charge and baryon-flux
constraints. Classical residual heat-flow equations can be solved explicitly,
or omitted as an explicitly selected approximation. No degeneracy switch or
thermal-diffusion prescription is inferred here.

Equations: Burgers momentum/heat balance, as written in Paxton et al. (2018),
equations (1), (2), (9), Appendix C; https://arxiv.org/abs/1710.08424.
Written from these equations, independently of MESA's implementation.
"""
from dataclasses import dataclass

import numpy as np

KB = 1.380649e-16
MU = 1.66053906660e-24
ME = 9.1093837015e-28


@dataclass(frozen=True)
class Solution:
    velocity_cm_s: np.ndarray
    residual_heat_velocity_cm_s: np.ndarray
    electric_force_dyn: float  # eE, outward positive
    scaled_backward_error: float
    scaled_condition_number: float


def solve(*, number_density, mass_u, charge, ion_pressure_gradient,
          gravity, temperature, temperature_gradient, resistance, z,
          zprime, zdoubleprime, heat_flow, ion_extra_force=None):
    """Solve local velocities in the frame of zero diffusive baryon flux.

    number_density: N positive populations in cm^-3, electrons last.
    mass_u: collision masses in m_u; integer baryonic A for ions, m_e/m_u
      for electrons. The mass constraint counts ions only, consistent with
      Ember's conserved baryonic mass coordinate. Electron mass remains in
      collision kinematics. Nuclear mass defects are not included here.
    ion_pressure_gradient: N-1 gradients dP_i/dr, dyn/cm^3, outward r.
    gravity: positive inward acceleration, cm/s^2.
    temperature_gradient: dT/dr, K/cm.
    resistance: symmetric K_ij in g/(cm^3 s), including self-collisions.
    z, zprime, zdoubleprime: symmetric, dimensionless collision integrals.
    ion_extra_force: optional outward force per ion, dyn. Nonideal forces
      require an independently validated thermodynamic derivation.
    heat_flow: 'classical' or 'suppressed'; no automatic physical selection.

    The classical heat equations assume Maxwellian populations. They are not
    accepted for degenerate electrons. Suppression is a comparison option,
    not a claim that thermal diffusion vanishes in the saved stellar model.
    Zero-abundance species must be removed explicitly by the caller. This
    function performs no abundance flooring, clipping or renormalization.
    """
    n, a, q = [np.asarray(v, dtype=float) for v in
               (number_density, mass_u, charge)]
    count = len(n)
    if count < 2 or any(v.shape != (count,) for v in (n, a, q)):
        raise ValueError('one or more ions plus electrons required')
    if not all(np.all(np.isfinite(v)) for v in (n, a, q)):
        raise ValueError('nonfinite population, mass or charge')
    if np.any(n <= 0) or np.any(a <= 0) or np.any(q[:-1] <= 0):
        raise ValueError('positive populations, masses and ion charges required')
    if q[-1] != -1 or abs(a[-1]/(ME/MU)-1) > 1e-10:
        raise ValueError('last population must be electrons')
    if np.any(a[:-1] != np.floor(a[:-1])):
        raise ValueError('integer baryonic ion masses required')
    if abs(n @ q) > 1e-12*(n @ np.abs(q)):
        raise ValueError('charge-neutral input required')
    if not (np.isfinite(gravity) and gravity >= 0 and
            np.isfinite(temperature) and temperature > 0 and
            np.isfinite(temperature_gradient)):
        raise ValueError('invalid gravity or temperature')
    if heat_flow not in ('classical', 'suppressed'):
        raise ValueError('explicit heat-flow approximation required')
    dp = np.asarray(ion_pressure_gradient, dtype=float)
    extra = np.zeros(count-1) if ion_extra_force is None else np.asarray(ion_extra_force, dtype=float)
    if any(v.shape != (count-1,) or not np.all(np.isfinite(v)) for v in (dp, extra)):
        raise ValueError('one finite pressure gradient and extra force per ion required')
    k, zz, zp, zpp = [np.asarray(v, dtype=float) for v in
                      (resistance, z, zprime, zdoubleprime)]
    for array in (k, zz, zp, zpp):
        if array.shape != (count, count) or not np.all(np.isfinite(array)):
            raise ValueError('finite square collision matrices required')
        if not np.array_equal(array, array.T):
            raise ValueError('collision matrices must be symmetric')
    if np.any(k < 0) or np.any(zp < 0) or np.any(zpp < 0):
        raise ValueError('negative resistance or heat collision coefficient')
    # A disconnected collision graph has undetermined relative velocities.
    reached = {0}
    while True:
        expanded = reached | {j for i in reached for j in range(count) if i != j and k[i,j] > 0}
        if expanded == reached:
            break
        reached = expanded
    if len(reached) != count:
        raise ValueError('disconnected collision graph')

    full = heat_flow == 'classical'
    size = (2*count if full else count)+1
    matrix, rhs = np.zeros((size, size)), np.zeros(size)
    # Ion momentum. Divide each equation by its own n_i, preserving trace
    # limits without imposing an abundance floor.
    for i in range(count-1):
        rhs[i] = dp[i]/n[i] + a[i]*MU*gravity - extra[i]
        matrix[i,-1] = q[i]
        for j in range(count):
            if i == j:
                continue
            drag = k[i,j]/n[i]
            matrix[i,i] -= drag
            matrix[i,j] += drag
            if full:
                matrix[i,count+i] += drag*zz[i,j]*a[j]/(a[i]+a[j])
                matrix[i,count+j] -= drag*zz[i,j]*a[i]/(a[i]+a[j])
    if full:
        # One heat equation for every population, including electrons.
        for i in range(count):
            row = count-1+i
            rhs[row] = 2.5*KB*temperature_gradient
            matrix[row,count+i] -= .4*k[i,i]*zpp[i,i]/n[i]
            for j in range(count):
                if i == j:
                    continue
                drag = k[i,j]/n[i]
                total = a[i]+a[j]
                matrix[row,i] += 2.5*drag*zz[i,j]*a[j]/total
                matrix[row,j] -= 2.5*drag*zz[i,j]*a[j]/total
                matrix[row,count+i] -= drag*((3*a[i]**2+a[j]**2*zp[i,j])
                                            +.8*a[i]*a[j]*zpp[i,j])/total**2
                matrix[row,count+j] += drag*a[i]*a[j]*(3+zp[i,j]-.8*zpp[i,j])/total**2
    # Baryonic mass and electric current, with no correction after solving.
    baryons = n[:-1]*a[:-1]
    matrix[-2,:count-1] = baryons/np.sum(baryons)
    matrix[-1,:count] = n*q/np.sum(n*np.abs(q))

    # Numerical scaling only, using physical force and drag units; no total
    # pressure or approximate hydrostatic relation enters the rescaling.
    force = max(np.max(np.abs(rhs[:-2])), np.max(a[:-1]*MU*gravity))
    if force == 0:
        force = KB*temperature  # reference force over an arbitrary 1 cm
    rates = k/n[:,None]
    np.fill_diagonal(rates, 0)
    speed = force/np.max(rates)
    columns = np.full(size, speed)
    columns[-1] = force
    scaled = matrix*columns[None,:]
    rows = np.max(np.abs(scaled), axis=1)
    if np.any(rows == 0) or not np.all(np.isfinite(rows)):
        raise ValueError('singular or overflowing equation scale')
    scaled /= rows[:,None]
    target = rhs/rows
    condition = float(np.linalg.cond(scaled))
    if not np.isfinite(condition) or condition > 1e13:
        raise ValueError('ill-conditioned diffusion equations')
    dimensionless = np.linalg.solve(scaled, target)
    backward = float(np.linalg.norm(scaled@dimensionless-target, ord=np.inf) /
                     max(np.linalg.norm(scaled, ord=np.inf)*np.linalg.norm(dimensionless, ord=np.inf)
                         + np.linalg.norm(target, ord=np.inf), np.finfo(float).tiny))
    if not np.all(np.isfinite(dimensionless)) or backward > 1e-12:
        raise ValueError('diffusion linear solve failed its residual check')
    value = columns*dimensionless
    return Solution(value[:count].copy(), value[count:2*count].copy() if full else np.zeros(count),
                    float(value[-1]), backward, condition)
