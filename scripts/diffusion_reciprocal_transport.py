"""Reciprocal form of the Maxwellian Burgers momentum and heat equations.

This is an explicitly classical kinetic model, not a degenerate-electron
prescription. Pair coefficients are supplied. Velocities satisfy mass/current
constraints before solving, and the resulting species/heat transport matrix
is conjugate to the reduced forces in diffusion_thermal_transport.py.
Equations follow Paxton et al. 2018 (1), (2), with baryonic ion mass conserved.
"""
from dataclasses import dataclass
import numpy as np
from diffusion_burgers import KB, MU, ME
from diffusion_composition_forces import velocity_basis


@dataclass(frozen=True)
class ClassicalFlux:
    velocity_cm_s: np.ndarray
    residual_heat_velocity_cm_s: np.ndarray
    independent_mass_flux_g_cm2_s: np.ndarray
    reduced_heat_flux_erg_cm2_s: float
    entropy_production_erg_cm3_s_K: float
    collision_entropy_erg_cm3_s_K: float


@dataclass(frozen=True)
class ClassicalOperator:
    velocity_basis: np.ndarray
    variable_response: np.ndarray
    flux_map: np.ndarray
    mobility: np.ndarray
    collision_matrix: np.ndarray
    temperature_K: float
    energy_scale_erg_g: float
    scaled_condition_number: float
    scaled_backward_error: float

    def solve(self, reduced_forces):
        f=np.asarray(reduced_forces,dtype=float)
        if f.shape!=(self.mobility.shape[0],) or not np.all(np.isfinite(f)):
            raise ValueError('finite independent-species and scaled thermal forces required')
        v=self.variable_response@f
        nc=self.velocity_basis.shape[1];j=self.flux_map@v
        return ClassicalFlux(self.velocity_basis@v[:nc],v[nc:].copy(),j[:-1].copy(),
            float(j[-1]*self.energy_scale_erg_g),float(j@f),
            float(v@self.collision_matrix@v/self.temperature_K))


def maxwellian_operator(*, density, temperature, mass_fractions, mass_numbers,
                        charges, independent_ions, reference_ion, resistance,
                        collision_z, collision_zprime, collision_zdoubleprime,
                        energy_scale):
    """Construct a reusable classical response to reduced chemical/heat forces.

    Ion populations must be positive, with fixed integer charges. Unselected
    ions are stationary but retain their collisional heat variables. Electrons
    are added as the final population. All supplied collision arrays include
    that electron row and self-collisions, using the Burgers convention.
    The final force is energy_scale*grad(1/T); chemical forces are in erg/g/K/cm.
    This function does not infer composition forces, gravity or degeneracy.
    """
    b=velocity_basis(mass_fractions=mass_fractions,mass_numbers=mass_numbers,
                     charges=charges,independent_ions=independent_ions,reference_ion=reference_ion)
    if (not np.all(np.isfinite([density,temperature,energy_scale])) or
            min(density,temperature,energy_scale)<=0):
        raise ValueError('positive density, temperature and energy scale required')
    x,ai,zi=[np.asarray(v,dtype=float) for v in (mass_fractions,mass_numbers,charges)]
    mass=np.append(ai,ME/MU);nion=density/MU*x/ai;n=np.append(nion,nion@zi)
    count=len(n);nc=b.shape[1]
    k,z,zp,zpp=[np.asarray(v,dtype=float) for v in
               (resistance,collision_z,collision_zprime,collision_zdoubleprime)]
    for v in [k,z,zp,zpp]:
        if v.shape!=(count,count) or not np.all(np.isfinite(v)) or not np.array_equal(v,v.T):
            raise ValueError('finite symmetric population collision arrays required')
    if np.any(k<0) or np.any(zp<0) or np.any(zpp<0):
        raise ValueError('negative resistance or heat coefficient')
    # Divide each Burgers heat equation by 5/2. This makes its velocity
    # coupling the transpose of the momentum equation's heat coupling.
    m=np.zeros((nc+count,nc+count))
    for i in range(count):
        m[nc+i,nc+i]+=.16*k[i,i]*zpp[i,i]
        for j in range(i+1,count):
            mi,mj=mass[i],mass[j];total=mi+mj;drag=k[i,j]
            diff=b[i]-b[j]
            m[:nc,:nc]+=drag*np.outer(diff,diff)
            m[:nc,nc+i]-=drag*z[i,j]*mj/total*diff
            m[:nc,nc+j]+=drag*z[i,j]*mi/total*diff
            m[nc+i,nc+i]+=.4*drag*(3*mi*mi+mj*mj*zp[i,j]+.8*mi*mj*zpp[i,j])/total**2
            m[nc+j,nc+j]+=.4*drag*(3*mj*mj+mi*mi*zp[i,j]+.8*mi*mj*zpp[i,j])/total**2
            cross=-.4*drag*mi*mj*(3+zp[i,j]-.8*zpp[i,j])/total**2
            m[nc+i,nc+j]+=cross;m[nc+j,nc+i]+=cross
    m[nc:,:nc]=m[:nc,nc:].T
    diagonal=np.diag(m)
    if np.any(diagonal<=0):raise ValueError('undetermined classical transport variables')
    scale=np.sqrt(diagonal);scaled=m/scale[:,None]/scale[None,:]
    condition=float(np.linalg.cond(scaled))
    if not np.isfinite(condition) or condition>1e12:
        raise ValueError('ill-conditioned classical transport')
    try:cholesky=np.linalg.cholesky(scaled)
    except np.linalg.LinAlgError as error:
        raise ValueError('collision coefficients do not give positive entropy production') from error
    g=np.zeros((nc+1,nc+count))
    g[np.arange(nc),np.arange(nc)]=density*x[np.asarray(independent_ions)]
    # Burgers r_s gives residual heat flux n_s kT r_s. Scaling energy leaves
    # every flux coordinate in g/cm^2/s and avoids disparate numerical units.
    g[-1,nc:]=n*KB*temperature/energy_scale
    target=(g/scale).T
    response=np.linalg.solve(cholesky.T,np.linalg.solve(cholesky,target))
    backward=float(np.linalg.norm(scaled@response-target,np.inf)/(
        np.linalg.norm(scaled,np.inf)*np.linalg.norm(response,np.inf)+np.linalg.norm(target,np.inf)))
    if not np.all(np.isfinite(response)) or backward>1e-12:
        raise ValueError('classical response residual failed')
    variables=temperature*response/scale[:,None]
    mobility=g@variables
    return ClassicalOperator(b,variables,g,mobility,m,temperature,energy_scale,condition,backward)
