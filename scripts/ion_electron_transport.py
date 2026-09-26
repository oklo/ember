"""Classical ion moments coupled to elastic finite-Fermi electron transport.

This is a comparison model: common-screening Born electron-ion scattering,
Maxwellian ions, a leading Brownian ion heat term, and supplied ion collision
moments. Electron-electron collisions, higher-order recoil and correlations
beyond the supplied ion moments are absent. It is not selected by Ember.
"""
from dataclasses import dataclass
import numpy as np
from diffusion_burgers import KB, MU
from diffusion_composition_forces import velocity_basis
from electron_elastic_mixture import elastic_mixture_form


def ion_moment_form(*, velocity_map, mass_numbers, resistance, collision_z,
                    collision_zprime, collision_zdoubleprime):
    """Ion-only quadratic collision cost on [independent velocities, ion r]."""
    b, mass, k, z, zp, zpp = [np.asarray(v, dtype=float) for v in
        (velocity_map, mass_numbers, resistance, collision_z,
         collision_zprime, collision_zdoubleprime)]
    if mass.ndim != 1 or not len(mass) or np.any(mass <= 0) or not np.all(np.isfinite(mass)):
        raise ValueError('positive finite ion masses required')
    count = len(mass)
    if b.ndim != 2 or b.shape[0] != count or not b.shape[1] or not np.all(np.isfinite(b)):
        raise ValueError('finite ion velocity map required')
    for a in [k,z,zp,zpp]:
        if a.shape != (count,count) or not np.all(np.isfinite(a)) or not np.array_equal(a,a.T):
            raise ValueError('finite symmetric ion collision arrays required')
    if np.any(k < 0) or np.any(zp < 0) or np.any(zpp < 0):
        raise ValueError('negative ion collision coefficient')
    nc = b.shape[1]
    m = np.zeros((nc+count,nc+count))
    for i in range(count):
        m[nc+i,nc+i] += .16*k[i,i]*zpp[i,i]
        for j in range(i+1,count):
            mi,mj = mass[i],mass[j]
            total = mi+mj
            drag = k[i,j]
            diff = b[i]-b[j]
            m[:nc,:nc] += drag*np.outer(diff,diff)
            m[:nc,nc+i] -= drag*z[i,j]*mj/total*diff
            m[:nc,nc+j] += drag*z[i,j]*mi/total*diff
            m[nc+i,nc+i] += .4*drag*(3*mi*mi+mj*mj*zp[i,j]+.8*mi*mj*zpp[i,j])/total**2
            m[nc+j,nc+j] += .4*drag*(3*mj*mj+mi*mi*zp[i,j]+.8*mi*mj*zpp[i,j])/total**2
            cross = -.4*drag*mi*mj*(3+zp[i,j]-.8*zpp[i,j])/total**2
            m[nc+i,nc+j] += cross
            m[nc+j,nc+i] += cross
    m[nc:,:nc] = m[:nc,nc:].T
    return m


@dataclass(frozen=True)
class CoupledResponse:
    velocity_map: np.ndarray
    variable_response: np.ndarray
    flux_map: np.ndarray
    mobility: np.ndarray
    collision_matrix: np.ndarray
    ion_collision_matrix: np.ndarray
    electron_collision_matrix: np.ndarray
    brownian_heat_diagonal: np.ndarray
    temperature: float
    energy_scale: float
    scaled_condition_number: float
    scaled_backward_error: float

    def solve(self, reduced_forces):
        force = np.asarray(reduced_forces, dtype=float)
        if force.shape != (len(self.mobility),) or not np.all(np.isfinite(force)):
            raise ValueError('finite composition and scaled thermal forces required')
        variables = self.variable_response@force
        flux = self.flux_map@variables
        nc = self.velocity_map.shape[1]
        return dict(velocity=self.velocity_map@variables[:nc],
                    heat_velocity=variables[nc:], mass_flux=flux[:-1],
                    reduced_heat_flux=float(flux[-1]*self.energy_scale),
                    entropy_from_forces=float(flux@force),
                    entropy_from_collisions=float(variables@self.collision_matrix@variables/self.temperature))


def coupled_elastic_operator(*, density, temperature, mass_fractions, mass_numbers,
                             charges, independent_ions, reference_ion,
                             ion_resistance, ion_z, ion_zprime, ion_zdoubleprime,
                             electron_moments, electron_full_response, energy_scale,
                             electron_relaxation='full_energy'):
    """Prepare reusable responses to matched reduced chemical/thermal forces.

The ion Brownian heat term is (6/5) K_ie r_i^2. It follows from the
Ornstein-Uhlenbeck ion operator with equilibrium electrons at the same T:
the third-degree heat polynomial relaxes three times as fast as drift and
has norm 2/5. Its K_ie is the *unrelaxed* displaced-electron resistance,
not a conductivity or a separately relaxed pair coefficient.
"""
    if not np.all(np.isfinite([density,temperature,energy_scale])) or min(density,temperature,energy_scale) <= 0:
        raise ValueError('positive finite density, temperature and energy scale required')
    x,mass,z = [np.asarray(v,dtype=float) for v in (mass_fractions,mass_numbers,charges)]
    b = velocity_basis(mass_fractions=x,mass_numbers=mass,charges=z,
                       independent_ions=independent_ions,reference_ion=reference_ion)
    ni = density/MU*x/mass
    ne = float(ni@z)
    nc,count = b.shape[1],len(x)
    ions = ion_moment_form(velocity_map=b[:-1],mass_numbers=mass,
                          resistance=ion_resistance,collision_z=ion_z,
                          collision_zprime=ion_zprime,collision_zdoubleprime=ion_zdoubleprime)
    electrons = elastic_mixture_form(velocity_map=b,ion_density=ni,ion_charges=z,
        collision_moments=electron_moments,full_energy_response=electron_full_response,
        relaxation=electron_relaxation)
    ion_part = np.zeros((nc+count+1,nc+count+1))
    ion_part[:-1,:-1] = ions
    electron_part = np.zeros_like(ion_part)
    index = np.r_[np.arange(nc), nc+count]
    electron_part[np.ix_(index,index)] = electrons.collision_matrix
    weights = ni*z*z
    brownian = 1.2*electrons.total_collision_prefactor*np.asarray(electron_moments)[0,0]*weights/sum(weights)
    electron_part[nc+np.arange(count),nc+np.arange(count)] += brownian
    m = ion_part+electron_part
    if not np.all(np.isfinite(m)) or np.any(np.diag(m) <= 0):
        raise ValueError('undetermined or nonfinite ion/electron transport')
    scale = np.sqrt(np.diag(m))
    scaled = m/scale[:,None]/scale[None,:]
    condition = float(np.linalg.cond(scaled))
    if not np.isfinite(condition) or condition > 1e12:
        raise ValueError('ill-conditioned ion/electron transport')
    try:
        chol = np.linalg.cholesky(scaled)
    except np.linalg.LinAlgError as error:
        raise ValueError('collision model is not positive') from error
    g = np.zeros((nc+1,len(m)))
    g[np.arange(nc),np.arange(nc)] = density*x[np.asarray(independent_ions)]
    g[-1,nc:] = np.r_[ni,ne]*KB*temperature/energy_scale
    target = (g/scale).T
    response = np.linalg.solve(chol.T,np.linalg.solve(chol,target))
    residual = float(np.linalg.norm(scaled@response-target,np.inf)/(
        np.linalg.norm(scaled,np.inf)*np.linalg.norm(response,np.inf)+np.linalg.norm(target,np.inf)))
    if not np.all(np.isfinite(response)) or residual > 1e-12:
        raise ValueError('ion/electron response residual failed')
    variables = temperature*response/scale[:,None]
    return CoupledResponse(b,variables,g,g@variables,m,ion_part,electron_part,
                           brownian,temperature,energy_scale,condition,residual)


def heat_decomposition(mobility, *, temperature, energy_scale):
    """Return h_transport and conductivity at zero represented element flux.

Q' = h_transport dot j - kappa grad T. Full material energy also contains
the EOS exchange enthalpies dot j. This identity does not justify substituting
a different conductivity while retaining unmatched cross coefficients.
"""
    l = np.asarray(mobility,dtype=float)
    if (l.ndim != 2 or l.shape[0] != l.shape[1] or len(l) < 2 or
            not np.all(np.isfinite(l)) or np.any(np.diag(l) <= 0) or
            not np.all(np.isfinite([temperature,energy_scale])) or min(temperature,energy_scale) <= 0):
        raise ValueError('positive finite response scales and mobility required')
    s = np.sqrt(np.diag(l))
    a = l/s[:,None]/s[None,:]
    if np.max(abs(a-a.T)) > 1e-12:
        raise ValueError('reciprocal mobility required')
    try:
        chol = np.linalg.cholesky((a+a.T)/2)
    except np.linalg.LinAlgError as error:
        raise ValueError('positive mobility required') from error
    # Cholesky's last pivot is the normalized Schur complement, avoiding
    # cancellation from subtracting almost equal heat responses.
    coefficient = np.linalg.solve(a[:-1,:-1],a[-1,:-1])
    h_transport = energy_scale*s[-1]*coefficient/s[:-1]
    kappa = (energy_scale/temperature*s[-1]*chol[-1,-1])**2
    return h_transport,float(kappa)
