"""Offline composition forces and a constrained momentum-only transport solve.

The chemical force follows Beznogov, Potekhin & Yakovlev (2016), equations
2 and 10, with the temperature derivative taken at fixed material pressure
and composition: https://doi.org/10.1093/mnras/stw751 . Its thermal-force
convention must be matched to the thermal transport coefficients. Away from
isothermality it is not automatically interchangeable with the ideal-ion
pressure gradients in Burgers equations.

Neutral-matter free-energy derivatives include electrons already. The velocity
basis enforces baryonic mass conservation and zero electric current before
solving, so no extra electron chemical force or electric field is added.
Fixed ions have zero velocity explicitly. This supports only composition
directions provided by the EOS; it does not supply missing metal derivatives.
No stellar integration, collision prescription or heat-flow closure is here.
"""
from dataclasses import dataclass
import numpy as np


def chemical_acceleration(*, temperature, hessian_phi, gradient_phi_lnrho,
                          material_delta, composition_gradient,
                          logarithmic_density_gradient,
                          logarithmic_temperature_gradient):
    """Return the outward gradient force per unit exchanged baryonic mass.

    phi = specific material Helmholtz free energy / T, in erg/g/K.
    Independent mass fractions replace a common reference ion at fixed metals.
    Spatial gradients use outward radius in cm. The result is in cm/s^2,
    positive for an increasing chemical potential (driving flux inward).
    material_delta = (dP/dlnT)/(dP/dlnrho) at fixed composition, without
    radiation. Temperature and composition derivatives at fixed P remove
    the explicit temperature part, leaving T [H dc/dr + phi_c,lnrho
    (dlnrho/dr + material_delta dlnT/dr)]. Thermal diffusion remains separate.
    """
    h,g,c=[np.asarray(x,dtype=float) for x in
           (hessian_phi,gradient_phi_lnrho,composition_gradient)]
    if g.ndim!=1 or not len(g) or h.shape!=(len(g),len(g)) or c.shape!=g.shape:
        raise ValueError('one Hessian and gradient for the active composition directions required')
    if not all(np.all(np.isfinite(x)) for x in (h,g,c)):
        raise ValueError('nonfinite composition derivative')
    if not np.array_equal(h,h.T):raise ValueError('composition Hessian must be symmetric')
    values=[temperature,material_delta,logarithmic_density_gradient,logarithmic_temperature_gradient]
    if not np.all(np.isfinite(values)) or temperature<=0:
        raise ValueError('finite physical state and gradients required')
    return temperature*(h@c+g*(logarithmic_density_gradient+
                               material_delta*logarithmic_temperature_gradient))


def velocity_basis(*, mass_fractions, mass_numbers, charges,
                   independent_ions, reference_ion):
    """Columns are velocities of independent ions; electrons occupy last row.

    Every supplied ion must have positive abundance; absent ions must be
    removed by the caller. Other ions are stationary. Charges are fixed
    positive integers; changing ionization stages are not represented.
    """
    x,a,z=[np.asarray(v,dtype=float) for v in (mass_fractions,mass_numbers,charges)]
    if x.ndim!=1 or len(x)<2 or any(v.shape!=x.shape for v in (a,z)):
        raise ValueError('two or more ion populations required')
    if not all(np.all(np.isfinite(v)) for v in (x,a,z)) or np.any(x<=0):
        raise ValueError('positive active mass fractions required; no abundance floor')
    if abs(sum(x)-1)>1e-12 or np.any(a<1) or np.any(a!=np.floor(a)):
        raise ValueError('normalized baryonic fractions and integer masses required')
    if np.any(z<1) or np.any(z!=np.floor(z)) or np.any(z>a):
        raise ValueError('fixed positive integer ion charges required')
    indices=np.asarray(independent_ions)
    if (indices.ndim!=1 or not len(indices) or indices.dtype.kind not in 'iu' or
            np.any(indices<0) or np.any(indices>=len(x)) or len(set(indices))!=len(indices) or
            not isinstance(reference_ion,(int,np.integer)) or not 0<=reference_ion<len(x) or
            reference_ion in indices):
        raise ValueError('distinct valid independent and reference ions required')
    ye=np.dot(x,z/a);basis=np.zeros((len(x)+1,len(indices)))
    for k,i in enumerate(indices):
        basis[i,k]=1.
        basis[reference_ion,k]=-x[i]/x[reference_ion]
        basis[-1,k]=x[i]*(z[i]/a[i]-z[reference_ion]/a[reference_ion])/ye
    return basis


@dataclass(frozen=True)
class ExchangeSolution:
    velocity_cm_s: np.ndarray
    independent_mass_flux_g_cm2_s: np.ndarray
    frictional_heating_erg_cm3_s: float
    driving_power_erg_cm3_s: float
    scaled_condition_number: float
    scaled_backward_error: float


def solve_momentum(*, density, mass_fractions, mass_numbers, charges,
                   independent_ions, reference_ion, resistance,
                   chemical_driving_acceleration):
    """Solve only the pair-drag balance in the conserved velocity subspace.

    K_ij is supplied in g/(cm^3 s), including electrons in the final row.
    The caller explicitly supplies chemical driving accelerations in the
    same independent-ion order. Positive chemical force gives negative flux.
    No residual heat flow or thermal diffusion term is inferred. A zero
    electron resistance is allowed when the constrained system is determined.
    """
    b=velocity_basis(mass_fractions=mass_fractions,mass_numbers=mass_numbers,
                     charges=charges,independent_ions=independent_ions,reference_ion=reference_ion)
    x=np.asarray(mass_fractions,dtype=float);k=np.asarray(resistance,dtype=float)
    force=np.asarray(chemical_driving_acceleration,dtype=float)
    if (not np.isfinite(density) or density<=0 or k.shape!=(len(x)+1,len(x)+1) or
            not np.all(np.isfinite(k)) or np.any(k<0) or not np.array_equal(k,k.T) or
            force.shape!=(b.shape[1],) or not np.all(np.isfinite(force))):
        raise ValueError('positive density, symmetric resistance and finite chemical forces required')
    i,j=np.triu_indices(len(x)+1,1);differences=b[i]-b[j]
    # Sum pair differences directly to avoid subtracting diagonal row sums.
    matrix=np.einsum('p,pi,pj->ij',k[i,j],differences,differences)
    diagonal=np.diag(matrix)
    if np.any(diagonal<=0):raise ValueError('undetermined constrained motion')
    scale=np.sqrt(diagonal);scaled=matrix/scale[:,None]/scale[None,:]
    condition=float(np.linalg.cond(scaled))
    if not np.isfinite(condition) or condition>1e12:
        raise ValueError('ill-conditioned constrained momentum equations')
    try:np.linalg.cholesky(scaled)
    except np.linalg.LinAlgError as error:raise ValueError('nonpositive constrained resistance') from error
    rhs=-density*x[np.asarray(independent_ions)]*force
    target=rhs/scale;v=np.linalg.solve(scaled,target)
    backward=float(np.linalg.norm(scaled@v-target,np.inf)/max(
        np.linalg.norm(scaled,np.inf)*np.linalg.norm(v,np.inf)+np.linalg.norm(target,np.inf),
        np.finfo(float).tiny))
    if not np.all(np.isfinite(v)) or backward>1e-12:
        raise ValueError('constrained momentum residual failed')
    coefficients=v/scale;velocity=b@coefficients
    flux=density*x[np.asarray(independent_ions)]*coefficients
    heat=float(np.dot(k[i,j],(velocity[i]-velocity[j])**2))
    power=float(-flux@force)
    return ExchangeSolution(velocity,flux,heat,power,condition,backward)
