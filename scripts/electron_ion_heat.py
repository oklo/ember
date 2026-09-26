"""Finite-degeneracy electron drift/heat moments for elastic fixed-ion scattering.

The two-mode Galerkin collision matrix and exact elastic Lorentz response
use the same prescribed collision kernel. They do not include electron-electron
collisions, ion recoil, ion heat variables, relativity or correlations. The
Born wrapper retains the existing screened potential; it selects no stellar
screening prescription. See docs/DIFFUSION_ELECTRON_HEAT.md for the derivation.
"""
from dataclasses import dataclass
import math
import numpy as np
from scipy.integrate import quad
from scipy.special import expit
from electron_ion_born import coulomb_bracket, KB, ME, HBAR, E2


def _integral(eta, function, *, tail, epsrel, epsabs=1e-100):
    low=max(-eta,-tail);high=max(-eta,0)+tail
    points=[v for v in [-20.,-5.,0.,5.,20.] if low<v<high]
    def integrand(y):
        x=max(0.,eta+y)
        return function(x,y)*expit(y)*expit(-y)
    return quad(integrand,low,high,points=points,epsrel=epsrel,epsabs=epsabs,limit=200)[0]


def _parameters(eta,tail,epsrel):
    if (not np.isfinite(eta) or not np.isfinite(tail) or tail<40 or tail>100 or
            not np.isfinite(epsrel) or not 1e-14<=epsrel<=1e-5):
        raise ValueError('finite chemical potential, controlled tail and quadrature required')


@dataclass(frozen=True)
class CurrentBasis:
    eta: float
    normalization: float
    centered_mean: float
    variance: float
    tail: float
    epsrel: float

    def heat_polynomial(self,y):
        return (y-self.centered_mean)/self.variance


def current_basis(eta, *, tail=50., epsrel=2e-11):
    """Normalize number current and a zero-current electron heat mode.

    Weight is x^(3/2) f(1-f). Mean and variance use x-eta to retain the
    narrow thermal response around a much larger degenerate electron energy.
    The mode makes the residual heat flux ne*kT*r_e at any degeneracy.
    """
    _parameters(eta,tail,epsrel)
    scale=max(eta,1.)
    w=lambda x,y:(x/scale)**1.5
    args=dict(tail=tail,epsrel=epsrel)
    norm=_integral(eta,w,**args)
    if not np.isfinite(norm) or norm<=0:raise ValueError('electron current normalization underflows')
    mean=_integral(eta,lambda x,y:w(x,y)*y,epsabs=norm*epsrel*.05,**args)/norm
    variance=_integral(eta,lambda x,y:w(x,y)*(y-mean)**2,epsabs=norm*epsrel*.05,**args)/norm
    normalization=norm*scale**1.5
    if not np.all(np.isfinite([mean,variance,normalization])) or variance<=0:
        raise ValueError('nonpositive or nonfinite electron heat variance')
    return CurrentBasis(float(eta),normalization,float(mean),float(variance),tail,epsrel)


def collision_response(eta,kernel, *, basis=None,tail=50.,epsrel=2e-11):
    """Two-mode resistance and exact fixed-ion response for a stated kernel.

    kernel(x) is x^(3/2) times collision frequency, with an arbitrary common
    positive frequency normalization. Born screening uses kernel=L(Bx).
    Returned two-mode collision moments integrate kernel*f(1-f)*p_a*p_b.
    The exact mobility integrates x^3/kernel*f(1-f)*(1,s) outer (1,s),
    divided by the squared current normalization, where s=x-eta-<x-eta>.
    Multiply resistance by the physical common prefactor, and divide mobility
    by it. This compares the same collision model at the same degeneracy.
    """
    _parameters(eta,tail,epsrel)
    if basis is None:basis=current_basis(eta,tail=tail,epsrel=epsrel)
    if basis.eta!=eta or basis.tail!=tail or basis.epsrel!=epsrel:
        raise ValueError('current basis must match the state and quadrature')
    def positive(x):
        v=float(kernel(x))
        if not np.isfinite(v) or v<0 or (x>0 and v==0):
            raise ValueError('positive finite scattering kernel required at positive energy')
        return v
    args=dict(tail=tail,epsrel=epsrel)
    m00=_integral(eta,lambda x,y:positive(x),**args)
    if not np.isfinite(m00) or m00<=0:raise ValueError('collision normalization underflows')
    m01=_integral(eta,lambda x,y:positive(x)*basis.heat_polynomial(y),epsabs=m00*epsrel*.05,**args)
    m11=_integral(eta,lambda x,y:positive(x)*basis.heat_polynomial(y)**2,epsabs=m00*epsrel*.05,**args)
    collision=np.array([[m00,m01],[m01,m11]])
    # Both currents are independent conserved moments, so the exact elastic
    # relaxation has a matrix of inverse-frequency rather than frequency moments.
    def inverse(x):return x**3/positive(x) if x>0 else 0.
    v00=_integral(eta,lambda x,y:inverse(x),**args)
    v01=_integral(eta,lambda x,y:inverse(x)*(y-basis.centered_mean),epsabs=v00*epsrel*.05,**args)
    v11=_integral(eta,lambda x,y:inverse(x)*(y-basis.centered_mean)**2,epsabs=v00*epsrel*.05,**args)
    exact=np.array([[v00,v01],[v01,v11]])/basis.normalization**2
    if not np.all(np.isfinite(collision)) or not np.all(np.isfinite(exact)):
        raise ValueError('nonfinite electron transport response')
    for matrix in [collision,exact]:
        if np.any(np.diag(matrix)<=0):raise ValueError('nonpositive electron transport diagonal')
        scale=np.sqrt(np.diag(matrix));normalized=matrix/scale[:,None]/scale[None,:]
        try:np.linalg.cholesky(normalized)
        except np.linalg.LinAlgError as error:raise ValueError('electron transport is not positive') from error
    trial=np.linalg.inv(collision)
    heat_exact=exact[1,1]-exact[0,1]**2/exact[0,0]
    heat_trial=1/m11
    if min(heat_exact,heat_trial)<=0:raise ValueError('nonpositive zero-drift heat transport')
    return dict(eta=float(eta),current_normalization=basis.normalization,
                centered_current_mean=basis.centered_mean,current_variance=basis.variance,
                collision_moments=collision.tolist(),two_mode_mobility=trial.tolist(),
                exact_elastic_mobility=exact.tolist(),
                zero_drift_heat_response_two_mode=float(heat_trial),
                zero_drift_heat_response_exact=float(heat_exact),
                two_mode_to_exact_zero_drift_heat=float(heat_trial/heat_exact))


def born_response(eta,b_thermal, *, basis=None,tail=50.,epsrel=2e-11):
    """Screened first-Born drift/heat moments, with explicit approximation error."""
    if not np.isfinite(b_thermal) or b_thermal<=0:
        raise ValueError('positive thermal screening parameter required')
    result=collision_response(eta,lambda x:float(coulomb_bracket(b_thermal*x)),
                              basis=basis,tail=tail,epsrel=epsrel)
    result['b_thermal']=float(b_thermal)
    return result


def physical_prefactor(ion_density,charge):
    """Multiply dimensionless Born resistance by this cgs factor."""
    if not np.all(np.isfinite([ion_density,charge])) or min(ion_density,charge)<=0:
        raise ValueError('positive ion density and charge required')
    return 2*ME**2*(charge*E2)**2*ion_density/(3*math.pi*HBAR**3)
