"""Offline three-moment Born/Pauli electron collision comparison.

Full scattering angles and final-state Pauli factors are retained. Static
Yukawa screening, nonrelativistic energies and a first Born amplitude are
prescribed approximations, not an accepted stellar transport model.
See docs/ELECTRON_ELECTRON_COLLISIONS.md for counting and normalization.
"""
from dataclasses import dataclass
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.special import expit
from electron_ion_heat import current_basis, _integral
from electron_ion_born import coulomb_bracket


@dataclass(frozen=True)
class EnergyBasis:
    eta: float
    normalization: float
    mean: float
    variance: float
    skew_ratio: float
    quadratic_norm: float
    statistics: str

    def polynomials(self, energy):
        s = np.asarray(energy)-self.eta-self.mean
        return np.stack([np.ones_like(s), s/self.variance,
                         (s*s-self.skew_ratio*s-self.variance)/self.quadratic_norm])


def energy_basis(eta, *, statistics='fermi'):
    if not np.isfinite(eta) or not -30 <= eta <= 4096:
        raise ValueError('source integration supports -30 <= eta <= 4096 only')
    if statistics == 'maxwell':
        return EnergyBasis(eta, .75*math.sqrt(math.pi)*math.exp(eta),
                           2.5-eta, 2.5, 2., math.sqrt(17.5), statistics)
    if statistics != 'fermi':
        raise ValueError('fermi or maxwell statistics required')
    b = current_basis(eta)
    def moment(n):
        return _integral(eta, lambda x,y:x**1.5*(y-b.centered_mean)**n,
                         tail=50.,epsrel=2e-11)/b.normalization
    skew, fourth = moment(3), moment(4)
    norm2 = fourth-skew*skew/b.variance-b.variance**2
    if not norm2 > 0: raise ValueError('nonpositive polynomial norm')
    return EnergyBasis(eta,b.normalization,b.centered_mean,b.variance,
                       skew/b.variance,math.sqrt(norm2),statistics)


def rule(count, low, high):
    x,w = leggauss(count)
    return low+(x+1)*(high-low)/2, w*(high-low)/2


def born_electron_ion_moments(basis, b_thermal):
    """Three-mode ion collision matrix with the existing dimensionless scale."""
    if not np.isfinite(b_thermal) or b_thermal <= 0:
        raise ValueError('positive finite screening parameter required')
    eta=basis.eta
    def factor(x):
        occupation = (expit(eta-x)*expit(x-eta) if basis.statistics=='fermi'
                      else math.exp(eta-x))
        return float(coulomb_bracket(b_thermal*x))*occupation
    matrix=np.empty((3,3))
    for i in range(3):
        for j in range(i,3):
            value=quad(lambda x:factor(x)*np.prod(basis.polynomials(x)[[i,j]]),
                       0,max(eta,0)+50,epsrel=2e-10,epsabs=1e-100,
                       points=[max(eta,0)+1],limit=200)[0]
            matrix[i,j]=matrix[j,i]=value
    return matrix


def collision_matrix(eta, b_thermal, *, orders=(48,24,24,32,32),
                     statistics='fermi', exchange=True, scattering='born', tail=40.):
    """Return M_ee / C_ref; C_ref = 2 m_e^2 e^4 n_e/(3 pi hbar^3).

    Tensor quadrature in pair energy A, mixing angle psi, cos(alpha),
    logarithmically spaced scattering angle and scattering azimuth. Energy
    modes are P0=1, P1=s/variance and an orthogonal quadratic P2. Pauli factors
    are evaluated on all four momenta. The matrix's first row/column vanish
    exactly by momentum conservation. The two other modes are coupled.

    For unpolarized identical electrons the amplitude is integrated on a
    hemisphere. With exchange=False, direct+exchanged cross sections on that
    hemisphere recover a distinguishable full-sphere integral. The constant
    cross-section option is a separate classical normalization control.
    """
    if not np.isfinite(b_thermal) or b_thermal<=0:
        raise ValueError('positive finite screening parameter required')
    if len(orders)!=5 or any(not isinstance(n,int) or not 8<=n<=160 for n in orders):
        raise ValueError('five bounded integer quadrature orders required')
    if not np.isfinite(tail) or not 35<=tail<=60:
        raise ValueError('controlled energy tail required')
    if scattering not in ['born','constant']:
        raise ValueError('born or constant scattering required')
    if scattering=='constant' and exchange:
        raise ValueError('constant cross section is a distinguishable control')
    basis=energy_basis(eta,statistics=statistics)
    na,npair,nu,nt,nphi=orders
    energies,wa=rule(na,0,max(eta,0)+tail)
    angles,wpair=rule(npair,0,math.pi/2)
    u,wu=rule(nu,-1,1);u=u[:,None,None]
    su=np.sqrt(1-u*u)
    phi=(np.arange(nphi)+.5)*(2*math.pi/nphi)
    cp=np.cos(phi)[None,None,:]
    unit,wt=rule(nt,0,1);unit=unit[None,:,None]
    angular_weight=wu[:,None,None]*wt[None,:,None]*(2*math.pi/nphi)
    screening=4/b_thermal
    result=np.zeros((2,2));evaluations=0
    for a,w_a in zip(energies,wa):
        r=math.sqrt(a)
        for psi,w_psi in zip(angles,wpair):
            co,si=math.cos(psi),math.sin(psi)
            G,g=r*co,2*r*si
            g2=g*g
            length=math.log1p(g2/(2*screening))
            t=unit*length
            y=(screening/g2)*np.expm1(t)
            ct=1-2*y;st=2*np.sqrt(y*(1-y))
            jacobian=2*screening/g2*np.exp(t)*length
            if scattering=='born':
                direct=np.exp(-t)/screening
                other=1/(screening+g2*(1-y))
                cross=direct*direct+other*other
                if exchange: cross=cross-direct*other
            else:
                cross=2*np.ones_like(t)
            b=G*g*u
            bp=G*g*(u*ct+su*st*cp)
            if statistics=='fermi':
                weight=expit(eta-a-b)*expit(eta-a+b)*expit(a+bp-eta)*expit(a-bp-eta)
            else:
                weight=math.exp(2*(eta-a))
            # Delta[c*c^2] = (G.g)g - (G.g')g'. Its squared norm is
            # expressed without cancellation at small scattering angles.
            d2=g2*((b-bp*ct)**2+(bp*st)**2)
            h=(b-bp)*(b+bp)
            z=2*(a-eta-basis.mean)-basis.skew_ratio
            f11=d2/basis.variance**2
            f12=(z*d2+2*h*h)/(basis.variance*basis.quadratic_norm)
            f22=(z*z*d2+4*(z+G*G)*h*h)/basis.quadratic_norm**2
            weighted=angular_weight*jacobian*cross*weight
            radial=8*a**2.5*co*co*si**3*w_a*w_psi
            result[0,0]+=radial*np.sum(weighted*f11)
            result[0,1]+=radial*np.sum(weighted*f12)
            result[1,1]+=radial*np.sum(weighted*f22)
            evaluations+=nu*nt*nphi
    result[1,0]=result[0,1]
    matrix=np.zeros((3,3))
    matrix[1:,1:]=3/(8*math.pi*basis.normalization)*result
    if not np.all(np.isfinite(matrix)):
        raise ValueError('nonfinite electron collision integral')
    np.linalg.cholesky(matrix[1:,1:])
    return dict(eta=eta,b_thermal=b_thermal,statistics=statistics,exchange=exchange,
                scattering=scattering,orders=list(orders),tail=tail,
                collision_matrix=matrix.tolist(),evaluations=evaluations,
                basis=vars(basis))
