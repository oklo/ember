"""H1/He3/He4 mixing terms and an offline tensor interpolation of F/T.

Z is fixed. FreeEOS treats all helium electronically as He4; the isotope
entropy is restored explicitly. Radiation is composition independent.
This module is an offline interpolation calculation, not a production EOS.
Version 2 corrects the isotope entropy reference to Ember's atomic masses;
version 1 and the source evaluations remain unchanged.
"""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import xlogy

RGAS=6.02214076e23*1.380649e-16
# Physical translational masses match nuclides[] in composition.hpp.
# Baryon counts in populations() still use the integer mass numbers.
ISOTOPE_MASS_RATIO=3.01602932/4.00260325


def populations(x,y,z=.02):
    if not np.isfinite(x+y+z) or min(x,y,z)<0 or x+y+z>=1:
        raise ValueError('finite nonnegative composition with positive He4 required')
    return x,y/3,(1-z-x-y)/4


def mixing_value(x,y,z=.02):
    n1,n3,n4=populations(x,y,z)
    return float(RGAS*(xlogy(n1,n1)+xlogy(n3,n3)+xlogy(n4,n4)
                       -1.5*n3*np.log(ISOTOPE_MASS_RATIO)))


def source_mixing_value(x,y,z=.02):
    n1,n3,n4=populations(x,y,z);he=n3+n4
    return float(RGAS*(xlogy(n1,n1)+xlogy(he,he)))


def isotope_value(x,y,z=.02):
    n1,n3,n4=populations(x,y,z);he=n3+n4
    return float(RGAS*(xlogy(n3,n3/he)+xlogy(n4,n4/he)-1.5*n3*np.log(ISOTOPE_MASS_RATIO)))


def mixing_derivatives(x,y,z=.02):
    n1,n3,n4=populations(x,y,z)
    if n1==0 or n3==0:
        raise ValueError('zero populations require an explicit active-species limit')
    # B maps (dX,dY3) to changes in nuclei per baryonic mass unit.
    b=np.array([[1.,0.],[0.,1/3],[-.25,-.25]])
    n=np.array([n1,n3,n4])
    grad=RGAS*(b.T@(np.log(n)+1)+np.array([0.,-.5*np.log(ISOTOPE_MASS_RATIO)]))
    hess=RGAS*((b.T/n)@b)
    return grad,hess


class PotentialSpline:
    """Tensor cubic interpolation of the smooth residual potential.

    CubicSpline uses a line for two nodes and a parabola for three. The same
    linear source weights can act on all thermodynamic derivative tables.
    Construction is deliberately offline; this is not a per-query design.
    """
    def __init__(self,x,y,residual):
        self.x=np.asarray(x);self.y=np.asarray(y);self.values=np.asarray(residual)
        if self.values.shape!=(len(x),len(y)) or not np.all(np.isfinite(self.values)):
            raise ValueError('finite rectangular potential required')
        self.bx=CubicSpline(self.x,np.eye(len(x)),axis=0,extrapolate=False)
        self.by=CubicSpline(self.y,np.eye(len(y)),axis=0,extrapolate=False)

    def jet(self,x,y):
        if not(self.x[0]<=x<=self.x[-1] and self.y[0]<=y<=self.y[-1]):
            raise ValueError('composition outside interpolation grid')
        w=[self.bx(x,k) for k in range(3)]
        v=[self.by(y,k) for k in range(3)]
        value=float(w[0]@self.values@v[0])+mixing_value(x,y)
        grad,hess=mixing_derivatives(x,y)
        grad=grad+np.array([w[1]@self.values@v[0],w[0]@self.values@v[1]])
        hess=hess+np.array([[w[2]@self.values@v[0],w[1]@self.values@v[1]],
                            [w[1]@self.values@v[1],w[0]@self.values@v[2]]])
        return value,grad,hess
