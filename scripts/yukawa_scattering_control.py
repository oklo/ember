"""Direct classical repulsive Yukawa scattering; independent of fitted data.

The turning radius replaces the impact parameter in the cross-section
integral. A Coulomb potential with the same value at the turning point is
subtracted analytically to avoid cancellation at small scattering angles.
No Stanton-Murillo fit coefficients are imported here.
"""
from functools import lru_cache
import math

import numpy as np
from scipy.special import lambertw, roots_genlaguerre, roots_legendre


@lru_cache(maxsize=16)
def gauss(order):
    nodes,weights=roots_legendre(order)
    return (nodes+1)/2,weights/2


def cross_sections(w, *, angle_order=128, radius_order=256, tail_radii=32.):
    if not math.isfinite(w) or w<=0:
        raise ValueError('positive reduced speed required')
    rc=float(lambertw(1/(w*w)).real)
    point,weight=gauss(radius_order)
    y=point*math.log1p(tail_radii/rc)
    r=rc*np.exp(y)
    u0=np.exp(-r)/(w*w*r)
    if np.any(u0<=0) or np.any(u0>=1):
        raise ValueError('invalid repulsive turning point')
    gamma=np.sqrt(1-u0)
    s,sw=gauss(angle_order)
    u=1-s*s
    d0=s[None,:]**2*(1+(1-u0[:,None])*u[None,:])
    delta=u0[:,None]*u[None,:]*(-np.expm1(-r[:,None]*s[None,:]**2/u[None,:]))
    root0,root=np.sqrt(d0),np.sqrt(d0+delta)
    difference=delta/(root0*root*(root+root0))
    theta=2*np.arcsin(u0/(2-u0))+4*gamma*np.sum(sw[None,:]*s[None,:]*difference,axis=1)
    jacobian=r*(r+.5*(r-1)*np.exp(-r)/(w*w))*math.log1p(tail_radii/rc)
    return np.array([np.dot(weight,jacobian*(2*np.sin(theta/2)**2)),
                     np.dot(weight,jacobian*np.sin(theta)**2)])


def collision_integrals(g, *, thermal_order=64, angle_order=128, radius_order=256):
    """All four CE moments with shared scattering evaluations, x=g*w^2."""
    if not math.isfinite(g) or g<=0:
        raise ValueError('positive collision strength required')
    nodes,weights=roots_genlaguerre(thermal_order,2)
    phi=np.array([cross_sections(math.sqrt(x/g),angle_order=angle_order,radius_order=radius_order) for x in nodes])
    return {(n,m):float(np.dot(weights*nodes**(m-1),phi[:,n-1])/(2*g*g))
            for n,m in ((1,1),(1,2),(1,3),(2,2))}
