#!/usr/bin/env python3
"""Finite-frequency plasma-transport comparisons and component bounds.

Keep the accepted comparison convention: absorption is divided by the
transverse refractive index and scattering is unchanged. This module does
not supply missing frequency tails, bound-electron dispersion, collision
damping, or a new scattering structure factor in the refractive medium.
"""
import math

import numpy as np
from scipy.special import roots_legendre

from audit_tops_electron_dispersion import index_squared, RW


def finite_means(u, remainder, known_scattering, cutoff, vstar2, order=16,
                 lower=None, upper=None, actual_absorption=None):
    """Bound the mean over a specified, fully supplied frequency interval.

    Remainder is total opacity minus known scattering, in cm2/g and including
    stimulated emission. Its unknown absorption fraction can lie anywhere
    from zero to one. All values are interpolated log-linearly in u; no
    source extrapolation is permitted. The normalization remains the full
    Rosseland weight RW, so returned means are finite-spectrum diagnostics.
    """
    u=np.asarray(u,dtype=float)
    components=np.column_stack([remainder,known_scattering])
    if (u.ndim!=1 or len(u)<2 or components.shape!=(len(u),2) or
            not np.isfinite(u).all() or np.any(u<=0) or np.any(np.diff(u)<=0) or
            not np.isfinite(components).all() or np.any(components<=0)):
        raise ValueError('nonpositive or unordered component data')
    if not math.isfinite(cutoff) or cutoff<=0 or not 0<=vstar2<.5:
        raise ValueError('invalid plasma parameters')
    if actual_absorption is not None:
        actual_absorption=np.asarray(actual_absorption,dtype=float)
        if (actual_absorption.shape!=u.shape or np.any(actual_absorption<0) or
                not np.isfinite(actual_absorption).all() or
                np.any(actual_absorption>components[:,0])):
            raise ValueError('absorption exceeds the unclassified opacity')
    lower=u[0] if lower is None else lower
    upper=u[-1] if upper is None else upper
    if not u[0]<=lower<upper<=u[-1] or cutoff>=upper:
        raise ValueError('requested integral leaves the supplied frequency interval')
    nodes,weights=roots_legendre(order)
    log_components=np.log(components)

    def integral(start,with_plasma):
        first=max(0,int(np.searchsorted(u,start,side='right'))-1)
        last=min(len(u)-2,int(np.searchsorted(u,upper,side='left'))-1)
        i=np.arange(first,last+1)
        lo=u[i].copy();hi=u[i+1].copy();lo[0]=start;hi[-1]=upper
        x=.5*((hi-lo)[:,None]*nodes+(hi+lo)[:,None])
        fraction=np.log(x/u[i,None])/np.log(u[i+1]/u[i])[:,None]
        k=np.exp(log_components[i,None,:]+fraction[:,:,None]*(log_components[i+1]-log_components[i])[:,None,:])
        wr=x**4*np.exp(-x)/(-np.expm1(-x))**2
        measure=.5*(hi-lo)[:,None]*weights*wr
        total=k[:,:,0]+k[:,:,1]
        if not with_plasma:
            return dict(inverse=float(np.sum(measure/total)),weight=float(measure.sum()))
        n=np.sqrt(index_squared(x,cutoff,vstar2))
        inverse_lower=n*n/total
        inverse_upper=n**3/(k[:,:,0]+n*k[:,:,1])
        result=dict(lower_inverse=float(np.sum(measure*inverse_lower)),
                    upper_inverse=float(np.sum(measure*inverse_upper)),
                    weight=float(measure.sum()))
        if actual_absorption is not None:
            # Interpolate the fraction linearly in log u. This remains in
            # [0,1] and preserves the component sum and pointwise bounds.
            fa=actual_absorption/components[:,0]
            f=fa[i,None]+fraction*(fa[i+1]-fa[i])[:,None]
            absorption=f*k[:,:,0]
            scattering=(1-f)*k[:,:,0]+k[:,:,1]
            inv=n**3/(absorption+n*scattering)
            if np.any(inv<inverse_upper*(1-1e-13)) or np.any(inv>inverse_lower*(1+1e-13)):
                raise ValueError('component bounds violated at a quadrature point')
            result['actual_inverse']=float(np.sum(measure*inv))
        return result

    base=integral(lower,False)
    plasma=integral(max(lower,cutoff),True)
    result=dict(uncut_finite_mean=RW/base['inverse'],
                scattering_extreme_mean=RW/plasma['lower_inverse'],
                absorption_extreme_mean=RW/plasma['upper_inverse'],
                finite_rosseland_weight_fraction=base['weight']/RW,
                propagating_rosseland_weight_fraction=plasma['weight']/RW,
                lower_frequency_u=float(lower),upper_frequency_u=float(upper),
                cutoff_u=float(cutoff),cutoff_below_supplied_interval=bool(cutoff<lower))
    if actual_absorption is not None:result['actual_component_mean']=RW/plasma['actual_inverse']
    result['ratio_lower']=result['scattering_extreme_mean']/result['uncut_finite_mean']
    result['ratio_upper']=result['absorption_extreme_mean']/result['uncut_finite_mean']
    result['component_interval_fraction']=result['ratio_upper']/result['ratio_lower']-1
    if result['ratio_lower']<1-1e-12 or result['ratio_upper']<result['ratio_lower']-1e-12:
        raise ValueError('invalid component interval or plasma enhancement')
    return result


def log_lagrange_weights(coordinates,target):
    x=np.log(np.asarray(coordinates,dtype=float));y=math.log(target)
    if np.any(np.diff(x)<=0) or not x[0]<=y<=x[-1]:
        raise ValueError('invalid thermodynamic interpolation coordinates')
    weights=np.ones(len(x))
    for i in range(len(x)):
        for j in range(len(x)):
            if i!=j:weights[i]*=(y-x[j])/(x[i]-x[j])
    if abs(weights.sum()-1)>1e-12:raise ValueError('interpolation does not preserve a constant')
    return weights


def interpolate_log_interval(weights,lower,upper):
    """Propagate a positive interval through a specified log interpolation.

    Cubic coefficients may be negative: simply interpolating the two endpoint
    models does not in general give bounds. Reverse their role for each
    negative coefficient instead. This bounds the interpolation operation,
    not its physical error between native states.
    """
    weights=np.asarray(weights,dtype=float)
    low=np.log(lower);high=np.log(upper)
    if (weights.shape!=low.shape or high.shape!=low.shape or
            not np.isfinite(weights).all() or not np.isfinite(low).all() or
            not np.isfinite(high).all() or np.any(high<low)):
        raise ValueError('invalid interpolation interval')
    lower_sum=np.sum(weights*np.where(weights>=0,low,high))
    upper_sum=np.sum(weights*np.where(weights>=0,high,low))
    return dict(lower=math.exp(lower_sum),upper=math.exp(upper_sum),
                weight_absolute_sum=float(np.abs(weights).sum()))
