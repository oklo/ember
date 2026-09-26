"""Partial-wave quantum scattering for a specified attractive screened potential.

The standalone C++ probe integrates the regular radial Schrodinger solution.
This module matches adjacent amplitudes to spherical Bessel functions and
sums the momentum-transfer cross section. No Born approximation, ion structure
factor, relativistic correction or physical screening selection is supplied.
"""
import math
import subprocess

import numpy as np
from scipy.special import spherical_jn, spherical_yn


def scatter_many(probe, queries):
    """Queries are (k*lambda, potential strength, radial step, rmax, lmax).

    lmax is the largest phase shift calculated; the transport sum has terms
    ell=0 through lmax-1. No missing partial-wave tail is silently appended.
    """
    text=''.join(' '.join(format(x,'.17g') if isinstance(x,float) else str(x) for x in q)+'\n'
                 for q in queries)
    result=subprocess.run([str(probe)],input=text,text=True,capture_output=True,
                          check=True,timeout=240)
    lines=result.stdout.splitlines()
    if len(lines)!=len(queries):raise ValueError('scattering probe returned incomplete results')
    out=[]
    for line,query in zip(lines,queries,strict=True):
        values=np.array(list(map(float,line.split())))
        a,c,h,rmax,ellmax=values[:5];ellmax=int(ellmax)
        if (len(values)!=5+2*(ellmax+1) or (a,c,rmax,ellmax)!=(query[0],query[1],query[3],query[4])):
            raise ValueError('scattering response differs from query')
        amplitudes=values[5:].reshape(-1,2);un,um=amplitudes.T
        ell=np.arange(ellmax+1);xn,xm=a*rmax,a*(rmax-h)
        jn,jm=xn*spherical_jn(ell,xn),xm*spherical_jn(ell,xm)
        yn,ym=xn*spherical_yn(ell,xn),xm*spherical_yn(ell,xm)
        # u proportional to jhat*cos(delta)-yhat*sin(delta).
        delta=np.arctan2(um*jn-un*jm,um*yn-un*ym)
        if not np.all(np.isfinite(delta)):raise ValueError('nonfinite matched phase')
        weights=np.arange(1,ellmax+1)*np.sin(np.diff(delta))**2
        cross=4*math.pi/a**2*float(np.sum(weights))
        out.append(dict(wavenumber=a,potential_strength=c,radial_step=h,
                        matching_radius=rmax,lmax=ellmax,phase_shifts=delta.tolist(),
                        transport_terms=weights.tolist(),cross_section_over_lambda2=cross))
    return out
