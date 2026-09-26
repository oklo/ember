"""Conditional trace-ion transport and conservative scalar redistribution.

Published Thoul collision coefficients, dimensional Burgers equations,
fully stripped ions, prescribed H/He background. Each metal is a separate
trace species; neither opacity nor the background responds to its motion.
"""
import math
import numpy as np
from scipy.linalg import solve_banded

from diffusion_burgers import KB, MU, ME, solve
from ion_collision_integrals import E2

ELEMENTS = [("C",12,6),("N",14,7),("O",16,8),("Na",23,11),
            ("Mg",24,12),("Si",28,14),("Ca",40,20),("Fe",56,26)]


def coefficients(lo, hi, mass_number, charge, heat_flow, trace_fraction=1e-10):
    """Return w = drift - D dln(X_trace)/dr in cm/s, cm²/s.

    States have columns m, r, rho, T, L, XH, X3. H/He fractions are
    normalized to one after removing the fixed metal background. The trace
    fraction is a limit control, not an abundance floor in the evolution.
    """
    lo,hi = np.asarray(lo),np.asarray(hi)
    radius=.5*(lo[1]+hi[1]);rho=.5*(lo[2]+hi[2]);T=math.sqrt(lo[3]*hi[3])
    dr=(hi[0]-lo[0])/(4*math.pi*radius**2*rho)
    g=6.6743e-8*.5*(lo[0]+hi[0])/radius**2
    assert dr>0 and min(lo[3],hi[3])>=2e6
    def background(p):
        return np.array([p[5],p[6],.98-p[5]-p[6]])/.98
    left,right=background(lo),background(hi)
    assert np.all(left>0) and np.all(right>0)
    f=np.r_[.5*(left+right)*(1-trace_fraction),trace_fraction]
    a=np.array([1.,3.,4.,float(mass_number),ME/MU])
    z=np.array([1.,2.,2.,float(charge),-1.])
    ni=rho/MU*f/a[:-1];n=np.r_[ni,ni@z[:-1]]
    # Same effective Coulomb logarithm as the corrected Thoul README;
    # dimensional constants and measured number densities are explicit.
    screening=max(math.sqrt(KB*T/(4*math.pi*E2*(n@(z*z)))),
                  (3/(4*math.pi*sum(ni)))**(1/3))
    argument=4*KB*T*screening/(E2*abs(np.outer(z,z)))
    cl=.81245*np.log1p(.18769*argument**1.2)
    reduced=MU*np.outer(a,a)/(a[:,None]+a[None,:])
    resistance=4/3*np.sqrt(2*np.pi*reduced)*E2**2*np.outer(n,n)*np.outer(z*z,z*z)/(KB*T)**1.5*cl
    resistance=.5*(resistance+resistance.T)
    gp=(math.log(hi[2]/lo[2])+math.log(hi[3]/lo[3]))/dr
    dpi=ni*KB*T*np.r_[gp+np.log(right/left)/dr,gp]
    kw=dict(number_density=n,mass_u=a,charge=z,ion_pressure_gradient=dpi,
            gravity=g,temperature=T,temperature_gradient=T*math.log(hi[3]/lo[3])/dr,
            resistance=resistance,z=np.full((5,5),.6),zprime=np.full((5,5),1.3),
            zdoubleprime=np.full((5,5),2.),heat_flow=heat_flow)
    zero=solve(**kw)
    slope=1e-8
    changed=dpi.copy();changed[3]+=ni[3]*KB*T*slope
    one=solve(**{**kw,"ion_pressure_gradient":changed})
    drift=float(zero.velocity_cm_s[3]);D=float(-(one.velocity_cm_s[3]-drift)/slope)
    assert math.isfinite(drift) and math.isfinite(D) and D>0
    return drift,D


def mass_weights(mass):
    mass=np.asarray(mass,dtype=float)
    assert np.all(np.diff(mass)>0) and mass[0]>0
    w=np.zeros(len(mass));w[0]=mass[0]
    w[:-1]+=.5*np.diff(mass);w[1:]+=.5*np.diff(mass)
    return w


def bernoulli(x):
    """x/(exp(x)-1), including zero and large positive/negative limits."""
    x=np.asarray(x,dtype=float);out=np.empty_like(x)
    small=abs(x)<1e-4;pos=x>50;neg=x < -50;mid=~(small|pos|neg)
    out[small]=1-x[small]/2+x[small]**2/12-x[small]**4/720
    out[pos]=x[pos]*np.exp(-x[pos])/(1-np.exp(-x[pos]))
    out[neg]=-x[neg]/(1-np.exp(x[neg]))
    out[mid]=x[mid]/np.expm1(x[mid])
    return out


def step(fractions, weights, forward, backward, regions, seconds):
    """Implicit conservative drift/diffusion with exact convective mixing.

    The outward face flux is forward*X_left - backward*X_right, in g/s.
    Coefficients must be nonnegative. Regions partition nodal cells; a
    multi-cell region represents instantaneous mass-weighted mixing.
    Reflecting center and surface: no tracer mass leaves the star.
    """
    x,w,u,v=[np.asarray(a,dtype=float) for a in (fractions,weights,forward,backward)]
    if x.shape!=w.shape or u.shape!=v.shape or u.shape!=(len(x)-1,):
        raise ValueError("inconsistent transport arrays")
    if not all(np.isfinite(a).all() for a in (x,w,u,v)) or min(x)<0 or min(w)<=0 or min(u)<0 or min(v)<0 or seconds<=0:
        raise ValueError("invalid composition, weights, rates or interval")
    if not regions or regions[0][0]!=0 or regions[-1][1]!=len(x):
        raise ValueError("incomplete mixing partition")
    for i,(a,b) in enumerate(regions):
        if b<=a or (i and a!=regions[i-1][1]):raise ValueError("invalid mixing partition")
    starts=np.array([a for a,b in regions]);ends=np.array([b for a,b in regions])
    masses=np.add.reduceat(w,starts);inventory=np.add.reduceat(w*x,starts)
    interfaces=ends[:-1]-1
    f=u[interfaces]*seconds;b=v[interfaces]*seconds
    diagonal=masses.copy();diagonal[:-1]+=f;diagonal[1:]+=b
    matrix=np.zeros((3,len(regions)));matrix[1]=diagonal
    matrix[0,1:]=-b;matrix[2,:-1]=-f
    answer=solve_banded((1,1),matrix,inventory,check_finite=False)
    out=np.repeat(answer,ends-starts)
    if not np.isfinite(out).all() or min(out)<0:raise ValueError("negative or nonfinite implicit result")
    total=float(w@x)
    mass_error=abs(float(w@out)-total)/max(total,np.finfo(float).tiny)
    if mass_error>1e-10:raise ValueError(f"tracer mass residual {mass_error}")
    flux=f*answer[:-1]-b*answer[1:]
    residual=masses*answer-inventory
    residual[:-1]+=flux;residual[1:]-=flux
    local_error=max(abs(residual))/max(total,np.finfo(float).tiny)
    if local_error>1e-10:raise ValueError(f"cell inventory residual {local_error}")
    return out,mass_error
