"""Repulsive screened-ion collision integrals, Stanton & Murillo (2016).

Equations (25), (34)-(35), (42) and Appendix C, Tables III-IV.
https://doi.org/10.1103/PhysRevE.93.043203
Public manuscript: https://www.osti.gov/servlets/purl/1414086
Only classical, positively charged ion pairs are represented. These fits
do not supply attractive or quantum electron-ion collision integrals.
"""
import math

import numpy as np

KB=1.380649e-16
MU=1.66053906660e-24
ME=9.1093837015e-28
HBAR=1.054571817e-27
E2=2.3070775523417355e-19

# All printed source digits are retained; these are numerical source data.
CROSS={
    1:((.30031,-.69161,.59607,-.39822,-.20685),
       (.48516,1.66045,-.88687,.55990,1.65798,-1.02457)),
    2:((.40688,-.86425,.77461,-.34471,-.27626),
       (.83061,1.05229,-.59902,1.41500,.78874,-.48155))}
COLLISION={
    (1,1):((1.4660,-1.7836,1.4313,-.55833,.061162),
           (.081033,-.091336,.051760,-.50026,.17044)),
    (1,2):((.52094,.25153,-1.1337,1.2155,-.43784),
           (.20572,-.16536,.061572,-.12770,.066993)),
    (1,3):((.30346,.23739,-.62167,.56110,-.18046),
           (.68375,-.38459,.10711,.10649,.028760)),
    (2,2):((.85401,-.22898,-.60059,.80591,-.30555),
           (.43475,-.21147,.11116,.19665,.15195))}


def cross_section_fit(w, order):
    """Return phi_n = sigma_n/(2 pi lambda^2), equations C15-C18."""
    if order not in CROSS or not math.isfinite(w) or w<=0:
        raise ValueError('positive finite reduced speed and order 1 or 2 required')
    c,d=CROSS[order];v=math.log(w)
    if w<1:
        result=(c[0]+v*(c[1]+v*(c[2]+v*c[3])))/(1+c[4]*v)
    else:
        correction=(d[0]+v*(d[1]+v*(d[2]+v)))/(d[3]+v*(d[4]+v*(d[5]+v)))
        result=order*.5*math.log1p(w*w)/w**4*correction
    if not math.isfinite(result) or result<=0:
        raise ValueError('cross-section fit outside numerical support')
    return result


def reduced_integral(g, order, moment):
    """Return dimensionless K_nm, not the dimensional Burgers resistance K_ij."""
    if (order,moment) not in COLLISION or not math.isfinite(g) or g<=0:
        raise ValueError('positive finite collision strength and supported moment required')
    a,b=COLLISION[order,moment];v=math.log(g)
    if g<1:
        poly=a[0]+g*(a[1]+g*(a[2]+g*(a[3]+g*a[4])))
        result=-order*.25*math.factorial(moment-1)*(v+math.log(poly))
    else:
        result=(b[0]+v*(b[1]+v*b[2]))/(1+g*(b[3]+g*b[4]))
    if not math.isfinite(result) or result<=0:
        raise ValueError('collision fit outside numerical support')
    return result


def screening_length(number_density, charge, temperature, *, electron_stiffness_erg=None):
    """Effective static screening for classical ions, with explicit ion charges.

    Optional electron_stiffness_erg is (dP_e/dln n_e)/(n_e) at fixed T.
    Without it, use the source's nonrelativistic approximation, equation 25.
    This function does not choose ionization states or describe neutral atoms.
    """
    n,q=np.asarray(number_density,dtype=float),np.asarray(charge,dtype=float)
    if n.ndim!=1 or len(n)==0 or q.shape!=n.shape or np.any(n<0) or np.any(q<=0):
        raise ValueError('nonnegative ion populations and positive charges required')
    if not np.all(np.isfinite(n)) or not np.all(np.isfinite(q)) or not math.isfinite(temperature) or temperature<=0:
        raise ValueError('finite ion state and positive temperature required')
    ne=float(n@q)
    if ne<=0:
        raise ValueError('empty plasma')
    energy=KB*temperature
    ef=HBAR**2*(3*math.pi**2*ne)**(2/3)/(2*ME)
    stiffness=math.hypot(energy,2*ef/3) if electron_stiffness_erg is None else electron_stiffness_erg
    if not math.isfinite(stiffness) or stiffness<=0:
        raise ValueError('positive electron compressibility required')
    radii=(3*q/(4*math.pi*ne))**(1/3)
    gamma=q*q*E2/(radii*energy)
    inv_e=4*math.pi*E2*ne/stiffness
    inv_ion=4*math.pi*E2*n*q*q/(energy*(1+3*gamma))
    return dict(length_cm=1/math.sqrt(inv_e+float(sum(inv_ion))),
                electron_density_cm3=ne,electron_stiffness_erg=stiffness,
                electron_screening='source_equation_25' if electron_stiffness_erg is None else 'supplied_pressure_derivative',
                ion_sphere_coupling=gamma.tolist(),ion_sphere_radius_cm=radii.tolist())


def ion_pair(*, temperature, number_i, number_j, mass_i_u, mass_j_u,
             charge_i, charge_j, screening_cm):
    """Dimensional Omega_nm and Burgers coefficients for a repulsive ion pair.

    Omega has cm^3/s. Resistance has g/(cm^3 s). Mapping to Burgers uses
    Paxton et al. (2015), equation 86: Sigma_nm ratios equal Omega_nm ratios,
    and K_ij=(16/3) n_i n_j mu_ij Omega_11. Electrons are not accepted here.
    """
    values=(temperature,number_i,number_j,mass_i_u,mass_j_u,charge_i,charge_j,screening_cm)
    if any(not math.isfinite(v) or v<=0 for v in values) or min(mass_i_u,mass_j_u)<1:
        raise ValueError('positive classical ion pair required; no electron-ion continuation')
    energy=KB*temperature
    reduced_mass=MU*mass_i_u*mass_j_u/(mass_i_u+mass_j_u)
    interaction=charge_i*charge_j*E2
    strength=interaction/(screening_cm*energy)
    prefactor=math.sqrt(2*math.pi/reduced_mass)*interaction**2/energy**1.5
    omega={f'{n}{m}':prefactor*reduced_integral(strength,n,m) for n,m in COLLISION}
    r12,r13,r22=[omega[k]/omega['11'] for k in ('12','13','22')]
    return dict(collision_strength=strength,omega_cm3_s=omega,
                resistance_g_cm3_s=16/3*number_i*number_j*reduced_mass*omega['11'],
                z=1-.4*r12,zprime=2.5-2*r12+.4*r13,zdoubleprime=r22,
                reduced_mass_g=reduced_mass,
                relative_debroglie_over_screening=HBAR/math.sqrt(2*reduced_mass*energy)/screening_cm)
