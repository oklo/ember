#!/usr/bin/env python3
"""Independent constant-pressure and reference-energy thermal transport checks.

Analytic ideal ions plus classical/polytropic electrons supply bulk enthalpy
and chemical potentials independently of the tested EOS derivative identity.
Transport matrices here are synthetic; no stellar transport is accepted.
"""
import argparse
from datetime import datetime, timezone
import hashlib, json
from pathlib import Path
import numpy as np
from audit_diffusion_composition_forces import AnalyticMixture, R, KB, MU
from diffusion_composition_forces import chemical_acceleration
from diffusion_thermal_transport import exchange_enthalpy, HeatFluxConvention


def bulk_enthalpy(eos, T, rho, x):
    n, ne, _ = eos.particles(T, rho, x)
    return (2.5 * KB * T * (sum(n) + ne) +
            eos.gamma / (eos.gamma - 1) * eos.K * ne**eos.gamma) / rho


def temperature_composition_gradient(eos, T, rho, x, indices, reference):
    _, ne, _ = eos.particles(T, rho, x)
    a, z = eos.a, eos.z
    alpha = z[indices]/a[indices] - z[reference]/a[reference]
    # Independent derivative of the analytic material energy at fixed T,rho.
    ec = 1.5 * R * T * (1/a[indices] - 1/a[reference] + alpha)
    ec += eos.gamma/(eos.gamma-1)*eos.K*ne**(eos.gamma-1)/MU*alpha
    return -ec/T


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=Path)
    args = parser.parse_args()
    checks, cases = [], []
    def check(name, actual, expected, scale, tolerance=2e-11, **meta):
        error=float(np.max(np.abs(np.asarray(actual)-expected))/scale)
        checks.append(dict(name=name, normalized_error=error, tolerance=tolerance,
                           passed=bool(np.isfinite(error) and error<=tolerance), **meta))
    rng=np.random.default_rng(68423)
    for masses, charges in [([1.,4.],[1.,2.]),([1.,3.,4.],[1.,2.,2.]),
                            ([1.,3.,4.,12.,56.],[1.,2.,2.,6.,26.])]:
        a,z=np.array(masses),np.array(charges);nions=len(a)
        for hydrogen in [.001,.15,.65]:
            x=np.arange(1,nions+1,dtype=float);x[0]=0
            x*=(1-hydrogen)/sum(x);x[0]=hydrogen
            T,rho=1.1e7,1300.
            for excess in [0.,.5,5.,50.,5000.]:
                eos=AnalyticMixture(a,z,T,rho,x,excess);P=eos.pressure(T,rho,x)
                for ref in [0,nions-1]:
                    ind=np.array([i for i in range(nions) if i!=ref]);nc=len(ind)
                    hessian,gr,delta,pr,pt=eos.derivatives(T,rho,x,ind,ref)
                    gt=temperature_composition_gradient(eos,T,rho,x,ind,ref)
                    enthalpy=exchange_enthalpy(temperature=T,gradient_phi_lnT=gt,
                                              gradient_phi_lnrho=gr,material_delta=delta)
                    scale=R*T+max(abs(enthalpy))
                    meta=dict(ions=nions,hydrogen=hydrogen,electron_excess=excess,reference=ref)
                    measured=[]
                    for factor in [1.,.5]:
                        row=[]
                        for i in ind:
                            step=factor*min(1e-4,.1*x[i],.1*x[ref])
                            def value(offset):
                                xx=x.copy();xx[i]+=offset;xx[ref]-=offset
                                rr=eos.rho_at_pressure(T,xx,P,rho)
                                return bulk_enthalpy(eos,T,rr,xx)
                            row.append((value(-2*step)-8*value(-step)+8*value(step)-value(2*step))/(12*step))
                        measured.append(np.array(row))
                    check('independent bulk enthalpy at fixed material pressure',measured[-1],enthalpy,scale,2e-6,**meta)
                    check('bulk-enthalpy composition spacing',measured[0],measured[1],scale,2e-6,**meta)
                    q=eos.chemical_exchange(T,rho,x,ind,ref)
                    for step in [1e-4,5e-5]:
                        values=[]
                        for sign in [1,-1]:
                            tt=T*np.exp(sign*step);rr=eos.rho_at_pressure(tt,x,P,rho)
                            values.append(eos.chemical_exchange(tt,rr,x,ind,ref))
                        hq=q-(values[0]-values[1])/(2*step)
                        check('independent fixed-pressure chemical temperature derivative',hq,enthalpy,scale,2e-6,step=step,**meta)
                    # An independent radial path checks force sign and factors
                    # of T, not only the matrix algebra below.
                    length=1e9;lt=-.15/length;lr=-.25/length
                    dc=.02*x[ind]*np.sin(ind+1)/length
                    dx=np.zeros(nions);dx[ind]=dc;dx[ref]=-sum(dc)
                    driving=chemical_acceleration(temperature=T,hessian_phi=hessian,
                        gradient_phi_lnrho=gr,material_delta=delta,composition_gradient=dc,
                        logarithmic_density_gradient=lr,logarithmic_temperature_gradient=lt)
                    distance=length*2e-5;phi_values=[]
                    for sign in [1,-1]:
                        tt=T*np.exp(sign*lt*distance);rr=rho*np.exp(sign*lr*distance)
                        xx=x+sign*distance*dx
                        phi_values.append(eos.chemical_exchange(tt,rr,xx,ind,ref)/tt)
                    full_force=-(phi_values[0]-phi_values[1])/(2*distance)
                    e0=scale;theta=-lt/T
                    convention=HeatFluxConvention(enthalpy,e0)
                    f=np.append(full_force,e0*theta)
                    reduced=convention.reduced_forces(f)
                    check('independent radial chemical force transform',reduced[:-1],-driving/T,R/length+max(abs(driving/T)),2e-6,**meta)
                    # An arbitrary reciprocal matrix with nonzero cross terms.
                    b=rng.normal(size=(nc+1,nc+1));l=b@b.T+np.eye(nc+1)*.2
                    lp=convention.reduced_matrix(l);j=l@f;jp=convention.reduced_fluxes(j)
                    check('coupled flux identity',lp@reduced,jp,max(abs(jp)),**meta)
                    check('entropy production invariance',jp@reduced,j@f,abs(j@f),**meta)
                    check('Onsager symmetry',lp,lp.T,max(abs(lp.flatten())),**meta)
                    check('nonnegative transformed entropy matrix',max(0.,-min(np.linalg.eigvalsh(lp))),0.,1.,**meta)
                    check('force inverse',convention.full_forces(reduced),f,max(abs(f)),**meta)
                    check('flux inverse',convention.full_fluxes(jp),j,max(abs(j)),**meta)
                    check('matrix inverse',convention.full_matrix(lp),l,max(abs(l.flatten())),**meta)
                    # Add species-specific constant energy and entropy references.
                    de=e0*rng.normal(size=nc);ds=R*rng.normal(size=nc)
                    shifted_gt=gt-de/T  # the constant entropy term has zero lnT response
                    shifted_h=exchange_enthalpy(temperature=T,gradient_phi_lnT=shifted_gt,
                        gradient_phi_lnrho=gr,material_delta=delta)
                    shifted=HeatFluxConvention(shifted_h,e0)
                    check('energy-reference shift of exchange enthalpy',shifted_h,enthalpy+de,scale,**meta)
                    fshift=f.copy();fshift[:-1]-=de*theta
                    jshift=j.copy();jshift[-1]+=de@j[:-1]/e0
                    check('reference-independent reduced force',shifted.reduced_forces(fshift),reduced,max(abs(reduced)),**meta)
                    check('reference-independent reduced heat flux',shifted.reduced_fluxes(jshift),jp,max(abs(jp)),**meta)
                    ar=np.eye(nc+1);ar[-1,:-1]=de/e0
                    check('reference-independent reduced kinetic matrix',shifted.reduced_matrix(ar@l@ar.T),lp,max(abs(lp.flatten())),**meta)
                    # Constant entropy reference cancels in q-T*dq/dT explicitly.
                    step=1e-4;tv=[]
                    for sign in [1,-1]:
                        tt=T*np.exp(sign*step);rr=eos.rho_at_pressure(tt,x,P,rho)
                        tv.append(eos.chemical_exchange(tt,rr,x,ind,ref)+de-tt*ds)
                    hshift=q+de-T*ds-(tv[0]-tv[1])/(2*step)
                    check('independent energy and entropy reference derivative',hshift,enthalpy+de,scale,2e-6,**meta)
                    # Change only numerical energy units, leaving physical flux fixed.
                    alt=HeatFluxConvention(enthalpy,7*e0)
                    f_alt=f.copy();f_alt[-1]*=7;j_alt=j.copy();j_alt[-1]/=7
                    check('energy-unit independent chemical forces',alt.reduced_forces(f_alt)[:-1],reduced[:-1],max(abs(reduced)),**meta)
                    check('energy-unit independent reduced heat',alt.reduced_fluxes(j_alt)[-1]*7,jp[-1],max(abs(jp)),**meta)
                    cases.append(meta)
    valid=dict(temperature=1e7,gradient_phi_lnT=[1.,2.],gradient_phi_lnrho=[3.,4.],material_delta=.1)
    guards=0
    for patch in [dict(temperature=0),dict(gradient_phi_lnT=[np.nan,2]),dict(gradient_phi_lnrho=[1]),dict(material_delta=np.inf)]:
        try:exchange_enthalpy(**{**valid,**patch})
        except ValueError:guards+=1
    for h,e0 in [([],1),([1,np.nan],1),([1],0)]:
        try:HeatFluxConvention(h,e0)
        except ValueError:guards+=1
    c=HeatFluxConvention([1.,2.],3.)
    for fn,value in [(c.reduced_fluxes,[1,2]),(c.reduced_forces,[1,2,np.nan]),(c.reduced_matrix,np.eye(2))]:
        try:fn(value)
        except ValueError:guards+=1
    check('invalid input guards',guards,10,1.,0.)
    passed=all(c['passed'] for c in checks)
    sources=[Path(__file__),Path('scripts/diffusion_thermal_transport.py'),Path('scripts/diffusion_composition_forces.py'),Path('scripts/audit_diffusion_composition_forces.py')]
    report=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed' if passed else 'failed',
                accepted_for_stellar_evolution=False,cases=cases,checks=checks,guard_count=guards,
                inputs_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                remaining_work=['Physical EOS enthalpy identities','Physical ion/electron thermal coefficients in the same convention','Conservative stellar species and energy transport'])
    with args.report.open('x') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(outcome=report['outcome'],cases=len(cases),checks=len(checks),failures=[c for c in checks if not c['passed']])))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
