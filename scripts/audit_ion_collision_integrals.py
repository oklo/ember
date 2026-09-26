#!/usr/bin/env python3
"""Check published repulsive-ion fits against independent scattering integrals.

Fit comparisons concern classical Yukawa scattering, not the accuracy of
that effective potential for a stellar plasma. No electron-ion rates or
stellar abundance evolution are accepted by these checks.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from fetch_tops_composition import digest
from ion_collision_integrals import (KB,MU,E2,CROSS,COLLISION,cross_section_fit,
                                     reduced_integral,screening_length,ion_pair)
from yukawa_scattering_control import cross_sections,collision_integrals


def integrated_fit(g,n,m):
    def integrand(x):
        return math.exp(-x)*x**(m+1)*cross_section_fit(math.sqrt(x/g),n)/(2*g*g)
    bounds=sorted({0.,min(g,80.),80.})
    parts=[quad(integrand,lo,hi,epsabs=1e-13,epsrel=1e-8,limit=160) for lo,hi in zip(bounds,bounds[1:])]
    return sum(p[0] for p in parts),sum(p[1] for p in parts)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('preserve collision comparisons')
    source=Path('/tmp/ember-diffusion-collision-literature-v1/stanton_murillo_2016.pdf')
    receipt=source.with_suffix('.pdf.receipt.json')
    assert digest(source)==json.loads(receipt.read_text())['sha256']
    paths=[Path(__file__),Path('scripts/ion_collision_integrals.py'),
           Path('scripts/yukawa_scattering_control.py'),source,receipt,
           Path('/tmp/ember-diffusion-collision-literature-v1/paxton_2015.pdf')]
    inputs={str(p.resolve()):digest(p) for p in paths}
    checks=[];cross=[];moments=[]
    def check(name,value,tolerance):
        checks.append(dict(name=name,error=float(value),tolerance=tolerance,
                           pass_check=bool(math.isfinite(value) and value<=tolerance)))
    for w in (.02,.05,.1,.2,.5,.8,1.,1.25,2.,5.,10.,16.,30.,60.):
        coarse=cross_sections(w,angle_order=128,radius_order=256)
        fine=cross_sections(w,angle_order=256,radius_order=512)
        tail=cross_sections(w,angle_order=128,radius_order=256,tail_radii=40.)
        for n in (1,2):
            fitted=cross_section_fit(w,n);reference=fine[n-1]
            row=dict(speed=w,order=n,direct=reference,fit=fitted,
                     quadrature_relative_change=abs(coarse[n-1]/reference-1),
                     tail_relative_change=abs(tail[n-1]/coarse[n-1]-1),
                     fit_relative_error=abs(fitted/reference-1))
            cross.append(row)
            check('direct cross-section convergence',row['quadrature_relative_change'],2e-5)
            check('direct cross-section tail',row['tail_relative_change'],2e-5)
            check('cross-section fit',row['fit_relative_error'],.005)
        print(json.dumps({'completed_cross_section_speed':w}),flush=True)
    for g in (.05,.3,1.,3.,30.):
        coarse=collision_integrals(g,thermal_order=64,angle_order=96,radius_order=192)
        fine=collision_integrals(g,thermal_order=128,angle_order=128,radius_order=256)
        for n,m in COLLISION:
            fitted=reduced_integral(g,n,m)
            scalar_integral,error=integrated_fit(g,n,m)
            row=dict(strength=g,order=n,moment=m,direct=fine[n,m],fit=fitted,
                     integrated_cross_section_fit=scalar_integral,
                     quadrature_relative_change=abs(coarse[n,m]/fine[n,m]-1),
                     fit_relative_error=abs(fitted/fine[n,m]-1),
                     integrated_cross_section_fit_relative_error=abs(scalar_integral/fine[n,m]-1))
            moments.append(row)
            check('direct collision-integral convergence',row['quadrature_relative_change'],2e-5)
            check('collision-integral fit',row['fit_relative_error'],.02)
            check('cross-section integral against direct scattering',row['integrated_cross_section_fit_relative_error'],.005)
        print(json.dumps({'completed_collision_strength':g}),flush=True)
    # Independent SI dimensional evaluation, with physical constants in SI.
    state=dict(temperature=1.3e7,number_i=1e25,number_j=1e27,
               mass_i_u=1.,mass_j_u=4.,charge_i=1.,charge_j=2.,screening_cm=6e-10)
    pair=ion_pair(**state)
    e_si=1.602176634e-19;eps0=8.8541878128e-12;kb_si=1.380649e-23
    mu_si=1.66053906660e-27*4/5
    interaction_si=2*e_si**2/(4*math.pi*eps0)
    energy_si=kb_si*state['temperature'];lambda_si=state['screening_cm']*.01
    g_si=interaction_si/(lambda_si*energy_si)
    omega_si=math.sqrt(2*math.pi/mu_si)*interaction_si**2/energy_si**1.5*reduced_integral(g_si,1,1)
    resistance_si=16/3*(state['number_i']*1e6)*(state['number_j']*1e6)*mu_si*omega_si
    check('SI-cgs collision-integral units',abs(omega_si*1e6/pair['omega_cm3_s']['11']-1),1e-8)
    check('SI-cgs resistance units',abs(resistance_si*.001/pair['resistance_g_cm3_s']-1),1e-8)
    # Independently use the source binary Chapman-Enskog diffusion equation,
    # then recover the resistance from K=n_i*n_j*kT/(n_total*D).
    omega=pair['omega_cm3_s']['11'];nt=state['number_i']+state['number_j']
    d=3*KB*state['temperature']/(16*nt*pair['reduced_mass_g']*omega)
    recovered=state['number_i']*state['number_j']*KB*state['temperature']/(nt*d)
    check('binary diffusion/resistance relation',abs(recovered/pair['resistance_g_cm3_s']-1),1e-12)
    weak=[reduced_integral(1e-100,n,m) for n,m in ((1,1),(1,2),(1,3),(2,2))]
    z,zp,zpp=1-.4*weak[1]/weak[0],2.5-2*weak[1]/weak[0]+.4*weak[2]/weak[0],weak[3]/weak[0]
    for name,value,limit in (('z',z,.6),('zprime',zp,1.3),('zdoubleprime',zpp,2.)):
        check('weak-coupling '+name,abs(value/limit-1),.01)
    screen=screening_length([1e27],[2.],1e7,electron_stiffness_erg=KB*1e7)
    ai=(3/(4*math.pi*1e27))**(1/3)
    lambda_i2=KB*1e7/(4*math.pi*4*E2*1e27)
    lambda_e2=KB*1e7/(4*math.pi*E2*2e27)
    independent=1/math.sqrt(1/lambda_e2+1/(lambda_i2+ai*ai))
    check('pure-ion screening limit',abs(screen['length_cm']/independent-1),1e-12)
    guards=[]
    for name,call in [
        ('electron mass',lambda:ion_pair(**{**state,'mass_i_u':.00055})),
        ('attractive charges',lambda:ion_pair(**{**state,'charge_i':-1.})),
        ('neutral particles',lambda:screening_length([1e20],[0.],1e7)),
        ('invalid screening response',lambda:screening_length([1e20],[1.],1e7,electron_stiffness_erg=-1.))]:
        try:call()
        except ValueError as error:guards.append(dict(name=name,reason=str(error)))
        else:raise ValueError('invalid collision input accepted: '+name)
    for path,sha in inputs.items():
        if digest(Path(path))!=sha:raise ValueError('source changed during comparison')
    passed=all(c['pass_check'] for c in checks)
    report=dict(scope=__doc__,status='pass' if passed else 'fail',accepted_for_stellar_evolution=False,
                cross_sections=cross,collision_integrals=moments,checks=checks,invalid_inputs=guards,
                source_tables={'cross_sections':CROSS,'collision_integrals':{str(k):v for k,v in COLLISION.items()}},
                source_typographical_interpretation='The undefined y in C16 is taken as w; numerical scattering checks test this interpretation.',
                maximum_cross_section_fit_error=max(r['fit_relative_error'] for r in cross),
                maximum_collision_integral_fit_error=max(r['fit_relative_error'] for r in moments),
                maximum_collision_quadrature_change=max(r['quadrature_relative_change'] for r in moments),
                dimensional_control_pair=pair,
                references={'repulsive_ion_integrals':'https://doi.org/10.1103/PhysRevE.93.043203',
                            'burgers_resistance_definition':'https://arxiv.org/abs/1506.03146'},
                input_sha256=inputs)
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ('status','maximum_cross_section_fit_error',
                      'maximum_collision_integral_fit_error','maximum_collision_quadrature_change')}),flush=True)
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
