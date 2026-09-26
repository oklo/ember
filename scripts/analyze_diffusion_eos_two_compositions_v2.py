#!/usr/bin/env python3
"""Correct the isotope entropy convention, reusing every completed source query.

Version 1 used 3/4 in a translational-mass logarithm when the production EOS
uses the atomic masses 3.01602932/4.00260325. Mixing source conventions created
spurious curvature in the three-plane comparison. Baryonic population counts
remain X/A with integer A. No stellar EOS or source output is altered.
"""
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path

import numpy as np

from composition_potential_v2 import (RGAS,PotentialSpline,mixing_value,
                                      source_mixing_value,mixing_derivatives,
                                      ISOTOPE_MASS_RATIO)
from audit_diffusion_eos_two_compositions import normalized_error


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root=Path(__file__).resolve().parents[1]
    original=root/'docs/results/diffusion_eos_two_compositions_v1.json'
    old=json.loads(original.read_text())
    for field in ['inputs_sha256','output_sha256']:
        for p,h in old[field].items():assert digest(p)==h
    work=Path('/tmp/ember-diffusion-eos-two-compositions-run-v1')
    sources={tuple(v['key']):v for v in map(json.loads,(work/'direct_source.jsonl').read_text().splitlines())}
    cached={tuple(row[:4]):row for name in [
        '/tmp/ember-diffusion-eos-spline-run-v1/stdout.txt',str(work/'runtime_stdout.txt')]
            for row in np.loadtxt(name)}
    family=root/'data/eos/numerical_electron_base_v1/freeeos300_gs98_z020.dat'
    axis=np.array([float(v) for v in family.read_text().splitlines()[1].split()[2:]])
    material=old['material_coordinates'];splines={}
    for mi,m in enumerate(material):
        ys=np.array([0.,.06,.12]);r=np.empty((len(axis),3))
        for i,x in enumerate(axis):
            for j,y in enumerate(ys):
                if y==.06:
                    prad=4*5.670374419e-5/2.99792458e10*m['T']**4/3
                    r[i,j]=(sources[0,mi,x,y]['phi']-source_mixing_value(x,y)
                            -prad/(m['rho']*m['T']))
                else:r[i,j]=cached[x,y,m['T'],m['rho']][4]-mixing_value(x,y)
        for label,indices in [('two_planes',[0,2]),('three_planes',[0,1,2])]:
            splines[mi,label]=PotentialSpline(axis,ys[indices],r[:,indices])
    def source_jet(t,bi):
        x,y,hx,hy=t['X'],t['Y3'],t['hx'],t['hy']
        def f(dx,dy):
            xx=x+dx*hx;yy=y+dy*hy
            return np.longdouble(sources[bi,t['material'],xx,yy]['phi'])-np.longdouble(source_mixing_value(xx,yy))
        def difference(s):
            middle=f(0,0)
            g=np.array([(f(s,0)-f(-s,0))/(2*s*hx),(f(0,s)-f(0,-s))/(2*s*hy)],dtype=float)
            xx=((f(s,0)-middle)+(f(-s,0)-middle))/(s*hx)**2
            yy=((f(0,s)-middle)+(f(0,-s)-middle))/(s*hy)**2
            xy=(f(s,s)-f(s,-s)-f(-s,s)+f(-s,-s))/(4*s*s*hx*hy)
            return g,np.array([[xx,xy],[xy,yy]],dtype=float)
        g0,h0=difference(1);g1,h1=difference(.5)
        gi,hi=mixing_derivatives(x,y)
        return gi+(4*g1-g0)/3,hi+(4*h1-h0)/3,hi+h1
    comparisons=[]
    for t in old['comparisons']:
        g0,h0,_=source_jet(t,0);g1,h1,fine=source_jet(t,1)
        precision=normalized_error(h0-h1,h1);spacing=normalized_error(fine-h1,h1)
        assert precision<1e-4 and spacing<1e-4
        result={k:t[k] for k in ['material','material_name','X','Y3','hx','hy']}
        result.update(source_gradient_over_Rgas=(g1/RGAS).tolist(),
                      source_Hessian_over_Rgas=(h1/RGAS).tolist(),
                      source_precision_error=precision,source_spacing_error=spacing,methods={})
        for name in ['two_planes','three_planes']:
            v,g,h=splines[t['material'],name].jet(t['X'],t['Y3'])
            err=normalized_error(h-h1,h1)
            result['methods'][name]={'gradient_over_Rgas':(g/RGAS).tolist(),
                'Hessian_over_Rgas':(h/RGAS).tolist(),'normalized_Hessian_error':err,
                'maximum_gradient_error_over_Rgas':float(np.max(np.abs(g-g1))/RGAS),
                'passes_0p5_percent_Hessian_criterion':err<.005}
        comparisons.append(result)
    result={'created_utc':datetime.now(timezone.utc).isoformat(),'scope':__doc__,
            'outcome':'completed_corrected_composition_derivative_pilot',
            'accepted_for_evolution':False,'material_coordinates':material,
            'comparisons':comparisons,'new_source_or_runtime_queries':0,
            'reused_direct_source_queries':len(sources),'reused_runtime_queries':old['new_runtime_queries']+old['reused_runtime_queries'],
            'isotope_translational_mass_ratio':ISOTOPE_MASS_RATIO,
            'source_builds':old['source_builds'],'criterion':old['criterion'],
            'supersedes':{'report':str(original),'sha256':digest(original),
                          'reason':'Incorrect isotope entropy mass convention in the offline comparison. Original data and report are preserved; the three-plane discrepancy was mostly artificial.'},
            'remaining_work':old['remaining_work'],
            'inputs_sha256':{str(p):digest(p) for p in [Path(__file__),original,
                root/'scripts/composition_potential_v2.py',root/'include/ember/composition.hpp',
                root/'src/eos_mixture.cpp',work/'direct_source.jsonl',work/'runtime_stdout.txt']}}
    out=root/'docs/results/diffusion_eos_two_compositions_v2.json'
    with out.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({name:max(c['methods'][name]['normalized_Hessian_error'] for c in comparisons)
                     for name in ['two_planes','three_planes']}))


if __name__=='__main__':main()
