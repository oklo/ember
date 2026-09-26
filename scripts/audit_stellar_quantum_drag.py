#!/usr/bin/env python3
"""Compare direct quantum and Born electron drag in representative stellar layers.

Both calculations use the same prescribed nonrelativistic, static Yukawa
model and displaced Fermi distribution. The only changed physical
approximation is the treatment of single-ion quantum scattering.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import expit

from electron_ion_born import KB, ME, HBAR, E2, drift_integral
from quantum_yukawa import scatter_many
from audit_quantum_yukawa_v2 import independent_phase
from fetch_tops_composition import digest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--workers',type=int,default=4)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('preserve completed thermal collision comparisons')
    if not 1<=a.workers<=4:raise ValueError('one through four workers allowed')
    paths=[Path(__file__),Path('scripts/quantum_yukawa.py'),Path('scripts/quantum_yukawa_probe_v2.cpp'),
           Path('scripts/electron_ion_born.py'),Path('scripts/audit_quantum_yukawa_v2.py'),a.probe,
           Path('docs/results/quantum_yukawa_v2.json'),Path('docs/results/stellar_electron_drag_3890gyr_v1.json')]
    inputs={str(p.resolve()):digest(p) for p in paths}
    checked=json.loads(paths[-2].read_text())
    if checked['status']!='pass' or checked['input_sha256'][str(a.probe.resolve())]!=digest(a.probe):
        raise ValueError('quantum scattering probe has not passed its controls')
    profile=json.loads(paths[-1].read_text())
    source=profile['records']
    chosen=[('center',source[0]),
            ('outer_nonconvective',max((r for r in source if not r['convective']),key=lambda r:r['helium_Born_parameter'])),
            ('outer_hot',max(source,key=lambda r:r['helium_Born_parameter']))]
    tasks=[(label,row,screen,charge) for label,row in chosen
           for screen in row['modes'][0]['cases'] for charge in (1,2)]

    def run(task):
        label,row,screen,charge=task
        T,eta=row['temperature_K'],row['eta_nonrelativistic'];length=screen['screening_length_cm']
        at=math.sqrt(2*ME*KB*T)*length/HBAR
        c=2*ME*charge*E2*length/HBAR**2
        bounds=[0.,.25,.5,1.,2.,3.,4.,6.,8.]
        if eta>0 and math.sqrt(eta) not in bounds:bounds=sorted(bounds+[math.sqrt(eta)])
        def nodes(order,edges):
            x,w=leggauss(order)
            return [(float((hi+lo)/2+(hi-lo)/2*v),float((hi-lo)/2*weight))
                    for lo,hi in zip(edges[:-1],edges[1:]) for v,weight in zip(x,w)]
        coarse=nodes(16,bounds);fine=nodes(24,bounds);tail=nodes(16,[8.,10.])
        points=coarse+fine+tail
        queries=[(float(at*t),float(c),.0005,32.,int(math.ceil(12*at*t+10))) for t,w in points]
        # Radial and independent ODE controls at the smallest thermal momentum;
        # this extends the earlier scattering audit into the low-energy tail.
        smallest=len(coarse)
        q=queries[smallest]
        queries.append((q[0],q[1],.00025,q[3],q[4]))
        raw=scatter_many(a.probe,queries)
        values=[]
        for (t,w),r in zip(points,raw[:-1],strict=True):
            value=w*t**5*r['cross_section_over_lambda2']*expit(eta-t*t)*expit(t*t-eta)
            values.append(float(value))
        nc,nf=len(coarse),len(fine)
        ic=sum(values[:nc]);iv=sum(values[nc:nc+nf]);it=sum(values[nc+nf:])
        born=math.pi*c*c/(4*at**4)*drift_integral(eta,4*at**2)
        radial_error=abs(raw[-1]['cross_section_over_lambda2']/raw[smallest]['cross_section_over_lambda2']-1)
        phase_errors=[]
        for ell in (0,1):
            exact=independent_phase(q[0],q[1],ell,32.)
            phase_errors.append(abs(math.asin(math.sin(raw[smallest]['phase_shifts'][ell]-exact))))
        quadrature_error=abs(ic/iv-1);tail_error=it/iv
        passed=quadrature_error<1e-3 and tail_error<1e-6 and radial_error<5e-5 and max(phase_errors)<2e-6
        result=dict(layer=label,zone=row['zone'],temperature_K=T,density_g_cm3=row['density_g_cm3'],
                    charge=charge,screening=screen['screening'],screening_length_cm=length,
                    eta_nonrelativistic=eta,thermal_wavenumber=at,potential_strength=c,
                    quantum_over_Born_drag=(iv+it)/born,
                    relative_quadrature_change=quadrature_error,relative_thermal_tail=tail_error,
                    low_momentum_radial_error=radial_error,low_momentum_phase_errors=phase_errors,
                    passed=passed,scattering_queries=len(queries),
                    samples=[dict(thermal_speed=t,quadrature_weight=w,cross_section_over_lambda2=r['cross_section_over_lambda2'])
                             for (t,w),r in zip(points,raw[:-1],strict=True)],
                    sample_counts=dict(coarse=nc,fine=nf,tail=len(tail)))
        print(json.dumps({k:result[k] for k in ['layer','charge','screening','quantum_over_Born_drag','relative_quadrature_change','passed']}),flush=True)
        return result

    with ThreadPoolExecutor(max_workers=a.workers) as pool:results=list(pool.map(run,tasks))
    for path,h in inputs.items():
        if digest(Path(path))!=h:raise ValueError('quantum thermal source changed')
    passed=all(r['passed'] for r in results)
    report=dict(scope=__doc__,status='pass' if passed else 'fail',accepted_for_stellar_evolution=False,
                cases=results,scattering_queries=sum(r['scattering_queries'] for r in results),
                maximum_quadrature_change=max(r['relative_quadrature_change'] for r in results),
                maximum_radial_change=max(r['low_momentum_radial_error'] for r in results),
                maximum_independent_phase_error=max(max(r['low_momentum_phase_errors']) for r in results),
                limitations=['Three saved stellar layers, two prescribed screening choices and H/He charges only.',
                             'Single-ion nonrelativistic Yukawa scattering; no ionic correlations, dynamic screening or electron heat perturbations.',
                             'No metal quantum correction, full stellar velocity, abundance evolution or accepted diffusion prescription.'],
                input_sha256=inputs)
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=report['status'],cases=len(results),scattering_queries=report['scattering_queries'],
                         maximum_quadrature_change=report['maximum_quadrature_change'])),flush=True)
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
