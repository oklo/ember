#!/usr/bin/env python3
"""Build a small table of directly integrated repulsive-ion collision moments.

The table covers a stated collision-strength interval, with independent
midpoints and quadrature controls. It is not an electron collision table or
an acceptance of a complete stellar diffusion prescription.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

from fetch_tops_composition import digest
from ion_collision_integrals import reduced_integral
from yukawa_scattering_control_v2 import collision_integrals

PAIRS=((1,1),(1,2),(1,3),(2,2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--report',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists() or a.report.exists():raise FileExistsError('preserve completed collision tables')
    paths=[Path(__file__),Path('scripts/yukawa_scattering_control.py'),Path('scripts/yukawa_scattering_control_v2.py'),
           Path('scripts/ion_collision_integrals.py'),Path('docs/results/yukawa_orbit_controls_v1.json'),
           Path('docs/results/ion_collision_integrals_v1.json')]
    inputs={str(p.resolve()):digest(p) for p in paths}
    def values(g,order=32,angle=128,radius=256,tail=8.):
        r=collision_integrals(g,thermal_order=order,angle_order=angle,radius_order=radius,thermal_max=tail)
        return np.array([r[k] for k in PAIRS])
    numerical=[]
    for g in (1e-4,1e-2,.05,.3,1.,3.,30.,100.):
        coarse=values(g,32);fine=values(g,64,192,384);tail=values(g,32,128,256,10.)
        numerical.append(dict(strength=g,coarse=coarse.tolist(),fine=fine.tolist(),
                              relative_quadrature_change=(coarse/fine-1).tolist(),
                              relative_tail_change=(tail/coarse-1).tolist(),
                              printed_fit_relative_difference=[reduced_integral(g,*p)/fine[i]-1 for i,p in enumerate(PAIRS)]))
        print(json.dumps({'quadrature_strength_complete':g}),flush=True)
    x=np.linspace(-4,2,97)
    table=np.array([values(10**v) for v in x])
    interpolator=CubicSpline(x,np.log(table),axis=0,extrapolate=False)
    comparisons=[]
    for xm in (x[:-1]+x[1:])/2:
        reference=values(10**xm)
        predicted=np.exp(interpolator(xm))
        comparisons.append(dict(log10_strength=float(xm),reference=reference.tolist(),
                                interpolated=predicted.tolist(),relative_error=(predicted/reference-1).tolist()))
    numerical_error=max(abs(e) for r in numerical for e in r['relative_quadrature_change'])
    tail_error=max(abs(e) for r in numerical for e in r['relative_tail_change'])
    interpolation_error=max(abs(e) for r in comparisons for e in r['relative_error'])
    passed=numerical_error<2e-5 and tail_error<2e-5 and interpolation_error<1e-3
    # Explicit bounds: callers must not clamp unsupported coupling values.
    bounds=[float(x[0]),float(x[-1])]
    if not np.all(np.isnan(interpolator([-4.01,2.01]))):raise ValueError('unexpected interpolation extrapolation')
    content=dict(scope=__doc__,collision_strength_min=1e-4,collision_strength_max=100.,
                 moment_pairs=PAIRS,log10_strength=x.tolist(),dimensionless_integrals=table.tolist(),
                 interpolation='natural variable log10(g), cubic not-a-knot interpolation of ln(K_nm), no extrapolation',
                 numerical_status='pass' if passed else 'fail',accepted_for_stellar_evolution=False,
                 input_sha256=inputs)
    a.output.write_text(json.dumps(content,indent=2,allow_nan=False)+'\n')
    for p,h in inputs.items():
        if digest(Path(p))!=h:raise ValueError('collision input changed during integration')
    report=dict(scope=__doc__,status='pass' if passed else 'fail',accepted_for_stellar_evolution=False,
                table=str(a.output),table_sha256=digest(a.output),nodes=len(x),midpoint_controls=len(comparisons),
                maximum_quadrature_change=numerical_error,maximum_tail_change=tail_error,
                maximum_interpolation_error=interpolation_error,
                quadrature_controls=numerical,interpolation_controls=comparisons,
                limitations=['Repulsive classical Yukawa scattering only; no attractive or quantum electron-ion treatment.',
                             'Screening prescription and ionization are external physical inputs.',
                             'A numerical cross-section table does not validate the effective potential at arbitrary plasma coupling.',
                             'No stellar diffusion velocity or abundance evolution has been calculated.'],
                input_sha256=inputs)
    a.report.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ('status','nodes','midpoint_controls','maximum_quadrature_change',
                      'maximum_tail_change','maximum_interpolation_error')}),flush=True)
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
