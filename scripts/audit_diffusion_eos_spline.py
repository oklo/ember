#!/usr/bin/env python3
"""Test smooth hydrogen interpolation after removing analytic ion mixing.

This is an offline Y3=0 comparison at three fixed stellar temperatures and
densities. It reuses the 36 existing composition planes and the independent
FreeEOS evaluations from the curvature diagnostic. It neither changes the
production EOS nor validates interpolation in helium-3 or selective metals.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import xlogy

from audit_diffusion_eos_curvature import RGAS


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root=Path(__file__).resolve().parents[1]
    reference=root/'docs/results/diffusion_eos_composition_curvature_v2.json'
    ref=json.loads(reference.read_text())
    for p,h in ref['EOS_inputs_sha256'].items():assert digest(p)==h
    raw=Path('/tmp/ember-diffusion-eos-curvature-run-v1/source/direct_source.json')
    assert digest(raw)==ref['inputs_sha256'][str(raw)]
    source=json.loads(raw.read_text())
    family=root/'data/eos/numerical_electron_base_v1/freeeos300_gs98_z020.dat'
    axis=np.array([float(v) for v in family.read_text().splitlines()[1].split()[2:]])
    coords=list(dict.fromkeys((c['T'],c['rho']) for c in ref['controls']))
    queries=[(t,rho,x) for t,rho in coords for x in axis]
    text=''.join(f'{x:.17g} 0 {t:.17g} {rho:.17g}\n' for t,rho,x in queries)
    probe=Path('/tmp/ember-diffusion-eos-curvature-run-v1/probe')
    p=subprocess.run([str(probe),str(family)],input=text,text=True,capture_output=True,timeout=90)
    work=Path('/tmp/ember-diffusion-eos-spline-run-v1')
    (work/'queries.txt').write_text(text);(work/'stdout.txt').write_text(p.stdout)
    (work/'stderr.txt').write_text(p.stderr)
    if p.returncode:raise RuntimeError(p.stderr)
    vals=np.array([list(map(float,l.split())) for l in p.stdout.splitlines()])
    assert vals.shape==(len(queries),10) and np.all(np.isfinite(vals))
    mixing=RGAS*(xlogy(axis,axis)+xlogy((.98-axis)/4,(.98-axis)/4))
    splines={}
    for i,c in enumerate(coords):
        data=vals[i*len(axis):(i+1)*len(axis)]
        residual=data[:,4]-mixing
        # This spline is linear in source values. The same composition
        # weights act on all T/rho derivatives of the shared potential.
        splines[c]=(CubicSpline(axis,residual,bc_type='not-a-knot',extrapolate=False),
                    CubicSpline(axis,data[:,5:7],bc_type='not-a-knot',extrapolate=False))
    comparisons=[]
    arad=4*5.670374419e-5/2.99792458e10
    for ci,c in enumerate(ref['controls']):
        x,t,rho=c['X'],c['T'],c['rho'];f,pe=splines[t,rho]
        curvature=float(f(x,2)/RGAS+1/x+1/(4*(.98-x)))
        actual=c['source_curvature_over_Rgas']
        direct=next(v for v in source if v['build']==1 and v['case']==ci and v['offset']==0)
        prad=arad*t**4/3
        targets=np.array([direct['P']+prad,direct['E']+3*prad/rho])
        changes=pe(x)/targets-1
        discrepancy=curvature/actual-1
        comparisons.append({'case':c['label'],'X':x,'T':t,'rho':rho,
                            'curvature_over_Rgas':curvature,'source_curvature_over_Rgas':actual,
                            'relative_curvature_difference':discrepancy,
                            'relative_pressure_difference':float(changes[0]),
                            'relative_energy_difference':float(changes[1]),
                            'within_comparison_criteria':bool(abs(discrepancy)<.005 and np.max(np.abs(changes))<.003)})
    result={'created_utc':datetime.now(timezone.utc).isoformat(),'scope':__doc__,
            'outcome':'completed_interpolation_pilot','comparisons':comparisons,
            'comparison_criteria':{'curvature':.005,'P_and_E':.003},
            'all_comparisons_pass':all(c['within_comparison_criteria'] for c in comparisons),
            'accepted_for_stellar_evolution':False,'new_stellar_models':0,
            'new_FreeEOS_source_queries':0,'existing_EOS_queries':len(queries),
            'new_scheme':'Interpolate Phi=F/T minus Rgas times [X ln X+(Y4/4) ln(Y4/4)] with a C2 not-a-knot cubic spline; restore the analytic mixing term. No density, temperature or composition extrapolation.',
            'remaining_work':'Verify full T/rho responses, both composition directions, boundaries and all source coverage; construct consistent chemical forces and thermal diffusion and couple abundance/energy fluxes. The pilot does not choose a metal-separation model.',
            'inputs_sha256':{str(p):digest(p) for p in [Path(__file__),reference,raw,probe,family]},
            'output_sha256':{str(p):digest(p) for p in [work/'queries.txt',work/'stdout.txt',work/'stderr.txt']}}
    out=root/'docs/results/diffusion_eos_hydrogen_spline_v1.json'
    with out.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'all_comparisons_pass':result['all_comparisons_pass'],
                      'max_curvature_error':max(abs(c['relative_curvature_difference']) for c in comparisons),
                      'max_pressure_error':max(abs(c['relative_pressure_difference']) for c in comparisons),
                      'max_energy_error':max(abs(c['relative_energy_difference']) for c in comparisons)}))


if __name__=='__main__':main()
