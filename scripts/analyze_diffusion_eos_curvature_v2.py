#!/usr/bin/env python3
"""Reuse the composition experiment, removing known ideal mixing before differences.

The first experiment's central finite difference was too wide for the X ln X
term. Subtract that analytically known term, difference the remaining free
energy, and restore its exact curvature. No source evaluation is repeated.
The original failed comparison and all 70 direct evaluations are retained.
This diagnoses the selected EOS; it does not install a diffusion prescription.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np

from audit_diffusion_eos_curvature import differences, RGAS


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root=Path(__file__).resolve().parents[1]
    work=Path('/tmp/ember-diffusion-eos-curvature-run-v1')
    raw=work/'source/direct_source.json'
    source=json.loads(raw.read_text())
    runtime=np.loadtxt(work/'source/ember_stdout.txt')
    receipt=json.loads((work/'receipt.json').read_text())
    assert receipt['outcome']=='failed' and receipt['jobs'][-1]['exit_code']==1
    for job in receipt['jobs']:
        for p,h in job['input_sha256'].items():assert digest(p)==h
    assert len(source)==70 and runtime.shape==(35,10)
    profile=np.genfromtxt(root/'docs/reports/2026-09-11/evolution_latest_profile.csv',
                          delimiter=',',names=True)
    cases=[]
    for ci,label,zone in [(0,'center',0),(1,'interior',250),
                          (2,'outer_nonconvective',396),
                          (3,'hydrogen_01625',0),(4,'hydrogen_065',0)]:
        center=next(v for v in source if v['build']==1 and v['case']==ci and v['offset']==0)
        adjacent=next(v for v in source if v['build']==1 and v['case']==ci and v['offset']==1)
        x=center['X'];h=adjacent['X']-x
        ideal=RGAS*(1/x+1/(4*(.98-x)))
        ds=[];original=[]
        for bi in range(2):
            subset=[v for v in source if v['build']==bi and v['case']==ci]
            original.append(differences({v['offset']:v['phi'] for v in subset},h))
            residual={}
            for v in subset:
                xx=np.longdouble(v['X']);n4=(np.longdouble('.98')-xx)/4
                mix=np.longdouble(RGAS)*(xx*np.log(xx)+n4*np.log(n4))
                residual[v['offset']]=np.longdouble(v['phi'])-mix
            d=differences(residual,np.longdouble(h))
            ds.append({k:float(v+ideal) if k!='slope' else float(v) for k,v in d.items()})
        ref=ds[1]['richardson_fine']
        step_error=abs(ds[1]['richardson']/ref-1)
        quad_error=abs(ds[0]['richardson_fine']/ref-1)
        dtable=differences(dict(zip([-2,-1,-.5,0,.5,1,2],runtime[ci*7:(ci+1)*7,4])),h)
        assert ref>0 and step_error<1e-4 and quad_error<1e-4
        assert abs(dtable['richardson_fine']/ref)<1e-5
        cases.append({'label':label,'zone':zone,'saved_X':float(profile[zone]['X']),
                      'X':x,'Y3':0.,'Z':.02,'h':h,'T':float(runtime[ci*7,2]),
                      'rho':float(runtime[ci*7,3]),'source_composition_interval':[x-4*h,x+4*h],
                      'source_original_differences':original,'source_analytic_mixing_removed':ds,
                      'source_step_relative_change':step_error,
                      'source_quadrature_relative_change':quad_error,
                      'table_curvature_over_Rgas':float(dtable['richardson_fine']/RGAS),
                      'source_curvature_over_Rgas':ref/RGAS,
                      'table_over_source_curvature':float(dtable['richardson_fine']/ref),
                      'ideal_ion_curvature_over_Rgas':ideal/RGAS,
                      'ideal_ion_only_relative_difference':ideal/ref-1})

    electron_probe=Path('/tmp/ember-ion-collision-integrals-run-v1/electron_probe')
    text=''.join(f'{c["rho"]:.17g} {c["T"]:.17g} {c["X"]:.17g} 0 {(.98-c["X"]):.17g} .02\n'
                 for c in cases)
    p=subprocess.run([str(electron_probe)],input=text,text=True,capture_output=True,timeout=20)
    if p.returncode:raise RuntimeError(p.stderr)
    erows=np.array([list(map(float,l.split())) for l in p.stdout.splitlines()])
    assert erows.shape==(5,4) and np.all(np.isfinite(erows))
    for c,(ne,pe,dpe,dpet) in zip(cases,erows,strict=True):
        # At fixed rho, replacing He4 by H1 increases ne by rho*NA/2 per X.
        derivative=.5*c['rho']*6.02214076e23
        chi_e=dpe/(c['rho']*c['T']*RGAS)*(derivative/ne)**2
        combined=c['ideal_ion_curvature_over_Rgas']+chi_e
        c.update(ideal_electron_curvature_over_Rgas=float(chi_e),
                 ideal_ions_electrons_relative_difference=float(combined/c['source_curvature_over_Rgas']-1))
    builds=[]
    for name in ['ember-freeeos-precision-builder-check-v1',
                 'ember-freeeos-quadrature13-builder-check-v1']:
        directory=Path('/tmp')/name
        b=json.loads((directory/'build_receipt.json').read_text())
        assert digest(directory/'probe')==b['probe_sha256']
        assert digest(b['library'])==b['library_sha256']
        builds.append({'directory':str(directory),'receipt':b,
                       'receipt_sha256':digest(directory/'build_receipt.json')})
    family=root/'data/eos/numerical_electron_base_v1/freeeos300_gs98_z020.dat'
    import shlex
    planes=[family.parent/shlex.split(line)[0] for line in family.read_text().splitlines()[3:] if line]
    result={'created_utc':datetime.now(timezone.utc).isoformat(),'scope':__doc__,
            'outcome':'completed_diagnostic','selected_EOS_accepted_for_diffusion_forces':False,
            'new_stellar_evolution':False,'controls':cases,'reused_direct_source_queries':70,
            'new_direct_FreeEOS_queries':0,'electron_pressure_queries':5,
            'original_failure':{'receipt':str(work/'receipt.json'),'sha256':digest(work/'receipt.json'),
                                'reason':'Central uncorrected finite-difference spacing failed the 1e-4 convergence criterion; no source evaluation failed. Criterion unchanged.'},
            'criteria':{'source_step_relative_change':1e-4,'source_quadrature_relative_change':1e-4,
                        'maximum_absolute_table_over_source_curvature':1e-5},
            'interpretation':'The selected composition interpolation lacks within-cell H1/He4 free-energy curvature. Restoring ideal ion mixing alone is inadequate at moderate/high X in dense matter. The ideal-ion-plus-electron comparison isolates a residual that still contains nonideal and ionization physics; it is not a complete correction or accepted physical model.',
            'electron_comparison_assumptions':'Fully stripped GS98 including K, representative baryonic nuclei, ideal relativistic Fermi electrons. FreeEOS omits K and treats ionization and nonideal interactions; these source differences are not fitted away.',
            'electron_probe_input':text,'electron_probe_output':p.stdout,
            'source_builds':builds,'EOS_inputs_sha256':{str(p):digest(p) for p in [family,*planes]},
            'inputs_sha256':{str(p):digest(p) for p in [Path(__file__),raw,work/'source/ember_stdout.txt',
                root/'scripts/audit_diffusion_eos_curvature.py',root/'scripts/diffusion_eos_potential_probe.cpp',
                root/'src/eos_mixture.cpp',electron_probe,root/'scripts/diffusion_electron_probe.cpp']}}
    output=root/'docs/results/diffusion_eos_composition_curvature_v2.json'
    with output.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'outcome':result['outcome'],'cases':len(cases),
                      'max_step_error':max(c['source_step_relative_change'] for c in cases),
                      'max_quadrature_error':max(c['source_quadrature_relative_change'] for c in cases),
                      'max_ion_only_error':max(abs(c['ideal_ion_only_relative_difference']) for c in cases),
                      'max_ions_electrons_error':max(abs(c['ideal_ions_electrons_relative_difference']) for c in cases)}))


if __name__=='__main__':main()
