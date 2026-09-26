#!/usr/bin/env python3
"""Check energy-mode convergence while preserving the three-mode reference."""
import concurrent.futures
from datetime import datetime,timezone
import hashlib,json,math,time
from pathlib import Path
import numpy as np
from electron_electron_energy_modes import PolynomialBasis,energy_mode_matrix
from electron_electron_collision import energy_basis,born_electron_ion_moments
from electron_ion_heat import born_response


def execute(task):
    name,args=task;start=time.monotonic();result=energy_mode_matrix(**args)
    result['elapsed_seconds']=time.monotonic()-start
    return name,result


def main():
    work=Path('/tmp/ember-electron-electron-energy-modes-run-v1')
    report=Path('docs/results/electron_electron_energy_modes_v1.json')
    assert work.is_dir() and not report.exists()
    source=Path('docs/results/electron_electron_collision_v1.json');old=json.loads(source.read_text())
    assert old['outcome']=='passed'
    checks=[];rows={};tasks=[]
    def check(name,error,tolerance,**meta):
        checks.append(dict(name=name,error=float(error),tolerance=tolerance,
                           passed=bool(np.isfinite(error) and abs(error)<=tolerance),**meta))
    for r in old['physical_comparisons']:
        eta=r['eta'];b=r['b_thermal'];zone=r['zone']
        q=PolynomialBasis(eta,7);basis=energy_basis(eta)
        x=np.linspace(0,max(eta,0)+35,100)
        compare=q.values(x)[:3]-basis.polynomials(x)
        check('independent low-degree basis',np.max(abs(compare))/(1+np.max(abs(basis.polynomials(x)))),2e-10,zone=zone)
        check('current and higher-mode orthogonality',q.gram_error,2e-9,zone=zone)
        coarse=q.ion_matrix(b,points=256);fine=q.ion_matrix(b,points=384)
        scales=np.sqrt(np.diag(fine))
        check('ion energy quadrature',np.max(abs(coarse-fine)/scales[:,None]/scales[None,:]),2e-6,zone=zone)
        expected=born_electron_ion_moments(basis,b)
        check('retained electron-ion collision matrix',np.max(abs(fine[:3,:3]-expected))/np.max(abs(expected)),2e-7,zone=zone)
        tasks.append((f'zone_{zone}',dict(eta=eta,b_thermal=b,degree=7)))
    # A dilute weak-screening control tests the same higher-mode procedure
    # toward the classical plasma limit, without assuming a Spitzer value.
    tasks.append(('dilute_weak_screening',dict(eta=-16.,b_thermal=1e6,degree=7)))
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
        futures=[pool.submit(execute,t) for t in tasks]
        for future in concurrent.futures.as_completed(futures):
            name,result=future.result();rows[name]=result
            (work/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n')
            print(name,'completed',flush=True)
    summaries=[]
    for name,r in rows.items():
        ee=np.asarray(r['collision_matrix']);ei=np.asarray(r['electron_ion_matrix'])
        if name.startswith('zone_'):
            prior=np.asarray(old['calculations'][name+'_high']['collision_matrix'])
            scale=np.sqrt(np.diag(prior)[1:])
            check('independent four-momentum versus reduced three-mode matrix',
                  np.max(abs(ee[1:3,1:3]-prior[1:,1:])/scale[:,None]/scale[None,:]),2e-9,case=name)
        # Fix U0=0 before minimizing over all zero-current higher modes.
        histories=[]
        exact=born_response(r['eta'],r['b_thermal'])['zero_drift_heat_response_exact']
        for modes in range(2,len(ei)+1):
            elastic=ei[1:modes,1:modes];combined=(ei+ee)[1:modes,1:modes]
            h0=float(np.linalg.inv(elastic)[0,0]);h1=float(np.linalg.inv(combined)[0,0])
            histories.append(dict(modes=modes,without_ee=h0,with_ee=h1,with_ee_over_full_elastic=h1/exact))
            check('collisions reduce heat response',max(0.,h1/h0-1),2e-10,case=name,modes=modes)
            check('variational Lorentz upper bound',max(0.,h0/exact-1),2e-6,case=name,modes=modes)
            if len(histories)>1:
                check('heat response increases with trial space',max(0.,histories[-2]['with_ee']/h1-1),2e-9,case=name,modes=modes)
        change=histories[-1]['with_ee']/histories[-2]['with_ee']-1
        check('last added energy mode',change,.005,case=name)
        summaries.append(dict(case=name,eta=r['eta'],b_thermal=r['b_thermal'],
                              heat_convergence=histories,last_mode_relative_change=change,
                              assumption='pure-H stationary scatterers; saved degeneracy and prescribed screening'))
    files=[Path(__file__),Path('scripts/electron_electron_energy_modes.py'),
           Path('scripts/electron_electron_collision.py'),source]
    output=dict(created_utc=datetime.now(timezone.utc).isoformat(),
        scope='Born/Pauli electron heat energy-mode convergence; no stellar mixture diffusion',
        outcome='passed' if all(c['passed'] for c in checks) else 'failed',
        accepted_for_stellar_evolution=False,checks=checks,calculations=rows,summaries=summaries,
        inputs_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        limitations=['Higher-mode quadrature refinement remains required.',
                     'A small last-mode change alone does not prove convergence to arbitrary distributions.',
                     'No dynamic screening, correlated mixture or full quantum scattering.',
                     'No stellar transport update or ultra-cold validation.'])
    report.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(dict(outcome=output['outcome'],checks=len(checks),failed=[c for c in checks if not c['passed']],summaries=summaries)),flush=True)
    if output['outcome']!='passed':raise SystemExit(1)


if __name__=='__main__':main()
