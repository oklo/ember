#!/usr/bin/env python3
"""Refine full pair quadrature, electron energy limits and trial-space size."""
import concurrent.futures
from datetime import datetime,timezone
import hashlib,json,math,time
from pathlib import Path
import numpy as np
from scipy.special import expit,roots_jacobi
from electron_electron_energy_modes import PolynomialBasis as PreviousBasis
from electron_electron_energy_modes_v2 import PolynomialBasis,energy_mode_matrix
from electron_ion_heat import born_response


def execute(task):
    name,args=task;start=time.monotonic();r=energy_mode_matrix(**args)
    r['elapsed_seconds']=time.monotonic()-start
    return name,r


def gram(basis,tail=100.):
    eta=basis.base.eta;length=max(eta,0)+tail
    u,w=roots_jacobi(640,0,1.5);x=(u+1)*length/2
    w=w*(length/2)**2.5*expit(eta-x)*expit(x-eta)/basis.base.normalization
    q=basis.values(x);q[1]*=math.sqrt(basis.base.variance)
    return (q*w)@q.T


def response(m):
    # Higher energy amplitudes carry neither of the retained currents.
    inverse=np.linalg.inv(m)
    return inverse[:2,:2]


def main():
    work=Path('/tmp/ember-electron-energy-refinement-run-v1')
    target=Path('docs/results/electron_electron_energy_refinement_v1.json')
    assert work.is_dir() and not target.exists()
    path=Path('docs/results/electron_electron_energy_modes_v1.json')
    previous=json.loads(path.read_text());assert previous['outcome']=='passed'
    checks=[];tasks=[];records={};diagnostics=[]
    def check(name,error,tolerance,**meta):
        checks.append(dict(name=name,error=float(error),tolerance=tolerance,
                           passed=bool(np.isfinite(error) and abs(error)<=tolerance),**meta))
    for name,r in previous['calculations'].items():
        eta,b=r['eta'],r['b_thermal'];basis=PolynomialBasis(eta,9)
        current=gram(basis);prior_gram=gram(PreviousBasis(eta,7))
        check('independent extended-tail polynomial norm',np.max(abs(current-np.eye(10))),2e-9,case=name)
        coarse=basis.ion_matrix(b);fine=basis.ion_matrix(b,points=768,tail=100.)
        scale=np.sqrt(np.diag(fine))
        check('momentum quadrature and electron tail',np.max(abs(coarse-fine)/scale[:,None]/scale[None,:]),2e-7,case=name)
        two=np.asarray(born_response(eta,b)['collision_moments'])
        scale=np.sqrt(np.diag(two))
        check('independent electron-ion lower block',np.max(abs(fine[:2,:2]-two)/scale[:,None]/scale[None,:]),2e-8,case=name)
        diagnostics.append(dict(case=name,old_eight_mode_norm_error_on_extended_domain=float(np.max(abs(prior_gram-np.eye(8)))),
                                new_ten_mode_norm_error_on_extended_domain=float(np.max(abs(current-np.eye(10))))))
        tasks.append((name,dict(eta=eta,b_thermal=b,degree=9,orders=(80,40,40,48,48),tail=50.,normalization_tail=80.)))
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
        for future in concurrent.futures.as_completed([pool.submit(execute,t) for t in tasks]):
            name,r=future.result();records[name]=r
            (work/(name+'.json')).write_text(json.dumps(r,indent=2)+'\n')
            print(name,'completed',flush=True)
    summaries=[]
    for name,r in records.items():
        old=previous['calculations'][name]
        ee,ei=np.asarray(r['collision_matrix']),np.asarray(r['electron_ion_matrix'])
        oe,oi=np.asarray(old['collision_matrix']),np.asarray(old['electron_ion_matrix'])
        # Compare the same eight-dimensional trial space after changing all
        # integration limits and orders. The effective two-current response
        # is insensitive to normalization of the eliminated polynomials.
        histories=[]
        for strength in [0.,.5,1.]:
            oldr=response(oi+strength*oe);newr=response((ei+strength*ee)[:8,:8])
            scale=np.sqrt(np.diag(newr));err=np.max(abs(newr-oldr)/scale[:,None]/scale[None,:])
            check('refined pair integration at fixed trial space',err,.0002,case=name,ee_prefactor_ratio=strength)
        scale=np.sqrt(np.diag(ee)[1:3])
        check('preserved normalized three-mode collision block',np.max(abs(ee[1:3,1:3]-oe[1:3,1:3])/scale[:,None]/scale[None,:]),2e-6,case=name)
        for count in range(2,11):
            combined=(ei+ee)[:count,:count]
            l=response(combined)
            heat=l[1,1]-l[0,1]**2/l[0,0]
            check('constrained heat response',heat/np.linalg.inv(combined[1:,1:])[0,0]-1,2e-10,case=name,modes=count)
            histories.append(dict(modes=count,heat=float(heat),two_current_response=l.tolist()))
        change=histories[-1]['heat']/histories[-3]['heat']-1
        check('two additional energy modes',change,.005,case=name)
        last=histories[-1]['heat']/histories[-2]['heat']-1
        exact=born_response(r['eta'],r['b_thermal'])['zero_drift_heat_response_exact']
        summaries.append(dict(case=name,heat_response_with_ee_over_exact_elastic=histories[-1]['heat']/exact,
                              eight_to_ten_mode_relative_change=change,last_mode_relative_change=last,convergence=histories))
    files=[Path(__file__),Path('scripts/electron_electron_energy_modes.py'),
           Path('scripts/electron_electron_energy_modes_v2.py'),Path('scripts/electron_electron_collision.py'),path]
    r=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed' if all(c['passed'] for c in checks) else 'failed',
           scope='Extended energy intervals and ten-mode Born/Pauli comparison; no selected stellar coefficients',
           accepted_for_stellar_evolution=False,checks=checks,polynomial_tail_diagnostics=diagnostics,
           calculations=records,summaries=summaries,reused_report=str(path),reused_checks=len(previous['checks']),
           inputs_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
           limitations=['Prescribed static Born potential; no correlated or relativistic collision physics.',
                        'These four sampled states do not establish accuracy throughout the full stellar track.',
                        'Strong degeneracy, nonideal heat conventions and conservative stellar evolution remain unresolved.'])
    target.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(dict(outcome=r['outcome'],checks=len(checks),failed=[c for c in checks if not c['passed']],polynomial_tail_diagnostics=diagnostics,
                         summaries=[{k:v for k,v in s.items() if k!='convergence'} for s in summaries])),flush=True)
    if r['outcome']!='passed':raise SystemExit(1)


if __name__=='__main__':main()
