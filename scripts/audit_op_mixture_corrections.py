#!/usr/bin/env python3
"""Check OP mixture corrections and apply them to retained native spectra.

The archived Fortran scattering routines are compiled in isolation as a
numerical reference after source inspection. They have no data-file access.
No downloaded executable, installation script, EOS or stellar code is run.
"""
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.integrate import quad

from op_native_spectra_v2 import read_mesh, read_spectra, native_rosseland
from compare_op_native_mixtures import interpolation
from op_mixture_corrections import (corrections, transport_shape,
    fermi_integrals, screened_momentum_difference, freefree_subtraction,
    scattering_factor)

ROOT = Path('/Users/greglaughlin/Projects/ember')
WORK = Path('/tmp/ember-op-mixture-corrections-run-v1')
DATA = Path('/tmp/ember-op-selected-spectra-v1/OP4STARS_1.3/mono')
OPCODE = Path('/tmp/ember-op-monochromatic-fetch-v2/inspection/OP4STARS_1.3/opserver/opac_ev.f')
OUTPUT = ROOT/'docs/results/op_mixture_corrections_v1.json'
BASE = ROOT/'docs/results/op_native_mixture_comparison_v1.json'
A0 = 5.29177210903e-9
MU = 1.66053906660e-24


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reference_program():
    source = OPCODE.read_text()
    first = source.index('      SUBROUTINE BRCKR(')
    last = source.index('      END SUBROUTINE FDF1F2', first)+len('      END SUBROUTINE FDF1F2')
    body = source[first:last]+'\n'
    driver = '''      program reference
      implicit none
      real t, ne, ions(28), u(0:100), f(0:100)
      integer n, i
      read(*,*) t, ne
      read(*,*) ions
      read(*,*) n
      read(*,*) (u(i), i=0,n)
      call brckr(t,ne,ions,28,u,n,f)
      do i=0,n
        write(*,'(ES26.17)') f(i)
      enddo
      end program reference
'''
    path = WORK/'op_scattering_reference.f'
    path.write_text(driver+body)
    compiler = shutil.which('gfortran')
    if compiler is None:
        raise ValueError('gfortran is required for the independent source reference')
    executable = WORK/'op_scattering_reference'
    command = [compiler, '-O2', '-fdefault-real-8', '-ffixed-line-length-none',
               str(path), '-o', str(executable)]
    result = subprocess.run(command, capture_output=True, text=True, timeout=30)
    (WORK/'compile.log').write_text(result.stdout+result.stderr)
    result.check_returncode()
    return executable, dict(source=str(OPCODE), source_sha256=digest(OPCODE),
                            extracted_source_sha256=digest(path), command=command,
                            executable_sha256=digest(executable))


def legacy_factor(executable, t, ne, ions, u):
    text = f'{t:.17g} {ne:.17g}\n'+' '.join(f'{v:.17g}' for v in ions)+'\n'
    text += str(len(u)-1)+'\n'+' '.join(f'{v:.17g}' for v in u)+'\n'
    p = subprocess.run([str(executable)], input=text, text=True,
                       capture_output=True, timeout=5)
    p.check_returncode()
    values = np.array([float(s) for s in p.stdout.split()])
    if len(values) != len(u) or not np.isfinite(values).all():
        raise ValueError('unexpected OP reference output')
    return values


def controls(executable):
    results = []
    def check(name, value, bound):
        results.append(dict(name=name, value=float(value), bound=bound,
                            passed=bool(abs(value) <= bound)))
    for d in np.logspace(-8, 10, 19):
        direct = quad(lambda m:(1+m*m)*(1-m)*d/(d+1-m), -1., 1.,
                      epsabs=1e-20, epsrel=2e-11, limit=200)[0]
        check(f'angular_integral_D{d:g}', float(transport_shape(d))/direct-1, 2e-10)
    for eta in (-20., -5., 0., 5., 20.):
        ph, mh = fermi_integrals(eta)
        h = 1e-4
        derivative = (fermi_integrals(eta+h)[0]-fermi_integrals(eta-h)[0])/(2*h)
        check(f'fermi_derivative_eta{eta:g}', derivative/mh-1, 3e-9)
    for alpha in (.0001, .1, 10.):
        for e, photon in ((.1,.01), (1.,1.), (10.,.1), (.001,10.)):
            qplus = math.sqrt(e+photon)+math.sqrt(e)
            qminus = photon/qplus
            direct = quad(lambda logq:(math.exp(logq)/(math.exp(logq)+alpha))**2-1,
                          math.log(qminus*qminus), math.log(qplus*qplus),
                          epsabs=1e-14, epsrel=1e-11)[0]
            value = float(screened_momentum_difference(alpha,e,photon))
            check(f'Born_momentum_{alpha:g}_{e:g}_{photon:g}', value/direct-1, 2e-9)
    ions = np.zeros(28); ions[0] = .3; ions[1] = .7
    u = np.array([.1,.3,1.,4.,10.,20.])
    ne = 1e21; t = 6e4; epa = 1.7
    full = corrections(t,ne,ions,epa,u,64)
    scaled = corrections(t,ne,3*ions,3*epa,u,64)
    for k in ('scattering','scattering_change','freefree_change'):
        check('reference_normalization_'+k, np.max(np.abs(scaled[k]/(3*full[k])-1)), 3e-13)
    # Independent adaptive thermal integration, with the momentum difference
    # itself integrated over log(q^2), at selected photons and ionic charges.
    rydt = t/157894.; aune = 1.48185e-25*ne
    alpha = 5.8804e-19*ne*(ions @ np.arange(1.,29.)**2)/(epa*t)
    for photon_u in (.1,1.,4.,20.):
        w = photon_u*rydt
        answer = 0.
        for z in (1,2):
            def integrand(x):
                if x == 0:return 0.
                e = x*rydt
                lo = math.sqrt(e); hi = math.sqrt(e+w)
                qplus = hi+lo; qminus = w/qplus
                born = quad(lambda s:(math.exp(s)/(math.exp(s)+alpha))**2-1,
                            2*math.log(qminus),2*math.log(qplus),epsabs=2e-13,epsrel=2e-10)[0]
                elwert = hi/lo*(-math.expm1(-6.283185*z/hi))/(-math.expm1(-6.283185*z/lo))
                return math.exp(-x)*born*elwert
            answer += ions[z-1]*z*z*quad(integrand,0.,70.,epsabs=1e-12,epsrel=3e-9,limit=200)[0]
        answer *= 1.7337*aune/math.sqrt(rydt)/w**3
        sampled = freefree_subtraction(t,ne,ions,epa,np.array([photon_u]),64)[0]
        check(f'adaptive_thermal_u{photon_u:g}',sampled/answer-1,.003)
    # The supplied fits and D>5 series are approximate. Record the difference,
    # while independently checked quadratures above establish numerical accuracy.
    reference = legacy_factor(executable,t,ne,ions,u)
    direct = scattering_factor(t,ne,ions,u)
    scattering_comparison = dict(temperature_K=t,electron_density_cm3=ne,
        u=u.tolist(),source_factors=reference.tolist(),direct_factors=direct.tolist(),
        relative_differences=(direct/reference-1).tolist())
    return results, scattering_comparison


def native_task(task):
    label, it, jn, sigma, ions, epa, atom_count = task
    dv, u = read_mesh(DATA/'m01.mesh')
    t = 10**(.025*it); ne = 10**(.25*jn)
    arrays = {}; means = []; invalid = []
    base = float(A0*A0/MU*native_rosseland(sigma,dv))
    for order in (3,32,64):
        correction = corrections(t,ne,ions,epa,u,order)
        total = sigma+correction['scattering_change']+correction['freefree_change']
        positive = bool(np.all(total > 0))
        if positive:
            kap = float(A0*A0/MU*native_rosseland(total,dv))
            means.append(dict(order=order,rosseland_cm2_g=kap,relative_to_uncorrected=kap/base-1))
        else:
            mask = total <= 0
            invalid.append(dict(order=order,points=int(mask.sum()),u_min=float(u[mask].min()),
                                u_max=float(u[mask].max()),minimum=float(total.min())))
        if order == 64:
            arrays=dict(total=total,free_electron_scattering=correction['scattering'],
                        freefree_change=correction['freefree_change'])
            response=correction['electron_response']
    if len(means)==3:
        convergence=means[2]['rosseland_cm2_g']/means[1]['rosseland_cm2_g']-1
    else:convergence=None
    row = dict(label=label,temperature_index=it,electron_index=jn,temperature_K=t,
               electron_density_cm3=ne,electrons_per_reference=epa,
               atom_count_per_reference=atom_count,positive_ion_counts_per_reference=ions.tolist(),
               electron_response=response,uncorrected_rosseland_cm2_g=base,
               means=means,invalid_corrected_spectra=invalid,quadrature_mean_change=convergence,
               used_in_mean_interpolation=not invalid,
               archival_screening_skip_would_apply=bool(5.8804e-19*ne*(ions @ np.arange(1.,29.)**2)/(epa*t*t)<5e-8))
    return row, arrays


def main():
    if OUTPUT.exists():raise FileExistsError(OUTPUT)
    inputs = {str(p):digest(p) for p in [Path(__file__).resolve(),
        ROOT/'scripts/op_mixture_corrections.py',ROOT/'scripts/op_native_spectra_v2.py',
        ROOT/'scripts/compare_op_native_mixtures.py',BASE,OPCODE,DATA/'m01.mesh']}
    previous=json.loads(BASE.read_text())
    for p,h in previous['input_sha256'].items():
        if digest(p)!=h:raise ValueError('retained mixture input changed: '+p)
    executable, compilation=reference_program()
    checks, reference=controls(executable)
    tasks=[]; source_states={}
    for case in previous['cases']:
        weights={int(z):v for z,v in case['included_atom_weights'].items()}
        artifact=case['spectral_artifact']; p=Path(artifact['path'])
        if digest(p)!=artifact['sha256']:raise ValueError('retained mixture spectra changed')
        inputs[str(p)]=digest(p)
        with np.load(p) as original:
            for plane in case['planes']:
                it=plane['temperature_index']
                selected=sorted(set(j for r in plane['density_interpolations'] for j in r['electron_indices']))
                ions={jn:np.zeros(28) for jn in selected}; electrons={jn:0. for jn in selected}
                for z, weight in weights.items():
                    key=(z,it)
                    if key not in source_states:
                        src=DATA/f'm{z:02d}.{it:03d}'; inputs[str(src)]=digest(src)
                        source_states[key]={s['electron_index']:s for s in read_spectra(src)['states']}
                    for jn in selected:
                        state=source_states[key][jn]
                        if not state['positive_cross_section']:raise ValueError('selected original state is masked')
                        electrons[jn]+=weight*state['electrons_per_atom']
                        for index,fraction in zip(state['ion_indices'],state['ion_fractions']):
                            charge=z-1-int(index)
                            if charge>0:ions[jn][charge-1]+=weight*fraction
                for jn in selected:
                    tasks.append((case['label'],it,jn,original[f't{it}_n{jn}'].copy(),ions[jn],electrons[jn],sum(weights.values())))
    records=[]; assets={}
    with ProcessPoolExecutor(max_workers=4) as pool:
        for row,arrays in pool.map(native_task,tasks):
            records.append(row)
            label=f"{row['label']}_t{row['temperature_index']}_n{row['electron_index']}"
            for name,value in arrays.items():assets[label+'_'+name]=value
    output_spectra=WORK/'corrected_native_spectra.npz'
    np.savez_compressed(output_spectra,**assets)
    lookup={(r['label'],r['temperature_index'],r['electron_index']):r for r in records}
    cases=[]
    for case in previous['cases']:
        mean_planes=[]; excluded=[]
        for plane in case['planes']:
            it=plane['temperature_index']; density_rows=[]
            for original in plane['density_interpolations']:
                native=[lookup[case['label'],it,j] for j in original['electron_indices']]
                if any(not r['used_in_mean_interpolation'] for r in native):
                    excluded.append(dict(temperature_index=it,density_points=original['density_points'],reason='nonpositive_corrected_source'))
                    continue
                values=[r['means'][-1]['rosseland_cm2_g'] for r in native]
                k=math.exp(interpolation(np.log(original['native_densities']),np.log(values),math.log(case['density_g_cm3'])))
                density_rows.append(dict(density_points=original['density_points'],rosseland_cm2_g=k))
            mean_planes.append(dict(temperature_K=plane['temperature_K'],rows=density_rows))
        final=[]
        for nt in (2,4):
            for nr in (2,4):
                planes=[mean_planes[i] for i in ([1,2] if nt==2 else [0,1,2,3])]
                available=[next((r for r in p['rows'] if r['density_points']==nr),None) for p in planes]
                if any(r is None for r in available):continue
                k=math.exp(interpolation(np.log([p['temperature_K'] for p in planes]),np.log([r['rosseland_cm2_g'] for r in available]),math.log(case['temperature_K'])))
                original=next(r['rosseland_cm2_g'] for r in case['interpolations'] if r['temperature_points']==nt and r['density_points']==nr)
                row=dict(temperature_points=nt,density_points=nr,rosseland_cm2_g=k,relative_to_uncorrected=k/original-1)
                if 'tops_on_OP_mesh_cm2_g' in case:row['relative_to_TOPS_on_OP_mesh']=k/case['tops_on_OP_mesh_cm2_g']-1
                final.append(row)
        cases.append(dict(label=case['label'],interpolations=final,excluded_interpolations=excluded))
    for p,h in inputs.items():
        if digest(p)!=h:raise ValueError('input changed during calculation: '+p)
    all_pass=all(r['passed'] for r in checks)
    result=dict(scope=__doc__,outcome='completed_conditional_comparison' if all_pass else 'numerical_checks_failed',
        accepted_for_stellar_opacity=False,input_sha256=inputs,compilation=compilation,
        controls=checks,source_scattering_comparison=reference,native_states=records,cases=cases,
        spectral_artifact=dict(path=str(output_spectra),bytes=output_spectra.stat().st_size,sha256=digest(output_spectra)),
        limitations=['Static long-wavelength ring scattering; finite-momentum exchange and full dynamic ring terms are not evaluated.',
            'Free-free correction retains OP ionic Debye screening, Born/Elwert approximation and classical thermal populations. No degenerate or correlated-plasma accuracy is established.',
            'Four trace elements remain absent; bound-electron scattering is not separately identified.',
            'Any nonpositive corrected spectrum is excluded as a complete state; no clipping, extrapolation or zero filling.',
            'The supplied first-frequency overwrite and coarse frequency interpolation are not used; source constants remain at their published precision.',
            'Means are on the finite OP mesh; no new plasma refraction, atmosphere or stellar selection is made.'])
    OUTPUT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],controls=len(checks),failed=[r for r in checks if not r['passed']],native_states=len(records),
        excluded_states=sum(not r['used_in_mean_interpolation'] for r in records),cases=cases,
        maximum_quadrature_mean_change=max((abs(r['quadrature_mean_change']) for r in records if r['quadrature_mean_change'] is not None),default=None))),flush=True)
    if not all_pass:raise SystemExit(1)


if __name__=='__main__':main()
