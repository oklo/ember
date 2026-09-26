#!/usr/bin/env python3
"""Check changing physical C++ collision coefficients against the Python model.

Saved cases reuse the existing physical reference. New states change T, rho,
composition and screening and recalculate the Python collision response.
No EOS or expensive electron-pair scattering integrals are repeated.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import time

import numpy as np
from scipy.integrate import quad
from scipy.special import expit
from electron_electron_energy_modes_v2 import PolynomialBasis
from electron_ion_born import chemical_potential, KB, ME, HBAR, E2
from electron_mixture_collision import reduce_electron_collision, symmetric
from electron_pair_table import ElectronPairTable
from ion_collision_table import ExtendedIonTable
from ion_collision_integrals import MU, screening_length
from ion_electron_transport import coupled_elastic_operator, heat_decomposition


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('probe', type=Path)
    parser.add_argument('table', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--scratch', type=Path, required=True)
    parser.add_argument('--reuse-raw', type=Path)
    args = parser.parse_args()
    assert not args.output.exists() and not args.scratch.exists()
    args.scratch.mkdir()
    root = Path(__file__).resolve().parents[1]
    inputs = {}
    def read(name):
        path = root/'docs/results'/name
        inputs[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        return json.loads(path.read_text())
    saved = read('stellar_diffusion_pair_transport_v1.json')
    hot = read('stellar_electron_heat_3890gyr_v1.json')
    material = read('diffusion_material_regime_3890gyr_v1.json')
    extension = read('ion_collision_table_extension_v1.json')
    pair_path = root/'docs/results/electron_pair_table_sources_v2.json'
    pair = ElectronPairTable(pair_path)
    ions = ExtendedIonTable(extension['table'])
    assert hashlib.sha256(Path(extension['table']).read_bytes()).hexdigest() == extension['table_sha256']
    header = root/'include/ember/gs98_mixture.hpp'
    metals = [tuple(map(float, s.split(','))) for s in re.findall(r'^\s*\{([^{}]+)\}, //', header.read_text(), re.M)]
    all_mass = np.array([1.,3.,4.]+[m[1] for m in metals])
    all_charge = np.array([1.,2.,2.]+[m[0] for m in metals])
    def mixture(state):
        T,rho,X,Y3,Z = state
        x = np.array([X,Y3,1-X-Y3-Z]+[Z*m[3] for m in metals])
        return x, rho/MU*x/all_mass
    def screened(state, include_ions):
        T,rho,X,Y3,Z=state
        x,n = mixture(state)
        ne = float(n@all_charge)
        eta = chemical_potential(ne,T)
        upper=math.sqrt(max(eta,0)+60)
        derivative=quad(lambda t:2*t*t*expit(eta-t*t)*expit(t*t-eta),0,upper,epsabs=1e-100,epsrel=2e-12)[0]
        scale=(2*ME*KB*T)**1.5/(2*math.pi**2*HBAR**3)
        stiffness=KB*T*ne/(scale*derivative)
        length=(screening_length(n,all_charge,T,electron_stiffness_erg=stiffness)['length_cm'] if include_ions
                else math.sqrt(stiffness/(4*math.pi*E2*ne)))
        return length,stiffness
    def reference(state,length):
        T,rho,X,Y3,Z=state
        x,_ = mixture(state)
        keep=x>0
        indices=np.flatnonzero(keep)
        active=[int(np.flatnonzero(indices==i)[0]) for i in [0,1] if keep[i]]
        outindices=[i for i in [0,1] if keep[i]]
        assert active
        ref=int(np.flatnonzero(indices==2)[0])
        x,mass,z=x[keep],all_mass[keep],all_charge[keep]
        n=rho/MU*x/mass
        eta=chemical_potential(float(n@z),T)
        B=8*ME*KB*T*length**2/HBAR**2
        basis=PolynomialBasis(eta,9)
        ei=symmetric(basis.ion_matrix(B)); ee=pair.matrices(eta,B)
        reduced=reduce_electron_collision(electron_ion_matrix=ei,electron_electron_matrix=ee,ion_density=n,ion_charges=z)
        reduced_mass=MU*(mass[:,None]*mass[None,:]/(mass[:,None]+mass[None,:]))
        interaction=np.outer(z,z)*E2
        moments=ions.moments(interaction/(KB*T*length))
        omega=np.sqrt(2*math.pi/reduced_mass)*interaction**2/(KB*T)**1.5*moments[:,:,0]
        op=coupled_elastic_operator(density=rho,temperature=T,mass_fractions=x,mass_numbers=mass,
            charges=z,independent_ions=active,reference_ion=ref,ion_resistance=16/3*np.outer(n,n)*reduced_mass*omega,
            ion_z=1-.4*moments[:,:,1]/moments[:,:,0],
            ion_zprime=2.5-2*moments[:,:,1]/moments[:,:,0]+.4*moments[:,:,2]/moments[:,:,0],
            ion_zdoubleprime=moments[:,:,3]/moments[:,:,0],electron_moments=ei[:2,:2],
            electron_full_response=reduced.retained_response,energy_scale=KB/MU*T)
        enthalpy,kappa=heat_decomposition(op.mobility,temperature=T,energy_scale=op.energy_scale)
        force=np.arange(1,len(outindices)+2)*.1
        sol=op.solve(force)
        speed=max(float(np.max(abs(sol['velocity']))),1e-100)
        current=abs((n*z)@sol['velocity'][:-1]-(n@z)*sol['velocity'][-1])/(2*(n@z)*speed)
        baryon=abs(x@sol['velocity'][:-1])/speed
        assert current<1e-12 and baryon<1e-12
        return dict(eta=eta,conductivity=kappa,enthalpy=enthalpy.tolist(),mobility=op.mobility.tolist(),indices=outindices+[2],
                    baryon_flux_error=baryon,electric_current_error=current)
    records=[];queries=[];expected=[]
    for sr,hr,mr in zip(saved['records'],hot['records'],material['records'],strict=True):
        assert sr['zone']==hr['zone']==mr['zone']
        state=[sr['temperature_K'],sr['density_g_cm3'],mr['X'],mr['Y3'],1-mr['X']-mr['Y3']-mr['Y4']]
        for sc,hc in zip(sr['cases'],hr['cases'],strict=True):
            assert sc['screening']==hc['screening']
            length=hc['screening_length_cm']
            queries.append(state+[length])
            records.append(dict(kind='saved',zone=sr['zone'],screening=sc['screening'],state=state))
            expected.append(dict(eta=hr['eta_nonrelativistic'],conductivity=sc['conductivity_zero_element_flux_cgs'],
                enthalpy=sc['transport_enthalpy_erg_g'],mobility=sc['mobility_scaled_heat'],indices=[0,1,2]))
    python_start=time.process_time()
    # Trial states are not on the saved trajectory. Screening follows a newly
    # integrated ideal nonrelativistic electron stiffness, passed explicitly.
    for zone in [0,80,160,240,300,360,396,410]:
        mr=material['records'][zone]
        for sign in [-1,1]:
            state=[mr['temperature_K']*math.exp(sign*.02),mr['density_baryonic_g_cm3']*math.exp(-sign*.015),
                   mr['X']*(1+sign*.04),mr['Y3']*(1-sign*.08),.02*(1+sign*.03)]
            for include_ions in [False,True]:
                length,stiffness=screened(state,include_ions)
                expected.append(reference(state,length));queries.append(state+[length,stiffness,int(include_ions)])
                records.append(dict(kind='changed_state',zone=zone,state=state,include_ions=include_ions,expected_screening_length=length))
    # Exact zeros remove species; compare with the independent active-ion
    # Python implementation. A separate pure-He control takes the trace limit.
    zero_indices=[]
    for X,Y3,Z in [(.01,0,.02),(0,.002,.02),(.01,.002,0),(.01,0,0),(0,.002,0)]:
        state=[1e7,3000,X,Y3,Z];length,stiffness=screened(state,True)
        expected.append(reference(state,length));queries.append(state+[length,stiffness,1])
        zero_indices.append(len(records));records.append(dict(kind='zero_species',state=state,expected_screening_length=length))
    trace_indices=[]
    for X in [0.,1e-10,1e-12]:
        state=[1e7,3000,X,0.,0.];length,stiffness=screened(state,True)
        expected.append(None);queries.append(state+[length,stiffness,1]);trace_indices.append(len(records))
        records.append(dict(kind='pure_helium_limit',state=state))
    invalid_indices=[]
    for state,length in [([1e7,3000,-.01,.001,.02],1e-9),([1e7,3000,.9,.2,.02],1e-9),
                         ([1e7,3000,.01,.001,.02],1e-20),([1e7,3000,.01,.001,.02],1e-4),
                         ([1e7,1e-10,.01,.001,.02],1e-9),([1e7,1e10,.01,.001,.02],1e-9),
                         ([0.,3000,.01,.001,.02],1e-9)]:
        invalid_indices.append(len(records));queries.append(state+[length]);expected.append(None)
        records.append(dict(kind='invalid',state=state))
    python_seconds=time.process_time()-python_start
    query_text=''.join(' '.join(format(float(x),'.17g') for x in q)+'\n' for q in queries)
    (args.scratch/'queries.txt').write_text(query_text)
    (args.scratch/'references.json').write_text(json.dumps(expected,indent=2)+'\n')
    if args.reuse_raw:
        assert (args.reuse_raw/'queries.txt').read_text()==query_text
        response_text=(args.reuse_raw/'response.jsonl').read_text()
        elapsed=None
    else:
        start=time.monotonic()
        completed=subprocess.run([str(args.probe),str(args.table)],input=query_text,text=True,capture_output=True,timeout=90,check=True)
        elapsed=time.monotonic()-start
        response_text=completed.stdout
    (args.scratch/'response.jsonl').write_text(response_text)
    actual=[json.loads(s) for s in response_text.splitlines()]
    assert len(actual)==len(queries)
    maxima=dict(mobility=0.,conductivity=0.,enthalpy=0.,eta=0.,screening=0.,backward_error=0.)
    failures=[]
    for index,(record,ref,value) in enumerate(zip(records,expected,actual,strict=True)):
        record['cpp']=value
        if index in invalid_indices:
            if 'error' not in value:failures.append(dict(index=index,reason='invalid query accepted'))
            continue
        if 'error' in value:
            failures.append(dict(index=index,reason=value['error']));continue
        maxima['backward_error']=max(maxima['backward_error'],value['backward_error'])
        if ref is None:continue
        record['reference']=ref
        ids=ref['indices'];matrix=np.asarray(value['mobility'])[np.ix_(ids,ids)]
        reference_matrix=np.asarray(ref['mobility'])
        diagonal=np.sqrt(np.diag(reference_matrix))
        normalized=reference_matrix/diagonal[:,None]/diagonal[None,:]
        assert np.max(abs(normalized-normalized.T))<1e-12
        # The reference matrix is G @ variables, with only roundoff asymmetry.
        # Its reciprocal part defines the quadratic form. Do not divide that
        # roundoff by the much smaller C++/Python difference being measured.
        reference_matrix=(reference_matrix+reference_matrix.T)/2
        chol=np.linalg.cholesky(reference_matrix)
        difference=np.linalg.solve(chol,matrix-reference_matrix)
        difference=np.linalg.solve(chol,difference.T).T
        assert np.max(abs(difference-difference.T))<1e-12
        operator_error=float(np.max(abs(np.linalg.eigvalsh((difference+difference.T)/2))))
        errors=dict(mobility=operator_error,
                    conductivity=abs(value['conductivity']/ref['conductivity']-1),
                    enthalpy=float(np.max(abs((np.array(value['enthalpy'])[ids[:-1]]-ref['enthalpy'])/
                                          np.maximum(abs(np.array(ref['enthalpy'])),1.)))),
                    eta=abs(value['eta']-ref['eta']))
        if 'expected_screening_length' in record:
            errors['screening']=abs(value['screening_length']/record['expected_screening_length']-1)
        record['errors']=errors
        for k,v in errors.items():maxima[k]=max(maxima[k],v)
        if max(errors.values())>2e-8:failures.append(dict(index=index,reason='comparison tolerance',errors=errors))
        absent=set([0,1])-set(ids)
        for a in absent:
            if any(value['mobility'][a]) or any(row[a] for row in value['mobility']):
                failures.append(dict(index=index,reason='absent species mobility is nonzero'))
    trace=[actual[i] for i in trace_indices]
    trace_error=None
    if all('error' not in t for t in trace):
        trace_error=max(abs(t['conductivity']/trace[0]['conductivity']-1) for t in trace[1:])
        if trace_error>1e-8:failures.append(dict(reason='pure helium trace limit',error=trace_error))
        assert trace[0]['active']==[0,0]
    source_files=['include/ember/collision_transport.hpp','src/collision_transport.cpp',
                  'scripts/collision_transport_probe.cpp','scripts/export_collision_transport.py',
                  'scripts/audit_collision_transport_cpp.py','scripts/electron_electron_energy_modes_v2.py',
                  'scripts/electron_mixture_collision.py','scripts/ion_electron_transport.py',
                  'scripts/electron_elastic_mixture.py','scripts/electron_pair_table.py','scripts/ion_collision_table.py',
                  'scripts/electron_ion_born.py','scripts/electron_ion_heat.py','scripts/ion_collision_integrals.py',
                  'scripts/diffusion_composition_forces.py']
    for path in [header,pair_path,args.probe,args.table,args.table.with_suffix('.json')]+[root/p for p in source_files]:
        inputs[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed_physical_evaluator_comparisons' if not failures else 'failed_comparisons',
        accepted_for_stellar_evolution=False,maximum_errors=maxima,failures=failures,comparison_tolerance=2e-8,
        records=records,saved_cases=844,changed_state_cases=32,zero_species_cases=5,pure_helium_cases=3,invalid_cases=len(invalid_indices),
        pure_helium_trace_conductivity_error=trace_error,python_new_state_cpu_seconds=python_seconds,
        cpp_process_elapsed_seconds=elapsed,reused_cpp_output_from=str(args.reuse_raw) if args.reuse_raw else None,
        cpp_successful_evaluation_seconds=sum(v.get('seconds',0) for v in actual),
        cpp_saved_evaluation_seconds=sum(v.get('seconds',0) for v in actual[:844]),new_EOS_queries=0,new_pair_collision_integrals=0,
        input_sha256=inputs,raw_artifacts={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in args.scratch.iterdir()},
        limitations=['Conditional fully stripped, stationary-metal kinetic model; changing ionization and correlation physics are not supplied.',
                     'Screening stiffness is an explicit input; changed-state controls use ideal nonrelativistic electron compressibility.',
                     'No stellar evolution selection, thermal/composition derivatives or atmosphere matching is established by these local tests.'])
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['records','input_sha256','raw_artifacts','limitations']},indent=2))
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
