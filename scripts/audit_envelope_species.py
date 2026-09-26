#!/usr/bin/env python3
"""Check full-star species exchange with a homogeneous outer reservoir.

Thermal structure, mass mesh and the convective boundary are fixed. No cold
microscopic collision law or full stellar heat closure is selected here.
"""
import argparse, csv, hashlib, json, subprocess, time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('probe',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--scratch',type=Path,required=True)
    p.add_argument('--reuse-raw',type=Path)
    a=p.parse_args();assert not a.output.exists() and not a.scratch.exists();a.scratch.mkdir()
    root=Path(__file__).resolve().parents[1]
    eos=Path('/tmp/ember-smooth-eos-refined-family-v1/freeeos300_gs98_z020.dat')
    collision=Path('/tmp/ember-collision-transport-table-v1.dat')
    profile=root/'docs/reports/2026-09-11/evolution_latest_profile.csv'
    prior_manifest=json.loads(Path('/tmp/ember-full-star-connection-v2/inputs.json').read_text())['input_sha256']
    # Reuse the exact native EOS family from the completed whole-star control.
    eos_inputs={k:v for k,v in prior_manifest.items() if '/eos/' in k or 'eos-refined-family' in k}
    for name,value in eos_inputs.items():assert sha(name)==value,name
    assert sha(profile)=='28ab0b571dbc495ff2f48a6c6e7cd79f09e5711f3b725a2f478d4b0390d7ef2f'
    assert sha(eos)=='7d6bfa0b26264834307dad651852bc5ff693f91f3861b4934e582dee047fd7f5'
    assert sha(collision)=='48f7328da79a6163a997ab524c623ecb334d964f47bb7fbdefc084e63b46e90f'
    sources=[Path(__file__),a.probe,profile,eos,collision,
             root/'scripts/envelope_species_probe.cpp',root/'src/species_transport.cpp',
             root/'include/ember/species_transport.hpp',root/'src/evolution.cpp',
             root/'src/screened_microscopic_transport.cpp',root/'src/eos_smooth_mixture.cpp',
             root/'src/nuclear_pp.cpp',Path('/tmp/ember-envelope-flux-build-v1/src/libember.a')]
    identities={**eos_inputs,**{str(p.resolve()):sha(p) for p in sources}}
    with profile.open() as f:rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
    keys=['mass_g','radius_cm','density_g_cm3','temperature_K','X','Y3']
    data=np.array([[r[k] for k in keys] for r in rows]);assert len(data)==512
    # Retained boundary: isolated nodes 0..396, homogeneous nodes 397..511.
    ends=list(range(1,398))+[512];starts=[0]+ends[:-1]
    cases=[dict(ions=i,years=y,transport=1) for i in [0,1] for y in [1e6,2.2e7]]
    cases.append(dict(ions=1,years=1e6,transport=0))
    queries=[]
    for c in cases:
        c['seconds']=c['years']*365.25*86400
        q=[512,len(ends),c['seconds'],1,c['ions'],-1,c['transport']]+ends+data.ravel().tolist()
        queries.append(' '.join(format(float(x),'.17g') for x in q))
    payload='\n'.join(queries)+'\n';(a.scratch/'queries.txt').write_text(payload)
    command=[str(a.probe),str(eos),str(collision)]
    (a.scratch/'inputs.json').write_text(json.dumps(dict(input_sha256=identities,command=command,cases=cases),indent=2)+'\n')
    begin=time.monotonic();outpath=a.scratch/'responses.jsonl';errpath=a.scratch/'stderr.txt'
    if a.reuse_raw:
        old=json.loads((a.reuse_raw/'inputs.json').read_text())
        assert old['command']==command and old['cases']==cases
        assert (a.reuse_raw/'queries.txt').read_text()==payload
        failure=json.loads((root/'docs/results/envelope_species_v1.json').read_text())
        for name,value in failure['artifacts_sha256'].items():assert sha(name)==value,name
        for name,value in old['input_sha256'].items():
            source=Path(name)
            if source==Path(__file__).resolve():source=Path(failure['original_auditor']['path'])
            assert sha(source)==value,name
        outpath.write_bytes((a.reuse_raw/'responses.jsonl').read_bytes())
        errpath.write_bytes((a.reuse_raw/'stderr.txt').read_bytes());code=0
    else:
        with (a.scratch/'queries.txt').open() as inp,outpath.open('x') as out,errpath.open('x') as err:
            try:code=subprocess.run(command,stdin=inp,stdout=out,stderr=err,timeout=200).returncode
            except subprocess.TimeoutExpired:code=124
    elapsed=time.monotonic()-begin
    records=[json.loads(line) for line in outpath.read_text().splitlines()]
    failures=[];summaries=[];maximum={}
    if code or len(records)!=len(cases):failures.append(dict(process_exit=code,replies=len(records),expected=len(cases)))
    def check(name,value,tolerance,case):
        value=float(value);maximum[name]=max(maximum.get(name,0.),value)
        if not np.isfinite(value) or value>tolerance:failures.append(dict(case=case,check=name,value=value,tolerance=tolerance))
    dm=np.diff(data[:,0]);weights=np.r_[data[0,0],np.zeros(511)];weights[:-1]+=.5*dm;weights[1:]+=.5*dm
    M=data[-1,0];w=weights/M
    for ci,(c,r) in enumerate(zip(cases,records)):
        if 'error' in r:failures.append(dict(case=ci,error=r['error']));continue
        x=np.array(r['composition']);nuclear=np.array(r['nuclear']);nj=np.array(r['nuclear_jacobian']);dt=c['seconds']
        nr=len(ends);matrix=lil_matrix((2*nr,2*nr));residual=np.zeros((nr,2))
        for region,(start,end) in enumerate(zip(starts,ends)):
            if not np.all(x[start:end]==x[start]):failures.append(dict(case=ci,check='region uniformity'))
            residual[region]=np.sum(w[start:end,None]*(x[start:end]-data[start:end,4:6]-dt*nuclear[start:end]),axis=0)
            matrix[2*region:2*region+2,2*region:2*region+2]=w[start:end].sum()*np.eye(2)-dt*np.sum(w[start:end,None,None]*nj[start:end],axis=0)
        total=np.array(r['total_faces']);assert total.shape==(511,2)
        for region,f in enumerate(r['faces']):
            assert f['index']==ends[region]-1
            assert np.array_equal(total[f['index']],np.array(f['rate']))
            rate=dt/M*np.array(f['rate']);left=dt/M*np.array(f['dleft']);right=dt/M*np.array(f['dright'])
            residual[region]+=rate;residual[region+1]-=rate
            matrix[2*region:2*region+2,2*region:2*region+2]+=left
            matrix[2*region:2*region+2,2*region+2:2*region+4]+=right
            matrix[2*region+2:2*region+4,2*region:2*region+2]-=left
            matrix[2*region+2:2*region+4,2*region+2:2*region+4]-=right
        correction=spsolve(matrix.tocsr(),-residual.ravel())
        check('independent_newton_correction',np.max(abs(correction)),5e-13,ci)
        # Independent all-cell continuity residual, including interior faces
        # that cancel out of the homogeneous reservoir's global equation.
        closed=np.vstack([np.zeros(2),total,np.zeros(2)]).astype(np.longdouble)
        changes=w.astype(np.longdouble)[:,None]*(x.astype(np.longdouble)-data[:,4:6].astype(np.longdouble)-np.longdouble(dt)*nuclear.astype(np.longdouble))
        local=changes+np.longdouble(dt/M)*np.diff(closed,axis=0)
        check('all_cell_mass_normalized_balance',np.max(abs(local)),2e-13,ci)
        check('global_species_balance',np.max(abs(changes.sum(axis=0))),2e-13,ci)
        check('reported_cell_balance_difference',np.max(abs(local-np.array(r['cell_balances']))),2e-15,ci)
        check('reported_region_balance_difference',np.max(abs(residual-np.array(r['region_balances']))),2e-15,ci)
        if (c['transport'] and r['maximum_requested_face']!=396) or (not c['transport'] and r['maximum_requested_face']!=0):
            failures.append(dict(case=ci,check='unexpected microscopic callback domain'))
        assert np.all(x>=0) and np.all(x.sum(axis=1)<.98)
        face_enthalpy=np.array(r['eos_face_enthalpy']);thermal_L=rows[-1]['luminosity_erg_s']
        hdot_total=np.sum(face_enthalpy*total,axis=1)
        hot=[]
        for f in r['hot_internal']:
            i=f['index'];micro=np.array(f['rate']);mix=total[i]-micro;h=face_enthalpy[i]
            hot.append(dict(face=i,total_rate=total[i].tolist(),microscopic_rate=micro.tolist(),required_redistribution=mix.tolist(),
                            eos_enthalpy_flux_total_over_L=float(h@total[i]/thermal_L),
                            microscopic_carried_over_L=f['carried']/thermal_L,
                            eos_redistribution_enthalpy_over_L=float(h@mix/thermal_L),
                            conditional_sum_over_L=float((f['carried']+h@mix)/thermal_L)))
        if c['transport']:assert [h['face'] for h in hot]==list(range(397,421))
        summaries.append(dict(case=ci,**c,elapsed_solve_seconds=r['seconds'],iterations=r['iterations'],
                              envelope_X=float(x[-1,0]),envelope_delta_X=float(x[-1,0]-data[-1,4]),
                              base_H_rate=float(total[396,0]),base_He3_rate=float(total[396,1]),
                              maximum_internal_eos_enthalpy_flux_over_L=float(np.max(abs(hdot_total[397:]))/thermal_L),
                              global_balance=[float(v) for v in changes.sum(axis=0)],hot_internal=hot))
    assert all(sha(name)==value for name,value in identities.items())
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed_fixed_structure_whole_star_species_controls' if not failures else 'failed',
                accepted_for_stellar_evolution=False,selected_star_unchanged=True,points=512,mixed_envelope_first_node=397,
                cases=cases,elapsed_seconds=elapsed,maxima=maximum,failures=failures,summaries=summaries,
                new_native_cases=0 if a.reuse_raw else len(cases),
                reused_raw_directory=str(a.reuse_raw) if a.reuse_raw else None,
                input_manifest=dict(path=str(a.scratch/'inputs.json'),sha256=sha(a.scratch/'inputs.json')),
                artifacts_sha256={str(p):sha(p) for p in a.scratch.iterdir()},
                limitations=['Fixed thermal structure, geometry and convective boundary; no stellar age or energy update.',
                             'Conditional fully stripped collision physics at hot region boundaries only; no cold collision law.',
                             'The inferred internal flux is total microscopic plus convective redistribution, obtained from species continuity.',
                             'EOS enthalpy times inferred redistribution is a diagnostic. Adding it to the existing heat equation requires independent physical accounting.',
                             'No partial-ionization error bound, finite convective mixing rate or timestep accuracy claim.'])
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ['outcome','elapsed_seconds','maxima','failures']},indent=2))
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
