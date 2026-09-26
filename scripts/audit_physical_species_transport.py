#!/usr/bin/env python3
"""Check implicit burning with changing physical diffusion coefficients.

Old zeros are conserved inputs, not EOS force queries. Interior Newton guesses
are conservative mixtures of the burn/mix predictor and its global mean.
Their blend is varied independently to establish that it does not determine
accepted abundances. Globally absent species are removed explicitly.
Thermal structure and mesh are fixed; no stellar age or energy update occurs.
"""
import argparse,csv,gzip,hashlib,json,math,subprocess,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reuse-raw',type=Path)
    p.add_argument('probe',type=Path);p.add_argument('output',type=Path);p.add_argument('--scratch',required=True,type=Path)
    args=p.parse_args();assert not args.output.exists() and not args.scratch.exists();args.scratch.mkdir()
    root=Path(__file__).resolve().parents[1];sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    family=Path('/private/tmp/ember-smooth-eos-refined-family-v1/freeeos300_gs98_z020.dat');table=Path('/tmp/ember-collision-transport-table-v1.dat')
    assert sha(family)=='7d6bfa0b26264834307dad651852bc5ff693f91f3861b4934e582dee047fd7f5'
    assert sha(table)=='48f7328da79a6163a997ab524c623ecb334d964f47bb7fbdefc084e63b46e90f'
    profile=root/'docs/reports/2026-09-11/evolution_latest_profile.csv'
    assert sha(profile)=='28ab0b571dbc495ff2f48a6c6e7cd79f09e5711f3b725a2f478d4b0390d7ef2f'
    with profile.open() as f:rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
    fields=['mass_g','radius_cm','density_g_cm3','temperature_K','X','Y3']
    star=np.array([[r[k] for k in fields] for r in rows[:397]])
    n=16;mass=np.linspace(1e27/(2*n),1e27,n);rho=500.;T=8e6
    patch=np.array([[m,(3*m/(4*math.pi*rho))**(1/3),rho,T,.08 if i<n//2 else 0.,0. if i<n//2 else .006] for i,m in enumerate(mass)])
    cases=[]
    def add(kind,data,years,burning,ions,blend,transport=1,ends=None):
        cases.append(dict(kind=kind,data=data.tolist(),years=years,burning=burning,ions=ions,blend=blend,transport=transport,ends=ends or list(range(1,len(data)+1))))
    for ions in [0,1]:
        for years in [1e3,1e5,1e6,2e7]:add('saved_star',star,years,1,ions,0.)
        for blend in [-1.,.1,1.]:add('saved_star',star,1e6,1,ions,blend)
        for years in [1e3,1e5,1e7]:
            for blend in [0.,.01,.1,1.]:add('initially_empty',patch,years,0,ions,blend)
        for species in [0,1]:
            data=patch.copy();data[:,4+species]=0.;data[:,5-species]=np.linspace(.04 if species else .003,.08 if species else .006,n)
            add('globally_absent_'+str(species),data,1e5,0,ions,.1)
            add('nuclear_activation_'+str(species),data,1e5,1,ions,.1)
        data=patch.copy();data[:,4:6]=0.;data[:,3]*=np.exp(np.linspace(-.02,.02,n))
        add('pure_helium',data,1e5,0,ions,1.)
        add('mixed_empty',patch,1e5,0,ions,.1,ends=[4,8,12,16])
        add('initially_empty',patch,1e5,0,ions,-1.)
    add('burn_only',star,1e6,1,0,0.,transport=0)
    queries=[]
    for c in cases:
        d=np.array(c['data']);dt=c['years']*365.25*86400
        q=[len(d),len(c['ends']),dt,c['burning'],c['ions'],c['blend'],c['transport']]+c['ends']+d.ravel().tolist()
        queries.append(' '.join(format(float(x),'.17g') for x in q)+'\n')
    payload=''.join(queries);(args.scratch/'queries.txt.gz').write_bytes(gzip.compress(payload.encode(),mtime=0))
    if args.reuse_raw:
        prior=json.loads((args.reuse_raw/'source_manifest.json').read_text())
        for name,identity in prior.items():assert sha(name)==identity,name
        assert gzip.decompress((args.reuse_raw/'queries.txt.gz').read_bytes()).decode()==payload
        stdout=gzip.decompress((args.reuse_raw/'responses.jsonl.gz').read_bytes()).decode()
        stderr=(args.reuse_raw/'stderr.txt').read_text();exit_code=0;elapsed=None
    else:
        start=time.monotonic()
        try:
            proc=subprocess.run([str(args.probe),str(family),str(table)],input=payload,text=True,capture_output=True,timeout=160)
            stdout,stderr,exit_code=proc.stdout,proc.stderr,proc.returncode
        except subprocess.TimeoutExpired as e:
            stdout=e.stdout or b'';stderr=e.stderr or b''
            if isinstance(stdout,bytes):stdout=stdout.decode()
            if isinstance(stderr,bytes):stderr=stderr.decode()
            exit_code=124
        elapsed=time.monotonic()-start
    (args.scratch/'responses.jsonl.gz').write_bytes(gzip.compress(stdout.encode(),mtime=0));(args.scratch/'stderr.txt').write_text(stderr)
    replies=[json.loads(x) for x in stdout.splitlines()];failures=[];summaries=[];groups={}
    maximum=dict(independent_correction=0.,global_balance=0.,reported_balance_difference=0.,guess_balance=0.,guess_dependence=0.)
    if exit_code or len(replies)!=len(cases):failures.append(dict(process_exit=exit_code,replies=len(replies),expected=len(cases)))
    def check(name,error,tolerance,case):
        maximum[name]=max(maximum[name],float(error))
        if not math.isfinite(error) or error>tolerance:failures.append(dict(case=case,check=name,error=float(error),tolerance=tolerance))
    for ci,(c,r) in enumerate(zip(cases,replies)):
        if 'error' in r:failures.append(dict(case=ci,kind=c['kind'],error=r['error']));continue
        data=np.array(c['data']);n=len(data);M=data[-1,0];dm=np.diff(data[:,0]);w=np.r_[data[0,0],np.zeros(n-1)];w[:-1]+=.5*dm;w[1:]+=.5*dm;w/=M
        x=np.array(r['composition']);source=np.array(r['nuclear']);source_j=np.array(r['nuclear_jacobian']);dt=c['years']*365.25*86400
        ends=c['ends'];starts=[0]+ends[:-1];nr=len(ends);jac=lil_matrix((2*nr,2*nr));res=np.zeros((nr,2));xr=[]
        for k,(begin,end) in enumerate(zip(starts,ends)):
            fraction=float(w[begin:end].sum());old=np.sum(w[begin:end,None]*data[begin:end,4:6],axis=0)
            xr.append(x[begin]);res[k]=fraction*x[begin]-old-dt*np.sum(w[begin:end,None]*source[begin:end],axis=0)
            jac[2*k:2*k+2,2*k:2*k+2]=fraction*np.eye(2)-dt*np.sum(w[begin:end,None,None]*source_j[begin:end],axis=0)
            if not np.all(x[begin:end]==x[begin]):failures.append(dict(case=ci,check='mixed_region_uniformity'))
        for k,f in enumerate(r['faces']):
            if f['index']!=ends[k]-1:failures.append(dict(case=ci,check='face_index'))
            flux=dt/M*np.array(f['rate']);a=dt/M*np.array(f['dleft']);b=dt/M*np.array(f['dright'])
            res[k]+=flux;res[k+1]-=flux
            jac[2*k:2*k+2,2*k:2*k+2]+=a;jac[2*k:2*k+2,2*k+2:2*k+4]+=b
            jac[2*k+2:2*k+4,2*k:2*k+2]-=a;jac[2*k+2:2*k+4,2*k+2:2*k+4]-=b
        correction=spsolve(jac.tocsr(),-res.ravel());check('independent_correction',float(np.max(abs(correction))),5e-13,ci)
        # Sum the mass changes directly, without cancelling face fluxes.
        mass_change=np.sum(w.astype(np.longdouble)[:,None]*(x-data[:,4:6]).astype(np.longdouble),axis=0)
        burn=dt*np.sum(w.astype(np.longdouble)[:,None]*source.astype(np.longdouble),axis=0)
        balance=mass_change-burn
        check('global_balance',float(np.max(abs(balance))),2e-13,ci)
        check('reported_balance_difference',float(np.max(abs(balance-np.array(r['balance'])))),2e-13,ci)
        check('guess_balance',max(abs(v) for v in r['guess_balance']),2e-14,ci)
        if not np.all(np.isfinite(x)) or np.any(x<0) or np.any(x.sum(axis=1)>=.98):failures.append(dict(case=ci,check='composition_domain'))
        for k,active in enumerate(r['active']):
            if active and not np.all(x[:,k]>0):failures.append(dict(case=ci,check='active_positive',species=k))
            if not active and np.any(x[:,k]!=0):failures.append(dict(case=ci,check='absent_stays_absent',species=k))
        if c['kind']=='initially_empty' and c['blend'] in [-1.,0.] and not r['automatic_seed']:failures.append(dict(case=ci,check='missing_automatic_seed'))
        if c['kind'] in ['initially_empty','mixed_empty'] and r['new_phi']>r['old_phi']+1e-5:failures.append(dict(case=ci,check='isothermal_free_energy'))
        key=(c['kind'],c['years'],c['burning'],c['ions'],c['transport'],tuple(c['ends']))
        if key in groups:check('guess_dependence',float(np.max(abs(x-groups[key]))),3e-12,ci)
        else:groups[key]=x
        summaries.append(dict(case=ci,kind=c['kind'],years=c['years'],ions=c['ions'],blend=c['blend'],active=r['active'],iterations=r['iterations'],seconds=r['seconds'],
            automatic_seed=r['automatic_seed'],maximum_abundance_change=float(np.max(abs(x-data[:,4:6]))),minimum_abundance=x.min(axis=0).tolist(),integrated_balance=[float(v) for v in balance],independent_correction=float(np.max(abs(correction))),old_phi=r['old_phi'],new_phi=r['new_phi']))
    sources=[Path(__file__),args.probe,profile,family,table]+[root/p for p in ['src/species_transport.cpp','include/ember/species_transport.hpp','src/screened_microscopic_transport.cpp','include/ember/screened_microscopic_transport.hpp','include/ember/microscopic_transport.hpp','scripts/physical_species_probe.cpp','src/evolution.cpp','src/eos_smooth_mixture.cpp','src/collision_transport.cpp','src/nuclear_pp.cpp']]
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed_changing_physical_species_transport' if not failures else 'failed',accepted_for_stellar_evolution=False,
        cases=len(cases),completed_cases=len(summaries),maxima=maximum,failures=failures,summaries=summaries,inputs=cases,
        elapsed_seconds=elapsed,probe_timing=stderr,physical_step_seconds=sum(s['seconds'] for s in summaries),
        input_sha256={str(p):sha(p) for p in sources},artifacts_sha256={str(p):sha(p) for p in args.scratch.iterdir()},
        limitations=['Fixed thermal structure and sealed boundaries; no energy, atmosphere, hydrostatic, convective-boundary or age update.',
          'Fully stripped conditional screened collision model on T>=2e6 K; selected star and physical inputs unchanged.',
          'An initially absent local species is handled by an interior NEW-time guess and the original OLD abundance RHS. The instantaneous flux at a zero NEW-time endpoint is not evaluated.',
          'Positive seed is only an initial iterate. Its blend is varied, its species inventory is conserved, and the final equations are independently reconstructed.',
          'Globally absent species are explicitly removed only after the burning predictor; nuclear production can activate them.'])
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k in ['outcome','cases','completed_cases','maxima','failures','physical_step_seconds','probe_timing']}))
    if failures:raise SystemExit(1)

if __name__=='__main__':main()
