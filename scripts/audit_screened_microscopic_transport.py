#!/usr/bin/env python3
"""Check physical hot-face transport, its Jacobians, and finite chemical parts.

One native EOS load serves every query. The positive-abundance face law is
compared with the independently implemented full material matrix. Exact-zero
endpoint transport is deliberately rejected until its discretization is ready.
The regular chemical potential needed for that discretization is checked at
zero directly and by one-sided differences, without an abundance floor.
"""
import argparse,csv,gzip,hashlib,json,math,subprocess,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('probe',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--scratch',required=True,type=Path);p.add_argument('--reuse-raw',type=Path);p.add_argument('--reuse-report',type=Path)
    a=p.parse_args();assert not a.output.exists() and not a.scratch.exists();a.scratch.mkdir()
    root=Path(__file__).resolve().parents[1];sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    family=Path('/private/tmp/ember-smooth-eos-refined-family-v1/freeeos300_gs98_z020.dat')
    table=Path('/tmp/ember-collision-transport-table-v1.dat')
    assert sha(family)=='7d6bfa0b26264834307dad651852bc5ff693f91f3861b4934e582dee047fd7f5'
    assert sha(table)=='48f7328da79a6163a997ab524c623ecb334d964f47bb7fbdefc084e63b46e90f'
    profile=root/'docs/reports/2026-09-11/evolution_latest_profile.csv'
    with profile.open() as f:rows=[{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]
    old_input=Path('/tmp/ember-smooth-eos-refined-physical-run-v1/retained/retained.input');old_output=old_input.with_suffix('.output')
    assert sha(old_input)=='850c2a243e1a2fe12ddaabcf1d466e3ecb761c9e3fdada10fd1b3f5751e84e77'
    assert sha(old_output)=='8bc24b90d0eba4a3933313178d600bac883765cd70b576cf3c3c7d9631e3b3a2'
    cache={}
    for q,line in zip(old_input.read_text().splitlines(),old_output.read_text().splitlines(),strict=True):
        v=list(map(float,q.split()));r=json.loads(line)
        if v[4]==1 and r['ok']:cache[tuple(v[:4])]=r['values']
    queries=[];faces=[];regular=[];stencils=[];rejects=[]
    def add(mode,q):
        i=len(queries);queries.append(mode+' '+' '.join(format(float(x),'.17g') for x in q)+'\n');return i
    def point(r):return [math.log(r['radius_cm']),math.log(r['density_g_cm3']),math.log(r['temperature_K']),r['luminosity_erg_s'],r['X'],r['Y3']]
    for ions in [0,1]:
        for i in range(421):
            q=[ions,1,1,rows[i]['mass_g'],rows[i+1]['mass_g']]+point(rows[i])+point(rows[i+1])
            j=add('face',q);faces.append(dict(index=j,query=q,kind='saved',face=i))
    # Perturb the material state as well as each independent composition.
    for ions in [0,1]:
        for i in [80,220,360]:
            q=[ions,1,1,rows[i]['mass_g'],rows[i+1]['mass_g']]+point(rows[i])+point(rows[i+1])
            for e in range(2):
                o=5+e*6;q[o+1]+=.012;q[o+2]+=.009;q[o+4]*=1.04;q[o+5]*=.93
            j=add('face',q);faces.append(dict(index=j,query=q,kind='altered',face=i))
    for ions in [0,1]:
        for i,y in [(80,.0001),(220,.003),(360,.03)]:
            q=[ions,1,1,rows[i]['mass_g'],rows[i+1]['mass_g']]+point(rows[i])+point(rows[i+1])
            q[10]=y;q[16]=1.1*y
            j=add('face',q);faces.append(dict(index=j,query=q,kind='enriched_He3',face=i))
    controls=[f for f in faces if (f['face'] in [0,40,160,300,396,420] or f['kind']!='saved')]
    # Reverse all endpoint data, with the same positive interval, to test the
    # local orientation identity. Radius need not increase in this local test.
    reverse=[]
    for f in controls:
        q=f['query'].copy();q[1:3]=[0,0];q[5:11],q[11:17]=q[11:17],q[5:11]
        reverse.append(dict(original=f['index'],index=add('face',q)))
    replay=[]
    for f in faces:
        q=f['query'].copy();q[1:3]=[0,0];replay.append(dict(original=f['index'],index=add('face',q)))
    for f in controls:
        for e in range(2):
            for k in range(6):
                q=f['query'];coordinate=5+6*e+k
                scale=1. if k<3 else (max(abs(q[coordinate]),1e25) if k==3 else q[coordinate])
                widths=[2e-3,1e-3] if f['kind']=='enriched_He3' and k==5 else [2e-4,1e-4]
                for width in widths:
                    ids=[]
                    for step in [-2,-1,1,2]:
                        trial=q.copy();trial[1:3]=[0,0];trial[coordinate]+=step*width*scale
                        ids.append(add('face',trial))
                    stencils.append(dict(kind='face',index=f['index'],endpoint=e,coordinate=k,scale=scale,width=width,level='fine' if width==widths[-1] else 'coarse',replies=ids))
    for i,row in enumerate(rows):
        q=[row['X'],row['Y3'],row['temperature_K'],row['density_g_cm3']]
        regular.append(dict(index=add('regular',q),query=q,kind='saved',zone=i))
    for T,rho in [(4157.,7.2e-5),(27183.,.037),(2341000.,13.7),(9731000.,260.7)]:
        for x,y in [(0.,.001),(.05,0.),(0.,0.)]:
            q=[x,y,T,rho];regular.append(dict(index=add('regular',q),query=q,kind='zero'))
    for r in regular:
        if r.get('zone') not in [0,160,396,493,511] and r['kind']!='zero':continue
        q=r['query']
        for k in range(4):
            scale=(q[k] if q[k]>0 else .01) if k<2 else 1.
            for width in [2e-4,1e-4]:
                steps=[0,1,2,3,4] if k<2 and q[k]==0 else [-2,-1,1,2]
                ids=[]
                for step in steps:
                    trial=q.copy()
                    if k<2:trial[k]+=step*width*scale
                    else:trial[k]*=math.exp(step*width)
                    ids.append(add('regular',trial))
                stencils.append(dict(kind='regular',index=r['index'],coordinate=k,scale=scale,width=width,replies=ids,steps=steps))
    for kind in ['zero_H','zero_He3','negative','cold','mass','nonfinite']:
        q=faces[220]['query'].copy();q[1:3]=[1,0]
        if kind=='zero_H':q[9]=0
        if kind=='zero_He3':q[10]=0
        if kind=='negative':q[9]=-.01
        if kind=='cold':q[7]=math.log(1e6)
        if kind=='mass':q[4]=q[3]
        if kind=='nonfinite':q[6]=710. # finite log, overflowing density
        rejects.append(dict(kind=kind,index=add('face',q)))
    payload=''.join(queries)
    reused={}
    if a.reuse_raw:
        assert a.reuse_report
        prior=json.loads(a.reuse_report.read_text())
        for name,identity in prior['input_sha256'].items():
            if Path(name)!=Path(__file__):assert sha(name)==identity,name
        for name,identity in prior['artifacts_sha256'].items():assert sha(name)==identity,name
        oldq=gzip.decompress((a.reuse_raw/'queries.txt.gz').read_bytes()).decode().splitlines(keepends=True)
        oldr=gzip.decompress((a.reuse_raw/'responses.jsonl.gz').read_bytes()).decode().splitlines(keepends=True)
        for q,r in zip(oldq,oldr,strict=True):reused.setdefault(q,r)
    missing=list(dict.fromkeys(q for q in queries if q not in reused))
    elapsed=None;stderr='all replies reused\n'
    if missing:
        t=time.monotonic();run=subprocess.run([str(a.probe),str(family),str(table)],input=''.join(missing),text=True,capture_output=True,timeout=160,check=True)
        elapsed=time.monotonic()-t;stderr=run.stderr
        fresh=run.stdout.splitlines(keepends=True);assert len(fresh)==len(missing)
        for q,r in zip(missing,fresh,strict=True):reused[q]=r
    stdout=''.join(reused[q] for q in queries)
    (a.scratch/'queries.txt.gz').write_bytes(gzip.compress(payload.encode(),mtime=0))
    (a.scratch/'responses.jsonl.gz').write_bytes(gzip.compress(stdout.encode(),mtime=0));(a.scratch/'stderr.txt').write_text(stderr)
    replies=[json.loads(line) for line in stdout.splitlines()];assert len(replies)==len(queries)
    failures=[];maximum=dict(full_flux=0.,scalar_replay=0.,orientation=0.,regular_identity=0.,fine_face_derivative=0.,coarse_face_derivative=0.,fine_regular_derivative=0.,coarse_regular_derivative=0.)
    checks=[];unresolved_trace_columns=[]
    def measure(label,error,tolerance,context):
        maximum[label]=max(maximum[label],float(error))
        if not math.isfinite(error) or error>tolerance:failures.append(dict(check=label,error=float(error),tolerance=tolerance,context=context))
    for i,r in enumerate(replies):
        if 'error' in r and i not in {x['index'] for x in rejects}:failures.append(dict(index=i,error=r['error']))
    for f in faces:
        r=replies[f['index']]
        if 'error' in r:continue
        rate=np.array(r['rate']);ref=np.array(r['reference_rate'])
        measure('full_flux',max(float(np.max(abs(rate-ref)/np.maximum(abs(ref),1.))),abs(r['carried']-r['reference_carried'])/max(abs(r['reference_carried']),1.)),1e-9,f['index'])
        if r['entropy']<0:failures.append(dict(check='entropy',index=f['index'],value=r['entropy']))
    for items,label,sign in [(replay,'scalar_replay',1),(reverse,'orientation',-1)]:
        for p in items:
            r,s=replies[p['original']],replies[p['index']]
            if 'error' in r or 'error' in s:continue
            v=np.r_[r['rate'],r['carried'],r['conductivity']];w=np.r_[np.array(s['rate'])*sign,s['carried']*sign,s['conductivity']]
            measure(label,float(np.max(abs(v-w)/np.maximum(abs(v),1.))),1e-10,p)
    gas=6.02214076e23*1.380649e-16
    for p in regular:
        r=replies[p['index']]
        if 'error' in r:continue
        if p['kind']=='saved':
            x,y,T,rho=p['query'];old=cache[x,y,T,rho]
            restored=np.r_[r['phi'],np.array(r['gradient'])+gas*np.log([x,y])/np.array([1.,3.]),
                (np.array(r['hessian'])+np.diag(gas/np.array([x,3*y]))).ravel(),r['gradient_T'],r['gradient_rho']]
            expected=np.array(old[21:32]);measure('regular_identity',float(np.max(abs(restored-expected)/np.maximum(abs(expected),1.))),2e-12,p['index'])
    for s in stencils:
        r=replies[s['index']];rs=[replies[i] for i in s['replies']]
        if 'error' in r or any('error' in x for x in rs):continue
        k=s['coordinate']
        if s['kind']=='face':
            side='lo' if s['endpoint']==0 else 'hi'
            if k<4:
                vec=lambda r:np.array([r['carried'],r['conductivity']])
                exact=np.array([r['dcarried_'+side][k],r['dconductivity_'+side][k]])*s['scale']
            else:
                vec=lambda r:np.array(r['rate'])
                exact=np.array(r['dleft' if side=='lo' else 'dright'])[:,k-4]*s['scale']
            values=list(map(vec,rs));reference=(values[0]-8*values[1]+8*values[2]-values[3])/(12*s['width'])
        else:
            vec=lambda r:np.array(r['gradient']);values=list(map(vec,rs))
            exact=(np.array(r['hessian'])[:,k] if k<2 else np.array(r['gradient_T' if k==2 else 'gradient_rho']))*s['scale']
            if s['steps'][0]==0:reference=(-25*values[0]+48*values[1]-36*values[2]+16*values[3]-3*values[4])/(12*s['width'])
            else:reference=(values[0]-8*values[1]+8*values[2]-values[3])/(12*s['width'])
        column_norm=np.maximum(abs(vec(r)),1.)+abs(exact)
        norm=column_norm
        column_error=float(np.max(abs(reference-exact)/column_norm))
        if s['kind']=='face' and k>=4:
            f=next(f for f in faces if f['index']==s['index'])
            q=f['query']
            # Norm of the complete species Jacobian row acting on relative
            # abundance changes. An individual trace cross term can lie below
            # the resolution of subtracting the native chemical potentials.
            norm=np.maximum(abs(vec(r)),1.)+np.sum(abs(np.array(r['dleft']))*np.array(q[9:11]),axis=1)+np.sum(abs(np.array(r['dright']))*np.array(q[15:17]),axis=1)
            if column_error>2e-6:
                unresolved_trace_columns.append(dict(**s,column_normalized_error=column_error,relative_column_size=(abs(exact)/norm).tolist()))
            if f['kind']=='enriched_He3' and column_error>2e-6:
                failures.append(dict(check='resolved_enriched_cross_derivative',error=column_error,context=s))
        error=float(np.max(abs(reference-exact)/norm))
        label=s.get('level','fine' if s['width']==1e-4 else 'coarse')+'_'+s['kind']+'_derivative'
        measure(label,error,2e-6 if s['kind']=='face' else 3e-6,s)
        checks.append(dict(**s,normalized_error=error,column_normalized_error=column_error))
    for p in rejects:
        if 'error' not in replies[p['index']]:failures.append(dict(check='missing_rejection',**p))
    sources=[Path(__file__),a.probe,family,table,profile,old_input,old_output]+[root/p for p in ['src/screened_microscopic_transport.cpp','include/ember/screened_microscopic_transport.hpp','src/eos_smooth_mixture.cpp','include/ember/eos_smooth_mixture.hpp','src/collision_transport.cpp','src/differential.hpp','src/material_flux.cpp','scripts/screened_microscopic_probe.cpp','src/eos_components.cpp']]
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed_hot_physical_faces_and_regular_chemical_parts' if not failures else 'failed',accepted_for_stellar_evolution=False,
        saved_faces=842,altered_faces=12,queries=len(queries),new_native_queries=len(missing),reused_queries=len(queries)-len(missing),derivative_checks=len(checks),regular_saved_states=512,regular_zero_states=12,rejections=rejects,
        maxima=maximum,failures=failures,unresolved_trace_columns=unresolved_trace_columns,checks=checks,faces=[dict(**f,response=replies[f['index']]) for f in faces],
        elapsed_seconds=elapsed,probe_timing=stderr,local_face_evaluation_seconds=sum(replies[f['index']].get('seconds',0.) for f in faces),
        input_sha256={str(p):sha(p) for p in sources},artifacts_sha256={str(p):sha(p) for p in a.scratch.iterdir()},
        limitations=['Conditional fully stripped screened collision model; minimum T=2e6 K is an explicit restriction, not an ionization proof.',
          'Both independent endpoint species must be positive. Exact-zero flux still requires a finite discretization; no floor or closed face is substituted.',
          'No stellar evolution, partial-ionization transport, or accepted screening prescription.',
          'Face T is geometric mean; radius, density, and composition are arithmetic means. Electron stiffness is relativistic ideal Fermi gas; collision electrons are nonrelativistic.',
          'Composition Jacobian comparisons use the row norm acting on relative abundance changes. Individually unresolved trace cross columns are retained explicitly; enriched-He3 controls separately require resolved column comparisons.',
          'Reported errors are numerical comparisons, not physical accuracy bounds; timings exclude EOS loading unless specified.'])
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in report.items() if k in ['outcome','queries','derivative_checks','maxima','failures','local_face_evaluation_seconds','probe_timing']}));
    if failures:raise SystemExit(1)

if __name__=='__main__':main()
