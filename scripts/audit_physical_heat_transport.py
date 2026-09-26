#!/usr/bin/env python3
"""Check physical heat on fixed-composition subspaces, including absent species.

Retained positive-state heat responses are reused only as independent reference
values. New thermal and stellar-equation derivatives are differenced directly.
The radiation law in the zone controls is analytic; no trajectory is advanced.
"""
import argparse,csv,gzip,hashlib,json,math,subprocess,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('probe',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--scratch',type=Path,required=True)
    p.add_argument('--reuse-report',type=Path);p.add_argument('--reuse-raw',type=Path)
    a=p.parse_args();assert not a.output.exists() and not a.scratch.exists()
    a.scratch.mkdir();root=Path(__file__).resolve().parents[1]
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    family=Path('/private/tmp/ember-smooth-eos-refined-family-v1/freeeos300_gs98_z020.dat')
    table=Path('/tmp/ember-collision-transport-table-v1.dat')
    profile=root/'docs/reports/2026-09-11/evolution_latest_profile.csv'
    with profile.open() as f:rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
    oldraw=Path('/tmp/ember-screened-microscopic-audit-v3')
    oldreport=root/'docs/results/screened_microscopic_transport_v3.json'
    old=json.loads(oldreport.read_text())
    for name,h in old['artifacts_sha256'].items():assert sha(name)==h,name
    saved={}
    for q,r in zip(gzip.decompress((oldraw/'queries.txt.gz').read_bytes()).decode().splitlines(),
                   gzip.decompress((oldraw/'responses.jsonl.gz').read_bytes()).decode().splitlines(),strict=True):
        saved.setdefault(q,json.loads(r))
    queries=[];metadata=[];stencils=[];comparisons=[];rejects=[]
    def point(r):return [math.log(r['radius_cm']),math.log(r['density_g_cm3']),math.log(r['temperature_K']),r['luminosity_erg_s'],r['X'],r['Y3']]
    def add(mode,q):
        i=len(queries);queries.append(mode+' '+' '.join(format(float(v),'.17g') for v in q));return i
    controls=[]
    for ions in [0,1]:
        for i in range(421):
            q=[ions,1,3,rows[i]['mass_g'],rows[i+1]['mass_g']]+point(rows[i])+point(rows[i+1])
            idx=add('heat',q)
            key='face '+' '.join(format(float(v),'.17g') for v in [ions,1,1]+q[3:])
            assert key in saved and 'error' not in saved[key],key
            comparisons.append(dict(index=idx,reference=saved[key],kind='saved_positive',face=i))
            if i in [160,396]:controls.append(dict(index=idx,query=q,mode='heat',kind='saved',face=i))
        for i in [0,160,300,360,396,420]:
            for kind in ['absent_H','absent_He3','pure_He']:
                q=[ions,1,3,rows[i]['mass_g'],rows[i+1]['mass_g']]+point(rows[i])+point(rows[i+1])
                for e in [0,1]:
                    if kind in ['absent_H','pure_He']:q[9+6*e]=0.
                    if kind in ['absent_He3','pure_He']:q[10+6*e]=0.
                idx=add('heat',q);metadata.append(dict(index=idx,query=q,kind=kind,face=i))
                mask=(1 if q[9]>0 else 0)+(2 if q[10]>0 else 0)
                reduced=q.copy();reduced[2]=mask;rid=add('face',reduced)
                comparisons.append(dict(index=idx,reference_index=rid,kind='explicit_reduced_subspace',face=i))
                scalar=q.copy();scalar[1]=0;sid=add('heat',scalar)
                comparisons.append(dict(index=idx,reference_index=sid,kind='scalar',face=i))
                reversed_q=q.copy();reversed_q[5:11],reversed_q[11:17]=q[11:17],q[5:11]
                reverse=add('heat',reversed_q)
                comparisons.append(dict(index=idx,reference_index=reverse,kind='reversed',face=i))
                if i in [160,300]:controls.append(dict(index=idx,query=q,mode='heat',kind=kind,face=i))
                if i in [396,420] and kind in ['absent_H','pure_He']:
                    for j in [idx,rid,sid,reverse]:rejects.append(dict(index=j,kind='unsupported_zero_subspace',expected='electron density outside pair-table degeneracy range'))
        for i in [160,300]:
            for kind in ['saved','absent_He3','pure_He']:
                q=[ions,1,3,rows[i]['mass_g'],rows[i+1]['mass_g']]+point(rows[i])+point(rows[i+1])
                if kind!='saved':
                    q[10]=q[16]=0.
                    if kind=='pure_He':q[9]=q[15]=0.
                idx=add('zone',q);controls.append(dict(index=idx,query=q,mode='zone',kind=kind,face=i))
    for c in controls:
        for e in [0,1]:
            for k in range(4):
                q=c['query'];coord=5+6*e+k
                unit=1. if k<3 else max(abs(q[coord]),1e25)
                for width in [2e-4,1e-4]:
                    ids=[]
                    for step in [-2,-1,1,2]:
                        trial=q.copy();trial[1]=0;trial[coord]+=step*unit*width
                        ids.append(add(c['mode'],trial))
                    stencils.append(dict(index=c['index'],mode=c['mode'],endpoint=e,coordinate=k,unit=unit,width=width,replies=ids))
    for kind in ['one_zero_H','one_zero_He3','removed_present','cold','mass','negative']:
        q=controls[0]['query'].copy()
        if kind=='one_zero_H':q[9]=0
        if kind=='one_zero_He3':q[10]=0
        if kind=='removed_present':q[2]=0
        if kind=='cold':q[7]=math.log(1e6)
        if kind=='mass':q[4]=q[3]
        if kind=='negative':q[10]=-.001
        rejects.append(dict(index=add('heat',q),kind=kind))
    # The default species interface must still reject both-endpoint absence.
    q=controls[0]['query'].copy();q[10]=q[16]=0
    rejects.append(dict(index=add('face',q),kind='strict_species_domain'))
    inputs=[Path(__file__),a.probe,family,table,profile,oldreport,
      root/'include/ember/microscopic_transport.hpp',root/'src/microscopic_transport.cpp',
      root/'include/ember/screened_microscopic_transport.hpp',root/'src/screened_microscopic_transport.cpp',
      root/'src/evolution.cpp',root/'src/mixing.cpp',root/'src/structure.cpp',
      root/'src/eos_smooth_mixture.cpp',root/'src/collision_transport.cpp',root/'scripts/physical_heat_probe.cpp']
    identities={str(p):sha(p) for p in inputs}
    assert identities[str(family)]=='7d6bfa0b26264834307dad651852bc5ff693f91f3861b4934e582dee047fd7f5'
    assert identities[str(table)]=='48f7328da79a6163a997ab524c623ecb334d964f47bb7fbdefc084e63b46e90f'
    payload='\n'.join(queries)+'\n'
    (a.scratch/'queries.txt.gz').write_bytes(gzip.compress(payload.encode(),mtime=0))
    reused={}
    if a.reuse_report:
        assert a.reuse_raw
        prior=json.loads(a.reuse_report.read_text())
        for name,h in prior['input_sha256'].items():
            if Path(name)!=Path(__file__):assert sha(name)==h,name
        for name,h in prior['artifacts_sha256'].items():assert sha(name)==h,name
        qq=gzip.decompress((a.reuse_raw/'queries.txt.gz').read_bytes()).decode().splitlines()
        rr=gzip.decompress((a.reuse_raw/'responses.jsonl.gz').read_bytes()).decode().splitlines()
        for q,r in zip(qq,rr,strict=True):reused.setdefault(q,r)
    missing=list(dict.fromkeys(q for q in queries if q not in reused))
    start=time.monotonic();stderr='all responses reused\n'
    if missing:
        missing_payload='\n'.join(missing)+'\n'
        assert missing_payload.count('\n')==len(missing)
        (a.scratch/'fresh_queries.txt.gz').write_bytes(gzip.compress(missing_payload.encode(),mtime=0))
        run=subprocess.run([str(a.probe),str(family),str(table)],input=missing_payload,text=True,capture_output=True,timeout=160,check=True)
        (a.scratch/'fresh_responses.jsonl.gz').write_bytes(gzip.compress(run.stdout.encode(),mtime=0))
        (a.scratch/'stderr.txt').write_text(run.stderr)
        rr=run.stdout.splitlines();assert len(rr)==len(missing)
        for q,r in zip(missing,rr,strict=True):reused[q]=r
        stderr=run.stderr
    elapsed=time.monotonic()-start
    stdout='\n'.join(reused[q] for q in queries)+'\n'
    (a.scratch/'responses.jsonl.gz').write_bytes(gzip.compress(stdout.encode(),mtime=0))
    (a.scratch/'stderr.txt').write_text(stderr)
    replies=[json.loads(x) for x in stdout.splitlines()]
    assert len(replies)==len(queries)
    failures=[];checks=[];maxima={}
    def measure(kind,error,tolerance,context):
        maxima[kind]=max(maxima.get(kind,0.),float(error))
        item=dict(kind=kind,error=float(error),tolerance=tolerance,context=context)
        checks.append(item)
        if not math.isfinite(error) or error>tolerance:failures.append(item)
    for i,r in enumerate(replies):
        if 'error' in r and i not in {x['index'] for x in rejects}:failures.append(dict(index=i,error=r['error']))
    def heatvec(r,derivatives=True):
        v=[r['carried'],r['conductivity']]
        if derivatives:
            for key in ['dcarried_lo','dcarried_hi','dconductivity_lo','dconductivity_hi']:v+=r[key]
        return np.array(v)
    for c in comparisons:
        r=replies[c['index']];ref=c.get('reference',replies[c['reference_index']] if 'reference_index' in c else None)
        if 'error' in r or 'error' in ref:continue
        derivative=c['kind'] not in ['scalar','reversed']
        v,w=heatvec(r,derivative),heatvec(ref,derivative)
        if c['kind']=='reversed':w[0]*=-1
        measure(c['kind'],float(np.max(abs(v-w)/np.maximum(abs(w),1.))),2e-11,c['index'])
    for c in controls:
        r=replies[c['index']]
        if 'error' in r:continue
        if c['mode']=='heat':
            if r['conductivity']<=0:failures.append(dict(index=c['index'],error='nonpositive conductivity'))
            for e in ['lo','hi']:
                if r['dcarried_'+e][3]!=0 or r['dconductivity_'+e][3]!=0:
                    failures.append(dict(index=c['index'],error='heat depends on luminosity'))
        elif not r['nonthermal_equations_unchanged']:
            failures.append(dict(index=c['index'],error='microscopic heat changed mass, hydrostatic or energy equation'))
    centers={c['index']:c for c in controls}
    for s in stencils:
        r=replies[s['index']];rr=[replies[i] for i in s['replies']]
        if 'error' in r or any('error' in z for z in rr):continue
        q=centers[s['index']]['query'];side='lo' if s['endpoint']==0 else 'hi';k=s['coordinate']
        units=np.array([1,1,1,max(abs(q[8]),1e25),1,1,1,max(abs(q[14]),1e25)])
        if s['mode']=='heat':
            jac=np.array([r['dcarried_lo']+r['dcarried_hi'],r['dconductivity_lo']+r['dconductivity_hi']])*units
            values=[np.array([z['carried'],z['conductivity']]) for z in rr]
        else:
            jac=np.c_[r['left'],r['right']]*units
            values=[np.array(z['residual']) for z in rr]
        expected=(values[0]-8*values[1]+8*values[2]-values[3])/(12*s['width'])
        actual=jac[:,4*s['endpoint']+k]
        scale=np.maximum(np.max(abs(jac),axis=1),1e-200)
        error=float(np.max(abs(actual-expected)/scale))
        level='fine' if s['width']==1e-4 else 'coarse'
        measure(s['mode']+'_'+level+'_jacobian',error,3e-6 if s['mode']=='heat' else 1e-5,s)
    for c in rejects:
        if 'error' not in replies[c['index']]:failures.append(dict(index=c['index'],error='invalid input accepted',kind=c['kind']))
        elif c.get('expected') and c['expected'] not in replies[c['index']]['error']:failures.append(dict(index=c['index'],error='wrong domain rejection',kind=c['kind']))
    for name,h in identities.items():assert sha(name)==h
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),
      outcome='passed_physical_thermal_subspaces_and_stellar_equations' if not failures else 'failed',
      accepted_for_stellar_evolution=False,queries=len(queries),positive_saved_faces=842,
      zero_subspace_faces=len(metadata)-8,unsupported_zero_subspace_faces=8,derivative_centers=len(controls),derivative_comparisons=len(stencils),
      rejection_controls=len(rejects),rejections=rejects,maxima=maxima,failures=failures,checks=checks,
      zero_states=metadata,elapsed_seconds=elapsed,native_timing=stderr,new_native_queries=len(missing),reused_native_queries=len(queries)-len(missing),
      input_sha256=identities,artifacts_sha256={str(p):sha(p) for p in a.scratch.iterdir()},
      limitations=['Conditional fully stripped collision model, T>=2e6 K.',
        'Thermal derivatives hold composition fixed; no derivative introducing an absent species is claimed.',
        'A species absent at just one evaluated endpoint remains rejected.',
        'Zone controls use analytic radiation and fixed states; no new stellar age or trajectory.'])
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ['outcome','queries','zero_subspace_faces','derivative_comparisons','rejection_controls','maxima','failures','elapsed_seconds']},indent=2))
    if failures:raise SystemExit(1)
if __name__=='__main__':main()
