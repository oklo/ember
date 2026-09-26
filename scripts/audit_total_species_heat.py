#!/usr/bin/env python3
"""Native controls for material enthalpy at prescribed total species rates."""
import argparse, csv, gzip, hashlib, json, math, subprocess, time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('probe',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--scratch',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists() and not a.scratch.exists();a.scratch.mkdir()
    root=Path(__file__).resolve().parents[1]
    profile=root/'docs/reports/2026-09-11/evolution_latest_profile.csv'
    eos=Path('/tmp/ember-smooth-eos-refined-family-v1/freeeos300_gs98_z020.dat')
    collision=Path('/tmp/ember-collision-transport-table-v1.dat')
    prior=json.loads(Path('/tmp/ember-full-star-connection-v2/inputs.json').read_text())['input_sha256']
    identities={k:v for k,v in prior.items() if '/eos/' in k or 'eos-refined-family' in k}
    for name,h in identities.items():assert sha(name)==h,name
    sources=[Path(__file__),a.probe,profile,eos,collision,
             root/'scripts/total_species_heat_probe.cpp',root/'src/screened_microscopic_transport.cpp',
             root/'src/microscopic_transport.cpp',root/'include/ember/microscopic_transport.hpp',
             root/'include/ember/screened_microscopic_transport.hpp',
             Path('/tmp/ember-total-heat-build-v1/src/libember.a')]
    identities.update({str(p.resolve()):sha(p) for p in sources})
    assert sha(profile)=='28ab0b571dbc495ff2f48a6c6e7cd79f09e5711f3b725a2f478d4b0390d7ef2f'
    assert sha(collision)=='48f7328da79a6163a997ab524c623ecb334d964f47bb7fbdefc084e63b46e90f'
    with profile.open() as f:rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
    queries=[];centers=[];stencils=[];rejects=[]
    def point(r):return [math.log(r['radius_cm']),math.log(r['density_g_cm3']),math.log(r['temperature_K']),r['luminosity_erg_s'],r['X'],r['Y3']]
    def add(q):
        index=len(queries);queries.append('total '+' '.join(format(float(v),'.17g') for v in q));return index
    for ions in [0,1]:
        for face in [0,160,300,396,410,420]:
            q=[ions,1,3,rows[face]['mass_g'],rows[face+1]['mass_g']]+point(rows[face])+point(rows[face+1])
            scale=rows[face]['mass_g']/rows[-1]['mass_g'];q += [5e12*scale,-1e11*scale]
            centers.append(dict(query=q,index=add(q),face=face,kind='saved'))
    for mask in [0,1,2]:
        q=centers[7]['query'].copy();q[2]=mask
        for k in [0,1]:
            if not mask&(1<<k):q[9+k]=q[15+k]=q[17+k]=0
        centers.append(dict(query=q,index=add(q),face=160,kind='absent_species',mask=mask))
    for control in centers:
        q=control['query'];refs={}
        for kind in ['scalar','zero','double','reverse']:
            trial=q.copy();trial[1]=0
            if kind=='zero':trial[17:19]=[0.,0.]
            if kind=='double':trial[17:19]=[2*v for v in q[17:19]]
            if kind=='reverse':
                trial[5:11],trial[11:17]=q[11:17],q[5:11];trial[17:19]=[-v for v in q[17:19]]
            refs[kind]=add(trial)
        control['references']=refs
        for endpoint in [0,1]:
            for coordinate in range(3):
                for width in [2e-4,1e-4]:
                    indices=[]
                    for step in [-2,-1,1,2]:
                        trial=q.copy();trial[1]=0;trial[5+6*endpoint+coordinate]+=step*width
                        indices.append(add(trial))
                    stencils.append(dict(center=control['index'],endpoint=endpoint,coordinate=coordinate,width=width,replies=indices))
    for kind in ['removed_flux','cold','mass','one_absent']:
        q=centers[-2]['query'].copy()
        if kind=='removed_flux':q[18]=1e10
        if kind=='cold':q[7]=math.log(1e6)
        if kind=='mass':q[4]=q[3]
        if kind=='one_absent':q[9]=0
        rejects.append(dict(index=add(q),kind=kind))
    payload='\n'.join(queries)+'\n'
    (a.scratch/'queries.txt.gz').write_bytes(gzip.compress(payload.encode(),mtime=0))
    command=[str(a.probe),str(eos),str(collision)]
    (a.scratch/'inputs.json').write_text(json.dumps(dict(input_sha256=identities,command=command,centers=centers,stencils=stencils,rejects=rejects),indent=2)+'\n')
    begin=time.monotonic();run=subprocess.run(command,input=payload,text=True,capture_output=True,timeout=240)
    elapsed=time.monotonic()-begin
    (a.scratch/'responses.jsonl.gz').write_bytes(gzip.compress(run.stdout.encode(),mtime=0))
    (a.scratch/'stderr.txt').write_text(run.stderr)
    replies=[json.loads(line) for line in run.stdout.splitlines()]
    failures=[];maxima={};checks=0
    if run.returncode or len(replies)!=len(queries):failures.append(dict(exit_code=run.returncode,replies=len(replies),expected=len(queries)))
    def check(kind,value,tolerance,context):
        nonlocal checks
        checks+=1;value=float(value);maxima[kind]=max(maxima.get(kind,0),value)
        if not math.isfinite(value) or value>tolerance:failures.append(dict(kind=kind,value=value,tolerance=tolerance,context=context))
    rejection_indices={r['index'] for r in rejects}
    for i,r in enumerate(replies):
        if 'error' in r and i not in rejection_indices:failures.append(dict(index=i,error=r['error']))
    for c in centers:
        r=replies[c['index']]
        if 'error' in r:continue
        h=np.array(r['total_rate_enthalpy']);flux=np.array(c['query'][17:19]);micro=np.array(r['microscopic_rate'])
        expected=r['microscopic_carried']+h@(flux-micro)
        scale=max(abs(r['carried']),abs(r['microscopic_carried']),abs(h@flux),1.)
        check('independent_heat_split',abs(r['carried']-expected)/scale,2e-11,c['index'])
        ref={k:replies[v] for k,v in c['references'].items()}
        if any('error' in v for v in ref.values()):continue
        check('rate_linearity',abs(ref['double']['carried']-ref['zero']['carried']-2*h@flux)/scale,2e-11,c['index'])
        check('scalar_identity',abs(ref['scalar']['carried']-r['carried'])/scale,2e-11,c['index'])
        check('reversal',abs(ref['reverse']['carried']+r['carried'])/scale,2e-11,c['index'])
        for v in ref.values():check('conductivity_independent_of_rate',abs(v['conductivity']/r['conductivity']-1),2e-11,c['index'])
        for side in ['lo','hi']:
            if r['dcarried_'+side][3]!=0 or r['dconductivity_'+side][3]!=0:failures.append(dict(index=c['index'],kind='luminosity_dependence'))
    for s in stencils:
        r=replies[s['center']];rr=[replies[i] for i in s['replies']]
        if 'error' in r or any('error' in v for v in rr):continue
        matrix=np.array([r['dcarried_lo']+r['dcarried_hi'],r['dconductivity_lo']+r['dconductivity_hi']])
        values=[np.array([v['carried'],v['conductivity']]) for v in rr]
        numerical=(values[0]-8*values[1]+8*values[2]-values[3])/(12*s['width'])
        expected=matrix[:,4*s['endpoint']+s['coordinate']]
        scale=np.maximum(np.max(abs(matrix),axis=1),1e-200)
        check('thermal_jacobian_'+str(s['width']),np.max(abs(numerical-expected)/scale),3e-6,s)
    for r in rejects:
        if 'error' not in replies[r['index']]:failures.append(dict(**r,error='invalid input accepted'))
    assert all(sha(name)==value for name,value in identities.items())
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed_native_total_species_heat_controls' if not failures else 'failed',
                accepted_for_stellar_evolution=False,queries=len(queries),centers=len(centers),derivative_comparisons=len(stencils),checks=checks,
                maxima=maxima,failures=failures,elapsed_seconds=elapsed,native_timing=run.stderr,input_sha256=identities,
                artifacts_sha256={str(p):sha(p) for p in a.scratch.iterdir()},
                limitations=['Conditional hot collision model with T>=2e6 K.',
                             'Supplied total species rates are held fixed in these thermal derivatives.',
                             'No cooler transport, moving convection boundary or trajectory is validated here.'])
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ['outcome','checks','maxima','failures','elapsed_seconds']},indent=2))
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
