"""Test the archived gravity extension on the stopped mass's own surface."""
import hashlib
import json
import math
from pathlib import Path
import shlex
import shutil
import subprocess
from coordinate import last, out, read, reference, write

def members(path):
    paths=[Path(shlex.split(s)[1]) for s in path.read_text().splitlines() if s.startswith('file ')]
    def fnv(p):
        h=14695981039346656037
        for b in p.read_bytes():h=((h^b)*1099511628211)&0xffffffffffffffff
        return format(h,'016x')
    return {'atmosphere_pressure_matched':fnv(path),**{f'atmosphere.private_inputs.{i}':fnv(p) for i,p in enumerate(paths)}}

def prepare(case,m):
    cfg=(case/'config.txt').read_text()
    old=Path(dict(shlex.split(line) for line in cfg.splitlines() if line.strip())['atmosphere_pressure_matched'])
    archive=out/'fiducial-goal-oct22/joyplot_cool_gravity_candidate_v4'
    new=archive/'family.dat'
    if new==old:return {'eligible':False,'reason':'gravity extension already used'}
    r=last(case/'evolution/history.jsonl')
    x=r['surface_X'];z=r['surface_Z'];te=r['Teff_K'];g=math.log10(1.3271244e26*(m/100)/(r['radius_Rsun']*6.957e10)**2)
    coord=dict(Teff_K=te,logg=g,Z_over_X=z/x,dx=x-(1-.27783842244159906)/(1+z/x))
    review=read(archive/'independent_reviews.json')
    assert hashlib.sha256(new.read_bytes()).hexdigest()==review['family_sha256']
    checks=[v['coordinates'] for v in review['independent_evidence'] if v['matching_passes']]
    edge=read(reference/'atmosphere/edge.json')
    if read(reference/'atmosphere/edge-validation.json')['passes']:
        ratio=edge['Z_over_X']
        checks.append(dict(Teff_K=edge['teff'],logg=edge['logg'],Z_over_X=ratio,
                           dx=edge['arms']['anchor']['X']-(1-.27783842244159906)/(1+ratio)))
    def distance_to(point):
        return max(abs(coord[k]-point[k])/v for k,v in review['locality_limits'].items())
    distance=min(distance_to(p) for p in checks)
    info=dict(eligible=False,surface=coord,nearest_distance=distance,source=str(new))
    if any(not v['matching_passes'] and distance_to(v['coordinates'])<=1 for v in review['independent_evidence']):
        return dict(info,reason='nearby independent check failed')
    if distance>4 or not (1000<=te<=1300 and 5.49<=g<=5.65):
        return dict(info,reason='outside the archived gravity study')
    raw=(case/'evolution/final.checkpoint').read_bytes().splitlines(keepends=True)
    n=9+int(raw[8]);N=int(raw[n].split()[0]);zone=raw[n+N].split()
    assert abs(float(zone[7])-x)<1e-12 and raw[-1]==b'END\n'
    args=read(out/'highz-native-opacity-oct12-v1/matching_depth_controls.json')['command'][1:]
    line=f'thin {te} {g} {x} {float(zone[8])} {z}\n'
    def probe(table):
        q=subprocess.run([str(out/'boundary-corrected-item1-oct13/probe'),*args,str(table)],input=line,text=True,capture_output=True,timeout=40,check=True)
        return list(map(float,q.stdout.split()[6:9]))
    a,b=probe(old),probe(new)
    assert len(a)==len(b)==3
    diff=[y/x-1 for x,y in zip(a,b)]
    if not all(abs(e)<=t for e,t in zip(diff,[.003,.01,.02])):return dict(info,reason='surface overlap failed',relative=diff)
    if distance>1:
        from_native=case/'native-atmosphere'
        if not (from_native/'assessment.json').is_file():
            return dict(info,reason='one independent atmosphere column needed',prediction=b,helium3=float(zone[8]))
        validated=read(from_native/'assessment.json')
        if not validated.get('passes') or validated['prediction']!=b:
            return dict(info,reason='independent atmosphere column failed or prediction changed')
    ids=dict(shlex.split(l.decode()) for l in raw[9:n])
    for k,v in members(old).items():assert ids.pop(k)==v
    ids.update(members(new))
    folder=case/'gravity-control';folder.mkdir(exist_ok=True)
    ck=folder/'start.checkpoint'
    ck.write_bytes(b''.join(raw[:8])+f'{len(ids)}\n'.encode()+''.join(f'{json.dumps(k)} {json.dumps(v)}\n' for k,v in sorted(ids.items())).encode()+b''.join(raw[n:]))
    shutil.copyfile(case/'evolution/final.checkpoint.lithium',str(ck)+'.lithium')
    (folder/'config.txt').write_text(cfg.replace(str(old),str(new)))
    return dict(info,eligible=True,relative=diff,reason='local comparison passed; conditional gravity extension')

def launch_column(case,m,review):
    import fcntl,os
    d=case/'native-atmosphere';d.mkdir()
    s=read(reference/'atmosphere/edge.json')
    c=review['surface']
    x=(1-.27783842244159906)/(1+c['Z_over_X'])+c['dx'];z=c['Z_over_X']*x
    y3=review['helium3']
    s.update(tag=f'_solargrid_m{m:02}_{case.name}',teff=c['Teff_K'],logg=c['logg'],mass_Msun=m/100,
             years=last(case/'evolution/history.jsonl')['years'],Z_over_X=c['Z_over_X'],
             checkpoint=str(case/'evolution/final.checkpoint'),
             purpose='Independent column at stopped solar-grid surface',
             arms={'anchor':dict(X=x,Z=z,He3=y3,He4=1-x-z-y3)})
    write(d/'spec.json',s);write(d/'prediction.json',review['prediction'])
    shutil.copyfile(reference/'atmosphere/solve.py',d/'solve.py')
    exe=out/'picaso-pilot-oct11-fable/condaenv/bin/python'
    with (d/'worker.lock').open('a') as lock,(d/'run.log').open('w') as log:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        subprocess.Popen([str(exe),'solve.py','anchor','spec.json','3600'],cwd=d,stdout=log,stderr=log,
                         stdin=subprocess.DEVNULL,start_new_session=True,pass_fds=(lock.fileno(),),
                         env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))

def score_column(d):
    s=read(d/'spec.json')
    p=d/('hhe_arm_anchor'+s['tag']+'.json')
    if not p.is_file():return dict(passes=False,reason='column terminated without a final receipt',prediction=read(d/'prediction.json'))
    r=read(p);assert r['spec']==s
    prediction=read(d/'prediction.json');v=r['fresh']['at_tau']['3']
    delta=[a/b-1 for a,b in zip(prediction,[v['T'],v['P_bar']*1e6,v['depth_cm_mu_He3corr']])]
    return dict(passes=bool(r['quality']['usable'] and all(abs(e)<=lim for e,lim in zip(delta,[.003,.01,.02]))),
                prediction=prediction,relative=delta,receipt_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
