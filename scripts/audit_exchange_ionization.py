#!/usr/bin/env python3
"""Direct retained-FreeEOS checks of ionization along the saved stellar profile.

Use the existing source probes and isotope/mass conversion. No source grid,
ionization formula or stellar input is changed. Electron quadratures are
compared independently; source ion fractions are not transport coefficients.
"""
import argparse,csv,hashlib,json,math,subprocess,time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
from pathlib import Path
from metal_eos_composition import mixture


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path)
    p.add_argument('--scratch',type=Path,required=True);a=p.parse_args()
    assert not a.output.exists() and not a.scratch.exists();a.scratch.mkdir()
    root=Path(__file__).resolve().parents[1]
    profile=root/'docs/reports/2026-09-11/evolution_latest_profile.csv'
    assert sha(profile)=='28ab0b571dbc495ff2f48a6c6e7cd79f09e5711f3b725a2f478d4b0390d7ef2f'
    with profile.open() as f:rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
    states=[]
    for i in [0,100,200,300,350,380,395,396,397,410,420,430,450,470,490,511]:
        row=rows[i];states.append(dict(kind='node',index=i,T=row['temperature_K'],rho=row['density_g_cm3'],X=row['X'],Y3=row['Y3']))
    for i in [396,420]:
        left,right=rows[i:i+2]
        states.append(dict(kind='face',index=i,T=math.sqrt(left['temperature_K']*right['temperature_K']),
                           rho=.5*(left['density_g_cm3']+right['density_g_cm3']),X=.5*(left['X']+right['X']),Y3=.5*(left['Y3']+right['Y3'])))
    base=next(s for s in states if s['kind']=='face' and s['index']==396)
    for coordinate in ['T','rho']:
        for factor in [.8,1.2]:states.append({**base,'kind':'altered_base','coordinate':coordinate,'factor':factor,coordinate:base[coordinate]*factor})
    builds=[]
    for directory in ['/tmp/ember-freeeos-precision-builder-check-v1','/tmp/ember-freeeos-quadrature13-builder-check-v1']:
        work=Path(directory);r=json.loads((work/'build_receipt.json').read_text())
        assert sha(work/'probe')==r['probe_sha256'] and sha(r['library'])==r['library_sha256']
        builds.append(dict(probe=str(work/'probe'),receipt=r,receipt_sha256=sha(work/'build_receipt.json')))
    identities={str(p):sha(p) for p in [Path(__file__),profile,root/'scripts/metal_eos_composition.py',root/'scripts/freeeos_probe.f90',
                                       root/'scripts/generate_nongrey_grid.py',root/'scripts/stellar_composition.py']}
    spec=dict(states=states,builds=builds,input_sha256=identities,source_options=[3,223,-2])
    (a.scratch/'specification.json').write_text(json.dumps(spec,indent=2)+'\n')
    groups={}
    for i,s in enumerate(states):groups.setdefault((s['X'],s['Y3']),[]).append(i)
    jobs=[(bi,x,y,indices) for bi in range(2) for (x,y),indices in groups.items()]
    def calculate(item):
        ji,(bi,x,y,indices)=item;mix=mixture(x,y);scale=mix['source_mass_scale']
        request=' '.join(format(v,'.17g') for v in mix['eps'])+'\n3 223 -2\n'
        request+=''.join(f"{math.log(states[i]['rho']*scale):.17g} {math.log(states[i]['T']):.17g}\n" for i in indices)
        prefix=a.scratch/f'group-{ji:02d}'
        prefix.with_suffix('.in').write_text(request);start=time.monotonic()
        try:
            r=subprocess.run([builds[bi]['probe']],input=request,text=True,capture_output=True,timeout=30)
            out,err,code=r.stdout,r.stderr,r.returncode
        except subprocess.TimeoutExpired as e:
            out=e.stdout or b'';err=e.stderr or b'';code=124
            if isinstance(out,bytes):out=out.decode()
            if isinstance(err,bytes):err=err.decode()
        prefix.with_suffix('.out').write_text(out);prefix.with_suffix('.err').write_text(err)
        result=dict(group=ji,build=bi,indices=indices,mixture=mix,elapsed_seconds=time.monotonic()-start,exit_code=code,
                    stdout=str(prefix.with_suffix('.out')),stderr=str(prefix.with_suffix('.err')),source_rows=[])
        for line in out.splitlines():result['source_rows'].append(list(map(float,line.split())))
        prefix.with_suffix('.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');return result
    begin=time.monotonic()
    with ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(calculate,enumerate(jobs)))
    elapsed=time.monotonic()-begin;failures=[];solutions={}
    charges=[1,2,6,7,8,10,11,12,13,14,15,16,17,18,20,22,24,25,26,28]
    for r in results:
        if r['exit_code'] or len(r['source_rows'])!=len(r['indices']):failures.append(dict(group=r['group'],failure='process or row count'));continue
        mix=r['mixture'];scale=mix['source_mass_scale']
        for i,row in zip(r['indices'],r['source_rows']):
            if len(row)!=22 or row[0]!=0 or not all(math.isfinite(v) for v in row):
                failures.append(dict(group=r['group'],state=i,failure='invalid source row',row=row));continue
            h,he1,he2=row[19:22]
            if min(h,he1,he2)<0 or max(h,he1,he2)>1 or he1+he2>1+1e-12:failures.append(dict(group=r['group'],state=i,failure='ion fraction range'))
            bare=scale*sum(e*z for e,z in zip(mix['eps'],charges))
            electrons=row[17]*scale
            solutions[(r['build'],i)]=dict(**states[i],H_ionized=h,He_singly_ionized=he1,He_doubly_ionized=he2,
                He_neutral=1-he1-he2,free_electrons_per_baryonic_mass=electrons,
                fully_stripped_free_electrons=bare,relative_electron_deficit=1-electrons/bare,
                source_iterations=int(row[1]))
    comparisons=[]
    for i in range(len(states)):
        if (0,i) not in solutions or (1,i) not in solutions:continue
        left,right=solutions[0,i],solutions[1,i]
        difference=max(abs(left[k]-right[k]) for k in ['H_ionized','He_singly_ionized','He_doubly_ionized','free_electrons_per_baryonic_mass'])
        comparisons.append(dict(state=i,maximum_absolute_difference=difference))
        if difference>1e-8:failures.append(dict(state=i,failure='source quadrature difference',difference=difference))
    assert all(sha(name)==value for name,value in identities.items())
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='completed_native_ionization_comparison' if not failures else 'failed',
                accepted_for_stellar_transport=False,new_source_queries=2*len(states),source_groups=len(jobs),elapsed_seconds=elapsed,
                failures=failures,quadrature_comparisons=comparisons,states=[solutions[1,i] for i in range(len(states)) if (1,i) in solutions],
                input_sha256=identities,artifacts_sha256={str(p):sha(p) for p in a.scratch.iterdir()},
                limitations=['FreeEOS pressure-ionization and atomic models, not an independent atomic calculation.',
                             'FreeEOS omits potassium; source He4 electronic physics represents both He isotopes.',
                             'Mean free-electron deficit does not determine the metals mean-square charge or a collision matrix.',
                             'No transport coefficients, energy closure or trajectory are selected by an ion-fraction query.'])
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ['outcome','elapsed_seconds','new_source_queries','failures']},indent=2))
    print(json.dumps([s for s in report['states'] if s['kind']=='face'],indent=2))
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
