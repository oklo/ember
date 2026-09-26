#!/usr/bin/env python3
"""Matched fully convective heat comparisons with the verified ion extension."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,subprocess,time,shutil
import numpy as np
from evolve_cn_transport import TABLES,audit,compare
from audit_cn_thermal import C0,N0,YEAR,physical


def main():
    scratch=Path('/tmp/ember-cn-convective-heat-v1');scratch.mkdir()
    output=Path('docs/results/cn_convective_heat_v1.json')
    if output.exists():raise ValueError('preserve completed comparisons')
    native=Path('/tmp/ember-cn-microscopic-v5/build/tests/cn_thermal_probe')
    collision=Path('docs/research/fable/results/fable-ion-table-extension-v1/collision_transport_table_candidate_v2.dat').resolve()
    acceptance=json.loads(Path('docs/results/fable_ion_extension_review_v1.json').read_text())
    assert acceptance['outcome']=='passed_independent_source_export_and_native_controls'
    assert hashlib.sha256(collision.read_bytes()).hexdigest()==acceptance['sha256'][str(collision.relative_to(Path.cwd()))]
    hashes={}
    for p in [Path(__file__),Path('scripts/evolve_cn_transport.py'),Path('scripts/audit_cn_thermal.py'),native,collision,Path('docs/results/fable_ion_extension_review_v1.json'),*[Path(t) for t in TABLES]]:
        hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
        if p.suffix=='.py':shutil.copy2(p,scratch/p.name)
    start=time.monotonic();stderr=(scratch/'native.stderr').open('w');command=[str(native),*TABLES,'.1','1',str(collision)]
    process=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=stderr,text=True)
    (scratch/'native_pid.txt').write_text(str(process.pid)+'\n');cases=[];failure=None
    def step(old,dt,mode,label):
        values=[1,1.,1.,dt,1e-13,len(old),mode,1e-6,*old.flat]
        query=' '.join(format(x,'.17g') for x in values)+'\n'
        (scratch/(label+'.query')).write_text(query);process.stdin.write(query);process.stdin.flush();raw=process.stdout.readline()
        (scratch/(label+'.json')).write_text(raw);result=json.loads(raw)
        model,checks,_=audit(old,result,dt)
        return model,result,checks
    try:
        for label,path,years in [('2500Tyr','out/evolution-metal-512-2500gyr-gas-checkpoint.json',1e7),('3300Tyr','/tmp/ember-cn-early-continuation-v4/checkpoint.json',1e7)]:
            p=Path(path);hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();d=json.loads(p.read_text())
            if 'profile' in d:
                m=np.array(d['profile'])[:,:7];old=np.column_stack((m,np.tile([C0,0,N0],(len(m),1))))
            else:old=np.array(d['model_record']['model'])
            endpoints=[];summaries=[]
            for mode in [2,1]:
                _,full,c1=step(old,years*YEAR,mode,f'{label}-{mode}-full')
                half,_,c2=step(old,years*YEAR/2,mode,f'{label}-{mode}-half1')
                fine,last,c3=step(half,years*YEAR/2,mode,f'{label}-{mode}-half2')
                norm,metrics=compare(full,last)
                if norm>1:raise ValueError(f'timestep control {label}, mode{mode}: {metrics}')
                endpoints.append(last);summaries.append(dict(mode=mode,error_norm=norm,criteria=metrics,conservation=[c1,c2,c3]))
            before,after=map(lambda r:np.array(r['model']),endpoints)
            response=dict(relative_L=float(after[-1,4]/before[-1,4]-1),relative_R=float(after[-1,1]/before[-1,1]-1),relative_Teff=float(endpoints[1]['Teff']/endpoints[0]['Teff']-1),maximum_log_structure=float(abs(np.log(after[:,1:4]/before[:,1:4])).max()),maximum_physical_abundance=float(abs(physical(after)-physical(before)).max()),same_mixing_regions=endpoints[0]['mixing_regions']==endpoints[1]['mixing_regions'])
            cases.append(dict(case=label,interval_years=years,response=response,controls=summaries))
            print(label,response,flush=True)
    except Exception as e:failure=str(e);print('FAILED',failure,flush=True)
    finally:process.stdin.close();code=process.wait(timeout=120);stderr.close()
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed' if not failure and code==0 else 'failed',failure=failure,wall_seconds=time.monotonic()-start,cases=cases,native_command=command,input_sha256=hashes,limitations=['Matched 10 Myr response at two fully convective structures; no full-interval bound implied.','Same fully stripped collision and fixed-GS98 material approximations; the extended table changes numerical domain only.','Material mode retains conduction and total composition enthalpy, with reduced microscopic heat omitted; comparison mode includes the hot kinetic prescription and its 2–3 MK join.'])
    output.write_text(json.dumps(report,indent=2)+'\n')
    return int(bool(failure) or code!=0)


if __name__=='__main__':raise SystemExit(main())
