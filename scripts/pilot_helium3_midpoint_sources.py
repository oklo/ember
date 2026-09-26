#!/usr/bin/env python3
"""Measure coverage, consistency and cost of the required Y3=.06 source plane.

The specification distinguishes the complete target domain from this finite
allocation of isotherms. All finite FreeEOS flags are preserved. A pilot is
not a complete table or acceptance for evolution. Outputs use baryonic units.
"""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

from metal_eos_composition import mixture


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root=Path(__file__).resolve().parents[1]
    plan=root/'data/eos/sources/helium3_midpoint_pilot_v1_specification.json'
    spec=json.loads(plan.read_text())
    work=Path('/tmp/ember-helium3-midpoint-pilot-run-v1')
    source=spec['source'];probe=Path(source['probe'])
    assert digest(probe)==source['probe_sha256']
    assert digest(source['library'])==source['library_sha256']
    for p,h in spec['inputs_sha256'].items():assert digest(p)==h
    start=time.monotonic()
    def run(item):
        idx=item['request_index'];x=item['hydrogen'];y=spec['helium3']
        t=item['logT'];qs=item['logQ'];mix=mixture(x,y);scale=mix['source_mass_scale']
        request=' '.join(format(v,'.17g') for v in mix['eps'])+'\n3 223 -2\n'
        request+=''.join(f'{math.log(scale)+math.log(10)*(q+1.5*(t-6)):.17g} {math.log(10)*t:.17g}\n' for q in qs)
        beginning=time.monotonic()
        p=subprocess.run([str(probe)],input=request,text=True,capture_output=True,timeout=90)
        try:data=[list(map(float,l.split())) for l in p.stdout.splitlines()]
        except ValueError:data=[]
        basic=p.returncode==0 and len(data)==len(qs) and all(len(r)==22 and all(map(math.isfinite,r)) for r in data)
        if not basic:
            failure={'request':item,'input':request,'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr}
            (work/f'request-{idx:03d}-failed.json').write_text(json.dumps(failure,indent=2)+'\n')
            return {'request_index':idx,'outcome':'source_process_or_output_failed','returncode':p.returncode}
        max_defect=0.;consistency=[];flagged=[];nonpositive=[]
        for i,row in enumerate(data):
            row[2]/=scale
            for k in [5,6,9,10,11,12,13,16,17]:row[k]*=scale
            if row[0]!=0:flagged.append({'density_index':i,'info':row[0]});continue
            rho,T,P=row[2:5];chit,Er,Et,Sr,St=row[8:13]
            if P<=0 or Et==0:
                nonpositive.append(i);continue
            defect=max(abs(rho*Er/P+chit-1),abs(T*St/Et-1),abs(rho*T*Sr/P+chit))
            max_defect=max(max_defect,defect)
            if defect>1e-7:consistency.append({'density_index':i,'maximum_absolute_defect':defect})
        record={'source_options':[3,223,-2],'actual_probe_sha256':source['probe_sha256'],
                'relative_electron_integral_target':source['relative_electron_integral_target'],
                'input':request,'input_sha256':hashlib.sha256(request.encode()).hexdigest(),
                'returncode':p.returncode,'stderr':p.stderr,'request':item,'mixture':mix,'data':data}
        path=work/f'request-{idx:03d}.json.gz'
        path.write_bytes(gzip.compress((json.dumps(record,separators=(',',':'),allow_nan=False)+'\n').encode(),mtime=0))
        out={'request_index':idx,'hydrogen':x,'helium3':y,'logT':t,'states':len(data),
             'outcome':'source_saved','elapsed_seconds':time.monotonic()-beginning,
             'source_file':str(path),'source_sha256':digest(path),'bytes':path.stat().st_size,
             'nonconvergence_flags':flagged,'first_law_failures':consistency,
             'nonpositive_pressure_or_zero_heat_response':nonpositive,'maximum_first_law_defect':max_defect}
        print(json.dumps({k:out[k] for k in ['request_index','states','elapsed_seconds','maximum_first_law_defect']}),flush=True)
        return out
    records=[]
    with (work/'completed.jsonl').open('x') as f:
        with ThreadPoolExecutor(max_workers=4) as pool:
            for result in pool.map(run,spec['pilot_requests']):
                records.append(result);f.write(json.dumps(result,allow_nan=False)+'\n');f.flush()
    saved=[r for r in records if r['outcome']=='source_saved']
    report={'created_utc':datetime.now(timezone.utc).isoformat(),'scope':__doc__,
            'outcome':'completed_source_pilot','accepted_for_evolution':False,
            'specification':str(plan),'specification_sha256':digest(plan),
            'full_target_source_states':spec['full_target_source_states'],
            'pilot_requests':len(records),'saved_requests':len(saved),
            'saved_states':sum(r['states'] for r in saved),'retained_bytes':sum(r['bytes'] for r in saved),
            'elapsed_seconds':time.monotonic()-start,
            'source_nonconvergence_states':sum(len(r['nonconvergence_flags']) for r in saved),
            'source_first_law_failures':sum(len(r['first_law_failures']) for r in saved),
            'source_nonpositive_pressure_or_zero_heat_response':sum(len(r['nonpositive_pressure_or_zero_heat_response']) for r in saved),
            'maximum_first_law_defect':max((r['maximum_first_law_defect'] for r in saved),default=None),
            'requests':records,'source':source,
            'inputs_sha256':{str(p):digest(p) for p in [Path(__file__),root/'scripts/metal_eos_composition.py',plan]}}
    output=root/'docs/results/helium3_midpoint_source_pilot_v1.json'
    with output.open('x') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:report[k] for k in ['saved_requests','saved_states','retained_bytes',
        'source_nonconvergence_states','source_first_law_failures','maximum_first_law_defect']}))


if __name__=='__main__':main()
