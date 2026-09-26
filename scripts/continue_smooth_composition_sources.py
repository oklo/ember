#!/usr/bin/env python3
"""Continue a bounded H/He3 source collection, reusing all saved isotherms.

Every source request is a local FreeEOS isotherm. Existing complete composition
planes are excluded explicitly; the off-grid pilot is reused for validation. This collector preserves source flags and first-law failures;
it does not assemble or accept an EOS. The allocation has fixed wall, CPU
and storage bounds and never changes an existing result or allocation clock.
"""
from collections import deque
from concurrent.futures import ThreadPoolExecutor,wait,FIRST_COMPLETED
from datetime import datetime,timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time

from metal_eos_composition import mixture


def now():return datetime.now(timezone.utc).isoformat()
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def encode(x):return (json.dumps(x,separators=(',',':'),allow_nan=False)+'\n').encode()


def main():
    root=Path(__file__).resolve().parents[1];os.chdir(root)
    reservation=Path(sys.argv[1]);r=json.loads(reservation.read_text())
    spec_path=Path(r['specification']);assert digest(spec_path)==r['specification_sha256']
    spec=json.loads(spec_path.read_text());source=spec['source']
    for p,h in spec['inputs_sha256'].items():assert digest(p)==h
    assert digest(source['probe'])==source['probe_sha256']
    assert digest(source['library'])==source['library_sha256']
    work=Path(r['work_directory']);work.mkdir(exist_ok=False)
    data_dir=work/'isotherms';data_dir.mkdir()
    start=time.monotonic();own_cpu_start=time.process_time()
    deadline=datetime.fromisoformat(r['deadline_utc']).timestamp()
    if time.time()>=deadline:raise ValueError('allocation deadline already passed')
    stop=threading.Event();lock=threading.Lock();active={};child_cpu=0.
    planes=spec['planes'];ts,qs=spec['logT'],spec['logQ'];cov=spec['source_coverage']
    parent=json.loads(Path(spec['parent_manifest']).read_text())
    present={(v['hydrogen'],v['helium3']) for v in parent['planes']}
    planned=[(v['hydrogen'],v['helium3']) for v in planes]
    if len(set(planned))!=len(planned) or set(planned)&present:
        raise ValueError('duplicate composition or attempted regeneration of completed plane')
    target={(x,y) for x in spec['target_hydrogen'] for y in spec['target_helium3']}
    if present|set(planned)!=target:raise ValueError('new and existing planes do not form the declared family')
    pilot=json.loads(Path(spec['pilot_report']).read_text())
    if pilot['outcome']!='refinement_pilot_passed':raise ValueError('refinement pilot has not passed')
    mixtures={(x,y):mixture(x,y) for x,y in planned}
    def qvalues(x,t):
        return [q for q in qs if q<=cov['original_logQ_max'] or
                (t>=cov['minimum_added_logT'] and x<=cov['maximum_added_hydrogen'])]
    # Only terminated allocations may supply reusable results. Their receipts,
    # costs, source manifests and actual isotherms retain separate identities.
    # A missing observer handle is not evidence that a writer has terminated.
    reused={};prior_cpu=0.;prior_source_bytes=0;allocation_ids=set();ancestry={}
    for previous in r['prior_allocations']:
        receipt_path=Path(previous['receipt'])
        if digest(receipt_path)!=previous['receipt_sha256']:
            raise ValueError('prior allocation receipt changed')
        prior=json.loads(receipt_path.read_text())
        if prior['status']!='released' or prior.get('active_tasks'):
            raise ValueError('cannot reuse an allocation with active or uncertain writers')
        identity=receipt_path.resolve()
        if identity in allocation_ids:raise ValueError('duplicate prior allocation accounting')
        allocation_ids.add(identity);prior_cpu+=prior['cpu_seconds']
        for ancestor in prior.get('prior_allocations',[]):
            path=Path(ancestor['receipt']).resolve();value=ancestor['receipt_sha256']
            if path in ancestry and ancestry[path]!=value:
                raise ValueError('conflicting prior allocation ancestry')
            ancestry[path]=value
        manifest_path=Path(prior['source_manifest'])
        if digest(manifest_path)!=prior['source_manifest_sha256']:
            raise ValueError('prior source manifest changed')
        manifest=json.loads(manifest_path.read_text())
        if (manifest['specification_sha256']!=r['specification_sha256'] or
                manifest['source']!=source):
            raise ValueError('prior source identity or specification differs')
        seen=set()
        for saved in manifest['requests']:
            key=saved['task_id']
            if key in seen:raise ValueError('duplicate saved task within a prior manifest')
            seen.add(key)
            ip,it=divmod(key,len(ts))
            if not 0<=ip<len(planes):raise ValueError('saved source task outside specification')
            x,y=planes[ip]['hydrogen'],planes[ip]['helium3']
            if (saved['outcome']!='source_saved' or saved['plane_index']!=ip or
                    saved['temperature_index']!=it or saved['hydrogen']!=x or
                    saved['helium3']!=y or saved['logT']!=ts[it] or
                    saved['states']!=len(qvalues(x,ts[it]))):
                raise ValueError('saved source coordinates differ')
            path=Path(saved['source_file'])
            if digest(path)!=saved['source_sha256'] or path.stat().st_size!=saved['bytes']:
                raise ValueError('saved source isotherm changed')
            if key in reused and saved['source_sha256']!=reused[key]['source_sha256']:
                raise ValueError('conflicting saved source values')
            reused.setdefault(key,saved)
    if not allocation_ids:raise ValueError('continuation requires a terminated prior allocation')
    supplied={Path(v['receipt']).resolve():v['receipt_sha256'] for v in r['prior_allocations']}
    if any(supplied.get(path)!=value for path,value in ancestry.items()):
        raise ValueError('all ancestor allocations must be included in cumulative accounting')
    prior_source_bytes=sum(v['bytes'] for v in reused.values())
    r.update(prior_allocations_cpu_seconds=prior_cpu,reused_source_bytes=prior_source_bytes,
             source_states_reused=sum(v['states'] for v in reused.values()))
    queue=deque((ip*len(ts)+it,ip,it) for ip in range(len(planes)) for it in range(len(ts))
                if ip*len(ts)+it not in reused)
    if not queue:raise ValueError('source collection is already complete; do not launch duplicate work')
    results=list(reused.values());failures=[];futures={};consecutive=0
    r.update(status='running',controller_pid=os.getpid(),controller_sha256=digest(__file__),
             start_utc=now(),completed_requests=len(results),completed_states=sum(v['states'] for v in results),
             reused_requests=len(results),new_completed_requests=0,failed_requests=[],active_tasks=[])
    def save():
        for p in [reservation,work/'receipt.json']:
            tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(r,indent=2)+'\n');os.replace(tmp,p)
    def terminate(signum,frame):stop.set();r['stop_reason']='external termination signal'
    signal.signal(signal.SIGTERM,terminate);signal.signal(signal.SIGINT,terminate)
    save()
    def run(job):
        nonlocal child_cpu
        key,ip,it=job;x=planes[ip]['hydrogen'];y=planes[ip]['helium3'];t=ts[it];q=qvalues(x,t);mix=mixtures[x,y];scale=mix['source_mass_scale']
        request=' '.join(format(v,'.17g') for v in mix['eps'])+'\n3 223 -2\n'
        request+=''.join(f'{math.log(scale)+math.log(10)*(a+1.5*(t-6)):.17g} {math.log(10)*t:.17g}\n' for a in q)
        stem=data_dir/f'request-{key:05d}'
        infile=stem.with_suffix('.input.tmp');outfile=stem.with_suffix('.stdout.tmp');errfile=stem.with_suffix('.stderr.tmp')
        infile.write_text(request);begin=time.monotonic()
        with infile.open() as fi,outfile.open('x') as fo,errfile.open('x') as fe:
            proc=subprocess.Popen([source['probe']],stdin=fi,stdout=fo,stderr=fe,start_new_session=True,
                                  env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1'})
            with lock:active[key]={'pid':proc.pid,'command':[source['probe']],
                                   'input_file':str(infile),'start_utc':now(),'hydrogen':x,'helium3':y,'temperature_index':it}
            terminated=None
            while True:
                pid,status,usage=os.wait4(proc.pid,os.WNOHANG)
                if pid:
                    proc.returncode=os.waitstatus_to_exitcode(status);cpu=usage.ru_utime+usage.ru_stime
                    with lock:child_cpu+=cpu;active.pop(key,None)
                    break
                if (stop.is_set() or time.monotonic()-begin>r['maximum_isotherm_seconds']) and terminated is None:
                    os.killpg(proc.pid,signal.SIGTERM);terminated=time.monotonic()
                elif terminated is not None and time.monotonic()-terminated>15:os.killpg(proc.pid,signal.SIGKILL)
                time.sleep(.05)
        stdout=outfile.read_text();stderr=errfile.read_text()
        base={'task_id':key,'plane_index':ip,'hydrogen':x,'helium3':y,'logT':t,'temperature_index':it,
              'exit_code':proc.returncode,'cpu_seconds':cpu,'elapsed_seconds':time.monotonic()-begin,
              'end_utc':now(),'input_sha256':hashlib.sha256(request.encode()).hexdigest()}
        try:data=[list(map(float,l.split())) for l in stdout.splitlines()]
        except ValueError:data=[]
        if proc.returncode or len(data)!=len(q) or not all(len(row)==22 and all(map(math.isfinite,row)) for row in data):
            return {**base,'outcome':'source_process_or_output_failed','input_file':str(infile),
                    'stdout_file':str(outfile),'stderr_file':str(errfile)}
        flagged=[];first_law=[];nonpositive=[];maximum=0.
        for i,row in enumerate(data):
            row[2]/=scale
            for k in [5,6,9,10,11,12,13,16,17]:row[k]*=scale
            if row[0]!=0:flagged.append({'density_index':i,'info':row[0]});continue
            expected_rho=10**(q[i]+1.5*(t-6));expected_T=10**t
            if abs(row[2]/expected_rho-1)>1e-9 or abs(row[3]/expected_T-1)>1e-10:
                return {**base,'outcome':'source_coordinate_mismatch','density_index':i,
                        'input_file':str(infile),'stdout_file':str(outfile),'stderr_file':str(errfile)}
            rho,T,P=row[2:5];chit,Er,Et,Sr,St=row[8:13]
            if P<=0 or Et==0:nonpositive.append(i);continue
            defect=max(abs(rho*Er/P+chit-1),abs(T*St/Et-1),abs(rho*T*Sr/P+chit))
            maximum=max(maximum,defect)
            if defect>1e-7:first_law.append({'density_index':i,'maximum_absolute_defect':defect})
        record={'source_options':[3,223,-2],'actual_probe_sha256':source['probe_sha256'],
                'relative_electron_integral_target':source['relative_electron_integral_target'],
                'input':request,'input_sha256':base['input_sha256'],'returncode':proc.returncode,'stderr':stderr,
                'request':{'request_index':key,'plane_index':ip,'hydrogen':x,'helium3':y,'logT':t,'temperature_index':it,'logQ':q},
                'mixture':mix,'data':data}
        output=stem.with_suffix('.json.gz')
        with output.open('xb') as f:f.write(gzip.compress(encode(record),mtime=0))
        result={**base,'outcome':'source_saved','states':len(data),'source_file':str(output),
                'source_sha256':digest(output),'bytes':output.stat().st_size,'origin':'new_composition_refinement_source',
                'nonconvergence_flags':flagged,'first_law_failures':first_law,
                'nonpositive_pressure_or_zero_heat_response':nonpositive,'maximum_first_law_defect':maximum}
        # Converted source values and the exact request are preserved above.
        # Remove only this collector's transient duplicate process streams.
        infile.unlink();outfile.unlink();errfile.unlink()
        return result
    next_check=start
    with (work/'completed.jsonl').open('x') as journal:
        for result in results:journal.write(json.dumps(result,allow_nan=False)+'\n')
        journal.flush()
        with ThreadPoolExecutor(max_workers=r['threads']) as pool:
            while queue or futures:
                if time.time()>=deadline:stop.set();r['stop_reason']='fixed allocation deadline'
                while queue and len(futures)<r['threads'] and not stop.is_set():
                    job=queue.popleft();futures[pool.submit(run,job)]=job
                done,_=wait(futures,timeout=.5,return_when=FIRST_COMPLETED)
                for future in done:
                    job=futures.pop(future)
                    try:result=future.result()
                    except Exception as error:result={'task_id':job[0],'outcome':'collector_exception','error':repr(error)}
                    journal.write(json.dumps(result,allow_nan=False)+'\n');journal.flush()
                    if result['outcome']=='source_saved':results.append(result);consecutive=0
                    else:
                        failures.append(result);consecutive+=1
                        if consecutive>=3:stop.set();r['stop_reason']='three consecutive source-process failures'
                if time.monotonic()>=next_check:
                    with lock:live=[{'task_id':k,**v} for k,v in active.items()];used_cpu=child_cpu+time.process_time()-own_cpu_start
                    size=sum(v['bytes'] for v in reused.values())
                    for p in work.rglob('*'):
                        try:
                            if p.is_file():size+=p.stat().st_size
                        except FileNotFoundError:pass
                    r.update(observed_utc=now(),completed_requests=len(results),new_completed_requests=len(results)-len(reused),
                             completed_states=sum(v['states'] for v in results),failed_requests=[v['task_id'] for v in failures],
                             pending_requests=len(queue),active_tasks=live,cpu_seconds=used_cpu,
                             elapsed_seconds=time.monotonic()-start,retained_and_transient_bytes=size,
                             source_nonconvergence_states=sum(len(v['nonconvergence_flags']) for v in results),
                             source_first_law_failures=sum(len(v['first_law_failures']) for v in results))
                    if used_cpu>=r['maximum_cpu_seconds']:stop.set();r['stop_reason']='CPU allocation limit'
                    if size>=r['data_cap_bytes']:stop.set();r['stop_reason']='source storage limit'
                    save();print(json.dumps({k:r[k] for k in ['completed_requests','completed_states','pending_requests',
                        'cpu_seconds','elapsed_seconds','retained_and_transient_bytes','source_first_law_failures']}),flush=True)
                    next_check=time.monotonic()+30
                if stop.is_set() and not futures:break
    manifest={'scope':__doc__,'specification':str(spec_path),'specification_sha256':digest(spec_path),
              'accepted_for_evolution':False,'source':source,'requests':sorted(results,key=lambda v:v['task_id']),
              'failed_requests':failures,'unstarted_task_ids':[j[0] for j in queue],
              'all_requested_sources_returned':not failures and not queue,
              'source_nonconvergence_states':sum(len(v['nonconvergence_flags']) for v in results),
              'source_first_law_failures':sum(len(v['first_law_failures']) for v in results),
              'source_nonpositive_pressure_or_zero_heat_response':sum(len(v['nonpositive_pressure_or_zero_heat_response']) for v in results)}
    target=work/'source_manifest.json';target.write_text(json.dumps(manifest,indent=2)+'\n')
    r.update(status='released',end_utc=now(),active_tasks=[],completed_requests=len(results),
             completed_states=sum(v['states'] for v in results),pending_requests=len(queue),
             failed_requests=[v['task_id'] for v in failures],
             source_nonconvergence_states=manifest['source_nonconvergence_states'],
             source_first_law_failures=manifest['source_first_law_failures'],
             cpu_seconds=child_cpu+time.process_time()-own_cpu_start,elapsed_seconds=time.monotonic()-start,
             source_manifest=str(target),source_manifest_sha256=digest(target),
             outcome='source_collection_complete_validation_pending' if not failures and not queue else 'incomplete_sources_preserved')
    save()
    if failures or queue:sys.exit(1)


if __name__=='__main__':main()
