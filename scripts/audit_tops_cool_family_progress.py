#!/usr/bin/env python3
"""Account for validated, pending and unstarted cold-family source requests.

Use terminal execution receipts and recheck every claimed completed source
hash. Failed observations remain pending; this does not infer whether an
unobservable remote calculation is still running. No request is resubmitted.
"""
import argparse
from collections import Counter
import json
from pathlib import Path

from fetch_tops_composition import digest
from reduce_tops_group_factors import add_inputs, verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    p.add_argument('--controller',action='append',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('preserve completed source inventories')
    inputs={str(Path(__file__).resolve()):digest(Path(__file__))}
    events={};original=None;completed={};controllers=[]
    for path in a.controller:
        r=json.loads(path.read_text())
        if r['status']!='released':raise ValueError('inventory expects terminal controllers')
        tf=Path(r['task_file']);result_path=Path(r['work_directory'])/'task_results.jsonl'
        if digest(tf)!=r['task_file_sha256'] or digest(result_path)!=r['task_results_sha256']:
            raise ValueError('execution archive changed')
        if original is None:original=json.loads(tf.read_text())
        for file in (path,tf,result_path):add_inputs(inputs,{str(file.resolve()):digest(file)})
        add_inputs(inputs,r['input_sha256'])
        for line in result_path.read_text().splitlines():
            e=json.loads(line);i=e.get('original_task_id',e['task_id'])
            if not 0<=i<len(original):raise ValueError('unknown original request')
            old=original[i]
            pp,op=Path(e['plan']),Path(old['plan'])
            plan,oldplan=json.loads(pp.read_text()),json.loads(op.read_text())
            job,oldjob=plan['requests'][e['request_index']],oldplan['requests'][old['request_index']]
            if e['kind']!=old['kind'] or {k:v for k,v in job.items() if k!='work'}!={k:v for k,v in oldjob.items() if k!='work'}:
                raise ValueError('recovery changed physical source coordinates')
            for field in ('baseline','baseline_manifest','baseline_manifest_sha256',
                          'hydrogen_atomic_mass_fraction','photon_boundaries','photon_min_keV','photon_max_keV'):
                if plan.get(field)!=oldplan.get(field):raise ValueError('recovery changed source mixture or frequencies')
            for file in (pp,op):add_inputs(inputs,{str(file.resolve()):digest(file)})
            events.setdefault(i,[]).append(e)
            if e['outcome']=='completed':
                if i in completed:raise ValueError('original physical source reported completed twice')
                completed[i]=e
        controllers.append({'receipt':str(path),'outcome':r['outcome'],
                            'completed_tasks':r['completed_tasks'],'failed_tasks':len(r['failed_tasks'])})
    rows=[];states=0
    for i,task in enumerate(original):
        history=events.get(i,[])
        if i in completed:
            e=completed[i];directory=Path(e['work']);rp=directory/'receipt.json'
            receipt=json.loads(rp.read_text())
            if digest(rp)!=e['source_receipt_sha256'] or receipt['sha256']!=e['source_sha256']:
                raise ValueError('completed source receipt changed')
            for file,expected in ((directory/receipt['file'],receipt['sha256']),
                                  (directory/receipt['request'],receipt['request_sha256'])):
                if digest(file)!=expected:raise ValueError('completed numeric source changed')
                add_inputs(inputs,{str(file.resolve()):expected})
            add_inputs(inputs,{str(rp.resolve()):digest(rp)})
            count=e['validated_source']['states']
            if count!=receipt['dimensions'][0]*receipt['dimensions'][1]:raise ValueError('source state count differs')
            states+=count
            row={'status':'validated','work':str(directory),'states':count}
        elif history:
            e=history[-1];directory=Path(e['work'])
            attempts=e.get('attempts',[]);last=attempts[-1][-1] if attempts and attempts[-1] else {}
            if last.get('phase')=='allocate':status='allocation_failed_fetch_not_started'
            elif (directory/'submit.html').exists():status='pending_existing_result'
            else:status='submission_or_execution_unconfirmed'
            row={'status':status,'work':str(directory),'last_outcome':e['outcome'],
                 'last_observation':e.get('stop_reason',e.get('error')),
                 'last_phase':last.get('phase')}
        else:row={'status':'not_started','work':task['work']}
        rows.append({'original_task_id':i,'kind':task['kind'],**row})
    verify(inputs)
    result={'scope':__doc__,'accepted_for_stellar_opacity':False,'all_source_work_complete':len(completed)==len(original),
        'source_requests':len(original),'status_counts':dict(Counter(r['status'] for r in rows)),
        'validated_source_states':states,'controllers':controllers,'requests':rows,'input_sha256':inputs}
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('source_requests','status_counts','validated_source_states')}),flush=True)


if __name__=='__main__':main()
