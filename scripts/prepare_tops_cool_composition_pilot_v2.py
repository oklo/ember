#!/usr/bin/env python3
"""Repair the pilot's source template, retaining all unsubmitted attempts.

The v1 template contained one temperature and failed local coordinate checks
before any opacity submission. Use the existing full native source template
at the same composition and explicitly verify its temperature list. Requested
physics, density coordinates and photon spacing are unchanged.
"""
import json
from pathlib import Path
import shutil
from datetime import datetime, timezone

from fetch_tops_composition import digest
from import_tops_composition import read
from reduce_tops_group_factors import add_inputs, verify


def write(path, value):
    with Path(path).open('x') as stream:
        stream.write(json.dumps(value, indent=2, allow_nan=False)+'\n')


def main():
    root=Path(__file__).resolve().parents[1]
    reservations=root/'docs/research/fable/coordination/reservations'
    for path in reservations.glob('*.json'):
        if json.loads(path.read_text()).get('status') not in ('released','complete','completed','cancelled','failed','superseded'):
            raise ValueError('inspect active reservation: '+str(path))
    prior_path=reservations/'primary_tops_cool_composition_pilot_v1.json'
    prior=json.loads(prior_path.read_text())
    if prior['completed_tasks'] != 0 or prior['status'] != 'released':
        raise ValueError('inspect completed or active original requests')
    old_comparison=root/'data/opacity/sources/tops_cool_composition_pilot_comparison_v1_specification.json'
    old=json.loads(old_comparison.read_text())
    inputs={str(p.resolve()):digest(p) for p in (prior_path,old_comparison,Path(__file__))}
    add_inputs(inputs,prior['input_sha256'])
    mf=root/'data/opacity/sources/hydrogen_exhaustion/family_manifest.json'
    manifest=json.loads(mf.read_text());add_inputs(inputs,{str(mf.resolve()):digest(mf)})
    work=Path('/tmp/ember-tops-cool-composition-pilot-v2')
    monitor=Path('/tmp/ember-tops-cool-composition-pilot-run-v2');monitor.mkdir()
    seeds={};tasks=[];plans=[];reductions=[];native_checks=[]
    for item in old['reductions']:
        key=item['X'],item['Z']
        if key not in seeds:
            entries=[r for r in manifest['planes'] if (r['X'],r['Z'])==key]
            if len(entries)!=1:raise ValueError('missing full source template')
            entry=entries[0]
            src,req=mf.parent/entry['file'],mf.parent/entry['request']
            tt,_,_,_=read(src,entry,dimensions=(50,71))
            if digest(req)!=entry['request_sha256']:raise ValueError('source request changed')
            directory=work/f'seed-x{round(key[0]*1000):03d}-z{round(key[1]*1000):03d}'
            directory.mkdir(parents=True)
            shutil.copyfile(src,directory/'source.txt');shutil.copyfile(req,directory/'request.json')
            receipt={**entry,'file':'source.txt','request':'request.json','dimensions':[50,71]}
            write(directory/'receipt.json',receipt)
            seeds[key]=directory,tt
            for p in (src,req,directory/'source.txt',directory/'request.json',directory/'receipt.json'):
                add_inputs(inputs,{str(p.resolve()):digest(p)})
        seed,tt=seeds[key]
        old_path=root/item['plan'];plan=json.loads(old_path.read_text())
        old_work=Path(plan['requests'][0]['work'])
        if any((old_work/n).exists() for n in ('source.txt','results.html','receipt.json','recipe.json')):
            raise ValueError('original request may have been submitted; inspect before retry')
        if any(t not in tt for job in plan['requests'] for t in job['temperatures_keV']):
            raise ValueError('requested source temperature is not native')
        path=Path(str(old_path).replace('_v1_specification.json','_v2_specification.json'))
        job={**plan['requests'][0],'work':str(work/old_work.name)}
        plan={**plan,'scope':plan['scope']+'\n'+__doc__,'baseline':str(seed),'requests':[job]}
        write(path,plan);plans.append(str(path.relative_to(root)))
        add_inputs(inputs,{str(path.resolve()):digest(path)})
        tasks.append({'kind':'groups','plan':str(path.relative_to(root)),'request_index':0,'work':job['work']})
        reductions.append({**item,'plan':str(path.relative_to(root)),
                           'report':item['report'].replace('_reduction_v1.json','_reduction_v2.json')})
        native_checks.append({'plan':str(path.relative_to(root)),'native_temperature_count':len(tt),
                              'requested_temperatures_native':True,'original_unsubmitted_work':str(old_work)})
    comparison=root/'data/opacity/sources/tops_cool_composition_pilot_comparison_v2_specification.json'
    write(comparison,{**old,'scope':old['scope']+'\n'+__doc__,'reductions':reductions})
    add_inputs(inputs,{str(comparison.resolve()):digest(comparison)})
    verify(inputs)
    tf,ef=monitor/'tasks.json',monitor/'existing_sources.json'
    write(tf,tasks);write(ef,{'prior_request_inventory':prior['existing_sources_file'],'native_template_checks':native_checks})
    reservation=reservations/'primary_tops_cool_composition_pilot_v2.json'
    write(reservation,{'owner':'primary','task_ids':['E-OPACITY-COOL-COMPOSITION-PILOT'],
        'utc':datetime.now(timezone.utc).isoformat(),'status':'reserved','threads':4,
        'expected_memory_bytes':1200000000,'data_cap_bytes':100000000,'maximum_elapsed_seconds':900,
        'maximum_attempts_per_task':3,'stop_after_consecutive_failed_tasks':3,
        'scratch_paths':[str(work),str(monitor)],'work_directory':str(monitor),
        'task_file':str(tf),'task_file_sha256':digest(tf),
        'existing_sources_file':str(ef),'existing_sources_sha256':digest(ef),
        'input_sha256':inputs,'plans':plans})
    controller=Path('/tmp/ember-tops-cool-frequency-scaled-run-v1.py').read_text().replace(
        'primary_tops_cool_frequency_scaled_v1.json',reservation.name)
    with Path('/tmp/ember-tops-cool-composition-pilot-run-v2.py').open('x') as stream:stream.write(controller)
    print(json.dumps({'requests':len(tasks),'source_states':2*len(tasks),'native_temperature_checks':len(native_checks)}))


if __name__=='__main__':main()
