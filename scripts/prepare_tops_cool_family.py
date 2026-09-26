#!/usr/bin/env python3
"""Prepare complete cold/transition source grids with explicit source reuse.

Cover the original cold atomic density range through 251.2 g/cm3, retaining
independent controls. Use separate fine grey and plasma-ratio composition
grids. Group-accompanying grey means replace duplicate grey requests. This
source specification does not accept its eventual interpolation tables.
"""
import json
import math
from pathlib import Path
import shutil
from datetime import datetime, timezone

from fetch_tops_composition import digest
from import_tops_composition import read
from reduce_tops_group_factors import add_inputs, verify
from tops_hydrogen_request import hydrogen_request


def write(path, value):
    with Path(path).open('x') as stream:
        stream.write(json.dumps(value, indent=2, allow_nan=False)+'\n')


def main():
    root=Path(__file__).resolve().parents[1]
    reservations=root/'docs/research/fable/coordination/reservations'
    for p in reservations.glob('*.json'):
        if json.loads(p.read_text()).get('status') not in ('released','complete','completed','cancelled','failed','superseded'):
            raise ValueError('inspect active reservation: '+str(p))
    inputs={}
    for name in ('tops_cool_composition_pilot_v3.json','tops_cool_full_spectra_v4.json',
                 'tops_cool_interpolation_components_v2.json'):
        p=root/'docs/results'/name;r=json.loads(p.read_text())
        add_inputs(inputs,r['input_sha256']);add_inputs(inputs,{str(p.resolve()):digest(p)})
    mf=root/'data/opacity/sources/hydrogen_exhaustion/family_manifest.json'
    manifest=json.loads(mf.read_text());add_inputs(inputs,{str(mf.resolve()):digest(mf)})
    mixture={(r['X'],r['Z']):r for r in manifest['planes']}
    xs=[0.,.025,.05,.075,.125,.15,.175,.2,.25,.275,.325,.4,.5,.6,.65,.75]
    zs=[.01,.03]
    tt=[.002,.0025,.003,.0035,.004,.005,.006,.007,.008,.009,.01,.0125,.015,.02,.025,.03,.04,.05,.06]
    cold=tt[:14]
    # Two dilute decades, six 0.4-dex coordinates, then the contiguous 0.05-dex
    # grid to the original cold density limit. Every request is logarithmic.
    dense=[10**(-1.95+.05*i) for i in range(88)]
    chunks=[[1e-10,1e-6],[10**(-4+.4*i) for i in range(6)],dense[:44],dense[44:]]
    rr=[r for chunk in chunks for r in chunk]
    if len(rr)!=96 or rr!=sorted(set(rr)) or any(len(chunk)>72 for chunk in chunks):
        raise ValueError('invalid source density partition')
    work=Path('/tmp/ember-tops-cool-family-v1')
    monitor=Path('/tmp/ember-tops-cool-family-run-v1');monitor.mkdir()
    seeds={};tasks=[];planes=[];plans=[]
    for x in xs:
        for z in zs:
            key=(x,z) if (x,z) in mixture else (.25,z)
            if key not in seeds:
                entry=mixture[key];src=mf.parent/entry['file'];req=mf.parent/entry['request']
                native,_,_,_=read(src,entry,dimensions=(50,71))
                if any(t not in native for t in tt) or digest(req)!=entry['request_sha256']:
                    raise ValueError('invalid full native source template')
                seed=work/f'seed-x{round(key[0]*1000):03d}-z{round(z*1000):03d}'
                seed.mkdir(parents=True)
                shutil.copyfile(src,seed/'source.txt');shutil.copyfile(req,seed/'request.json')
                write(seed/'receipt.json',{**entry,'file':'source.txt','request':'request.json','dimensions':[50,71]})
                for p in (src,req,seed/'source.txt',seed/'request.json',seed/'receipt.json'):
                    add_inputs(inputs,{str(p.resolve()):digest(p)})
                seeds[key]=seed
            seed=seeds[key];entry=json.loads((seed/'receipt.json').read_text())
            label=f'x{round(x*1000):03d}-z{round(z*1000):03d}'
            temperature_plans=[]
            for t in tt:
                requests=[{'temperatures_keV':[t],'densities_g_cm3':chunk,
                           'density_spacing':'log','work':str(work/label/f't{round(t*1e6):06d}'/f'density-{i}')}
                          for i,chunk in enumerate(chunks)]
                plan={'scope':__doc__,'baseline':str(seed),'hydrogen_atomic_mass_fraction':x,
                      'photon_min_keV':1e-8,'photon_max_keV':800*t,'photon_boundaries':993,
                      'maximum_response_bytes':20000000,'maximum_group_states':100000,'requests':requests}
                _,expected=hydrogen_request(plan,entry,json.loads((seed/'request.json').read_text()))
                if (expected['X'],expected['Z'])!=(x,z):raise ValueError('source target differs')
                p=root/'data/opacity/sources'/f'tops_cool_family_{label}_t{round(t*1e6):06d}_v1_specification.json'
                write(p,plan);add_inputs(inputs,{str(p.resolve()):digest(p)})
                temperature_plans.append({'plan':str(p.relative_to(root)),'sha256':digest(p)})
                plans.append(str(p.relative_to(root)))
                tasks.extend({'kind':'groups','plan':str(p.relative_to(root)),'request_index':i,'work':job['work']}
                             for i,job in enumerate(requests))
            p=root/'data/opacity/sources'/f'tops_cool_family_{label}_v1_specification.json'
            write(p,{'scope':__doc__,'X':x,'Z':z,'temperatures_keV':tt,
                     'densities_atomic_g_cm3_requested':rr,'temperature_plans':temperature_plans})
            add_inputs(inputs,{str(p.resolve()):digest(p)})
            planes.append({'X':x,'Z':z,'plan':str(p.relative_to(root)),'sha256':digest(p),
                           'reduction_report':f'docs/results/tops_cool_family_{label}_reduction_v1.json',
                           'output_directory':f'/tmp/ember-cool-family-reduced-v1/{label}'})
    # Existing group calculations have a different frequency grid or density.
    # Explicitly check matching-frequency plans to prevent overlapping requests.
    prior_states=[]
    for p in (root/'data/opacity/sources').glob('tops*_specification.json'):
        if p.name.startswith('tops_cool_family_'):continue
        spec=json.loads(p.read_text())
        if spec.get('photon_boundaries')!=993 or 'baseline' not in spec:continue
        b=json.loads((Path(spec['baseline'])/'receipt.json').read_text())
        x=spec.get('hydrogen_atomic_mass_fraction',b['X']);z=b['Z']
        if x not in xs or z not in zs:continue
        for job in spec['requests']:
            for t in job['temperatures_keV']:
                if t not in tt or (spec['photon_min_keV'],spec['photon_max_keV'])!=(1e-8,800*t):continue
                for rho in job['densities_g_cm3']:
                    if any(abs(rho/r-1)<5e-5 for r in rr):raise ValueError('reusable group coordinate: '+str(job))
                    prior_states.append({'X':x,'Z':z,'T_keV':t,'rho':rho,'work':job['work']})
    retained_grey=[]
    for label in ('010','030','070'):
        p=root/'docs/results'/f'tops_cool_grey_refined_x{label}-z020_v1.json'
        r=json.loads(p.read_text());add_inputs(inputs,r['input_sha256']);add_inputs(inputs,r['output_sha256'])
        add_inputs(inputs,{str(p.resolve()):digest(p)})
        retained_grey.append({'X':r['X'],'Z':r['Z'],'report':str(p.relative_to(root)),'sha256':digest(p)})
    supplied_by_groups=set(mixture)&{(x,z) for x in xs for z in zs}
    retained_keys={(r['X'],r['Z']) for r in retained_grey}
    missing=sorted(set(mixture)-supplied_by_groups-retained_keys)
    grey_requests=[]
    for x,z in missing:
        for offset in (1,2,3):
            densities=[10**(-2+.05*offset+.2*i) for i in range(22)]
            for it in range(0,len(cold),2):
                grey_requests.append({'X':x,'Z':z,'temperatures_keV':cold[it:it+2],
                    'densities_g_cm3':densities,'work':str(work/f'grey-x{round(x*1000):03d}-z{round(z*1000):03d}'/f'offset-{offset}-t{it:02d}')})
    grey_plan=root/'data/opacity/sources/tops_cool_family_grey_refinement_v1_specification.json'
    write(grey_plan,{'scope':__doc__,'baseline_manifest':str(mf),'baseline_manifest_sha256':digest(mf),
        'maximum_response_bytes':20000000,'maximum_states_per_request':72,'requests':grey_requests,
        'group_supplied_compositions':[list(k) for k in sorted(supplied_by_groups)],'retained_grey':retained_grey})
    add_inputs(inputs,{str(grey_plan.resolve()):digest(grey_plan)})
    tasks.extend({'kind':'grey','plan':str(grey_plan.relative_to(root)),'request_index':i,'work':j['work']}
                 for i,j in enumerate(grey_requests));plans.append(str(grey_plan.relative_to(root)))
    family=root/'data/opacity/sources/tops_cool_family_v1_specification.json'
    write(family,{'scope':__doc__,'ratio_X':xs,'ratio_Z':zs,'temperatures_keV':tt,
        'cold_temperatures_keV':cold,'density_spacing_dex':.05,'maximum_density_atomic_g_cm3_requested':rr[-1],
        'ratio_planes':planes,'grey_baseline_plan':'data/opacity/sources/tops_native_uncut_family_batched_form_v4_specification.json',
        'hot_grey_refinement_plan':'data/opacity/sources/tops_grey_density_refinement_v1_specification.json',
        'grey_refinement_plan':str(grey_plan.relative_to(root)),'retained_grey':retained_grey,
        'independent_controls':'docs/results/tops_cool_ratio_hydrogen_comparison_v1.json',
        'accepted_for_stellar_opacity':False})
    add_inputs(inputs,{str(family.resolve()):digest(family)})
    for p in (Path(__file__),*(root/'scripts'/name for name in
              ('tops_hydrogen_request.py','prepare_tops_form_lease_v4.py','fetch_tops_form_plan_v4.py'))):
        add_inputs(inputs,{str(p.resolve()):digest(p)})
    verify(inputs)
    tf,ef=monitor/'tasks.json',monitor/'existing_sources.json';write(tf,tasks);write(ef,prior_states)
    reservation=reservations/'primary_tops_cool_family_v1.json'
    write(reservation,{'owner':'primary','task_ids':['E-OPACITY-COOL-FAMILY'],'utc':datetime.now(timezone.utc).isoformat(),
        'status':'reserved','threads':4,'expected_memory_bytes':1600000000,'data_cap_bytes':6000000000,
        'maximum_elapsed_seconds':7200,'maximum_attempts_per_task':3,'stop_after_consecutive_failed_tasks':3,
        'scratch_paths':[str(work),str(monitor)],'work_directory':str(monitor),
        'task_file':str(tf),'task_file_sha256':digest(tf),'existing_sources_file':str(ef),
        'existing_sources_sha256':digest(ef),'input_sha256':inputs,'plans':plans})
    controller=Path('/tmp/ember-tops-cool-frequency-scaled-run-v1.py').read_text().replace(
        'primary_tops_cool_frequency_scaled_v1.json',reservation.name).replace(
        'prepare_tops_form_lease_v3.py','prepare_tops_form_lease_v4.py').replace('fetch_tops_form_plan_v3.py','fetch_tops_form_plan_v4.py')
    with Path('/tmp/ember-tops-cool-family-run-v1.py').open('x') as stream:stream.write(controller)
    summary={'group_planes':len(planes),'group_states':len(xs)*len(zs)*len(tt)*len(rr),
        'group_requests':4*len(xs)*len(zs)*len(tt),'grey_refinement_compositions':len(missing),
        'grey_states':sum(len(j['temperatures_keV'])*len(j['densities_g_cm3']) for j in grey_requests),
        'grey_requests':len(grey_requests),'total_requests':len(tasks),
        'prior_matching_frequency_states_inspected':len(prior_states),
        'grey_planes_reusing_groups':len(supplied_by_groups),'grey_planes_reusing_completed_refinement':len(retained_grey)}
    write(root/'docs/results/tops_cool_family_preflight_v1.json',{'scope':__doc__,'prepared':True,
        'accepted_for_stellar_opacity':False,'summary':summary,'input_sha256':inputs})
    print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
