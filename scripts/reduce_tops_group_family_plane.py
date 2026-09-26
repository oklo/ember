#!/usr/bin/env python3
"""Integrate a complete ratio plane from temperature-scaled source requests.

Each temperature has its own photon range. Source mixtures, rectangular
coordinates, frequency coverage and adaptive quadrature remain explicit.
Outputs are candidates for interpolation checks, not accepted stellar inputs.
"""
import argparse
from collections import Counter
import json
import math
from pathlib import Path

from fetch_tops_composition import digest
from reduce_tops_group_factors import (
    add_inputs, verify, adaptive_state, validate_coordinates, table_text,
    QUADRATURE_CRITERION, UNCUT_RECOVERY_CRITERION)
from tops_groups import source_groups
from tops_hydrogen_request import hydrogen_request


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('plan',type=Path)
    p.add_argument('output_directory',type=Path)
    p.add_argument('report',type=Path)
    a=p.parse_args()
    if a.output_directory.exists() or a.report.exists():
        raise FileExistsError('preserve completed plane reductions')
    spec=json.loads(a.plan.read_text())
    inputs={str(a.plan.resolve()):digest(a.plan)}
    for name in ('reduce_tops_group_family_plane.py','reduce_tops_group_factors.py',
                 'tops_groups.py','tops_hydrogen_request.py','reduce_tops_group_plan.py',
                 'audit_tops_electron_dispersion.py','audit_tops_spectral_means.py',
                 'import_tops_composition.py','fetch_tops_composition.py'):
        path=Path(__file__).with_name(name);add_inputs(inputs,{str(path.resolve()):digest(path)})
    jobs=[]
    for item in spec['temperature_plans']:
        path=Path(item['plan'])
        if digest(path)!=item['sha256']:raise ValueError('temperature plan changed')
        add_inputs(inputs,{str(path.resolve()):digest(path)})
        plan=json.loads(path.read_text());seed=Path(plan['baseline'])
        baseline=json.loads((seed/'receipt.json').read_text())
        request=json.loads((seed/baseline['request']).read_text())
        if digest(seed/baseline['request'])!=baseline['request_sha256']:
            raise ValueError('baseline request changed')
        _,mixture=hydrogen_request(plan,baseline,request)
        if (mixture['X'],mixture['Z'])!=(spec['X'],spec['Z']):
            raise ValueError('source plane mixture differs')
        for q in (seed/'receipt.json',seed/baseline['request']):
            add_inputs(inputs,{str(q.resolve()):digest(q)})
        for job in plan['requests']:
            if len(job['temperatures_keV'])!=1:raise ValueError('expected one native temperature per source')
            t=job['temperatures_keV'][0]
            if (plan['photon_min_keV'],plan['photon_max_keV'],plan['photon_boundaries'])!=(1e-8,800*t,993):
                raise ValueError('source frequency grid differs from verified prescription')
            jobs.append((job,mixture,plan))
    if len({j['work'] for j,_,_ in jobs})!=len(jobs):raise ValueError('duplicate source ownership')
    missing=[j['work'] for j,_,_ in jobs if not (Path(j['work'])/'receipt.json').is_file()]
    if missing:raise ValueError('required source incomplete: '+str(missing[:3]))
    verify(inputs)
    records={}
    for job,mixture,plan in jobs:
        directory=Path(job['work']);receipt=json.loads((directory/'receipt.json').read_text())
        if any(receipt[k]!=v for k,v in mixture.items()):raise ValueError('returned mixture differs')
        for name in ('receipt.json','request.json','source.txt','recipe.json'):
            path=directory/name;add_inputs(inputs,{str(path.resolve()):digest(path)})
        recipe=json.loads((directory/'recipe.json').read_text())
        if recipe['job']!=job:raise ValueError('returned recipe differs from source plan')
        add_inputs(inputs,recipe.get('input_sha256',{}))
        boundaries,groups,grey=source_groups(directory)
        if len(boundaries)!=993 or abs(boundaries[0]/plan['photon_min_keV']-1)>1e-14 or abs(boundaries[-1]/plan['photon_max_keV']-1)>1e-14:
            raise ValueError('actual frequency bounds differ')
        part=[]
        for (t,rho),data in groups.items():
            row=adaptive_state(boundaries,data,t,rho,grey[t,rho],directory)
            row['refractive_ratio']=row['rosseland_atomic_cm2_g']/row['uncut_recombined_rosseland']
            row['native_uncut_times_ratio_atomic_cm2_g']=row['uncut_source_rosseland']*row['refractive_ratio']
            if not math.isfinite(row['refractive_ratio']) or row['refractive_ratio']<1-1e-12:
                raise ValueError('invalid plasma correction')
            key=t,rho
            if key in records:raise ValueError('duplicate source state')
            records[key]=row;part.append(row)
        validate_coordinates(part,job)
        print(json.dumps({'states':len(records),'source':str(directory)}),flush=True)
    tt=sorted({t for t,r in records});rr=sorted({r for t,r in records})
    intended=spec['densities_atomic_g_cm3_requested']
    if tt!=spec['temperatures_keV'] or len(rr)!=len(intended) or any(abs(r/q-1)>5e-5 for r,q in zip(rr,intended,strict=True)):
        raise ValueError('complete source axes differ from plane plan')
    if set(records)!={(t,r) for t in tt for r in rr}:
        raise ValueError('complete source plane has missing or substituted states')
    texts={}
    for field,name,description in (
        ('rosseland_atomic_cm2_g','group_transport.dat','ATOMIC refractive group transport'),
        ('refractive_ratio','factor.dat','DIMENSIONLESS refractive ratio; not standalone opacity'),
        ('native_uncut_times_ratio_atomic_cm2_g','normalized_transport.dat','ATOMIC native uncut grey times plasma ratio')):
        texts[name],prefixes=table_text(records,tt,rr,spec['X'],spec['Z'],field,description)
    verify(inputs)
    a.output_directory.mkdir(parents=True)
    for name,text in texts.items():(a.output_directory/name).write_text(text)
    maximum=max(abs(r['uncut_recovery_relative_error']) for r in records.values())
    result={'scope':__doc__,'accepted_for_stellar_opacity':False,'X':spec['X'],'Z':spec['Z'],
        'states':len(records),'temperatures_keV':tt,'densities_atomic_g_cm3':rr,'density_prefix_sizes':prefixes,
        'quadrature_relative_criterion':QUADRATURE_CRITERION,
        'maximum_quadrature_relative_change':max(r['quadrature_relative_change'] for r in records.values()),
        'new_quadrature_orders':dict(Counter(r['quadrature_order'] for r in records.values())),
        'uncut_recovery_relative_criterion':UNCUT_RECOVERY_CRITERION,
        'maximum_uncut_recovery_relative_error':maximum,'uncut_recovery_check_passed':maximum<=UNCUT_RECOVERY_CRITERION,
        'normalization_is_an_independent_accuracy_check':False,'records':list(records.values()),
        'input_sha256':inputs,'output_sha256':{str((a.output_directory/n).resolve()):digest(a.output_directory/n) for n in texts},
        'limitations':['Full-spectrum and final interpolation comparisons are separate requirements.',
                       'Quadrature convergence does not establish opacity physics or group-frequency accuracy.',
                       'Uncut-mean recovery errors are retained; normalization is not an independent test.',
                       'No conductive opacity or stellar trajectory is changed by this reduction.']}
    a.report.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('X','Z','states','maximum_quadrature_relative_change',
        'maximum_uncut_recovery_relative_error','uncut_recovery_check_passed')}),flush=True)


if __name__=='__main__':main()
