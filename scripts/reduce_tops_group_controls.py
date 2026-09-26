#!/usr/bin/env python3
"""Integrate independent group controls without constructing an opacity table.

Control coordinates need not form an interpolation grid. Preserve every source
and quadrature result, including failed native-grey recovery diagnostics.
"""
import argparse
import json
from pathlib import Path
from fetch_tops_composition import digest
from tops_groups import source_groups
from reduce_tops_group_factors import adaptive_state, validate_coordinates, verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('plan',type=Path);p.add_argument('report',type=Path)
    a=p.parse_args()
    if a.report.exists():raise FileExistsError('preserve completed control reductions')
    plan=json.loads(a.plan.read_text())
    if plan.get('retained_requests'):raise ValueError('this control reducer expects explicit source requests')
    inputs={str(a.plan.resolve()):digest(a.plan)}
    for name in ('reduce_tops_group_controls.py','tops_groups.py','reduce_tops_group_factors.py',
                 'reduce_tops_group_plan.py','audit_tops_electron_dispersion.py',
                 'audit_tops_spectral_means.py','import_tops_composition.py','fetch_tops_group_plan.py'):
        path=Path(__file__).with_name(name);inputs[str(path.resolve())]=digest(path)
    baseline=json.loads((Path(plan['baseline'])/'receipt.json').read_text())
    inputs[str((Path(plan['baseline'])/'receipt.json').resolve())]=digest(Path(plan['baseline'])/'receipt.json')
    rows=[];seen=set()
    for job in plan['requests']:
        root=Path(job['work']);r=json.loads((root/'receipt.json').read_text())
        if any(r[k]!=baseline[k] for k in ('X','Z','metals')):raise ValueError('group mixture differs')
        for n in ('receipt.json','request.json','source.txt','recipe.json'):
            path=root/n;inputs[str(path.resolve())]=digest(path)
        boundaries,groups,grey=source_groups(root)
        if len(boundaries)!=plan['photon_boundaries']:raise ValueError('frequency boundary count differs')
        part=[]
        for (t,rho),data in groups.items():
            if (t,rho) in seen:raise ValueError('duplicate control state')
            seen.add((t,rho));row=adaptive_state(boundaries,data,t,rho,grey[t,rho],root)
            row['refractive_ratio']=row['rosseland_atomic_cm2_g']/row['uncut_recombined_rosseland']
            row['native_uncut_times_ratio_atomic_cm2_g']=row['uncut_source_rosseland']*row['refractive_ratio']
            part.append(row)
        validate_coordinates(part,job);rows.extend(part)
    verify(inputs)
    result={'scope':__doc__,'accepted_for_stellar_opacity':False,'X':baseline['X'],'Z':baseline['Z'],
            'states':len(rows),'photon_boundaries':plan['photon_boundaries'],
            'maximum_uncut_recovery_relative_error':max(abs(r['uncut_recovery_relative_error']) for r in rows),
            'uncut_recovery_check_passed':all(abs(r['uncut_recovery_relative_error'])<=1e-4 for r in rows),
            'records':rows,'input_sha256':inputs}
    a.report.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('scope','records','input_sha256')}))


if __name__=='__main__':main()
