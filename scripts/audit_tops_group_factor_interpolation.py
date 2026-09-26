#!/usr/bin/env python3
"""Test composition interpolation of the refractive opacity ratio separately.

F is the group opacity with refraction divided by the same groups' uncut
opacity. A candidate can multiply an independently tabulated UNCUT grey mean
by F, retaining the source's averaging convention when F tends to one. This
must never multiply the current conditional, plasma-cutoff grey tables.
"""
import argparse
import json
import math
from pathlib import Path

from audit_tops_spectral_means import digest
from reduce_tops_group_plan import uncut_mean
from tops_groups import source_groups


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('compositions', 'base', 'hot', 'output'):
        p.add_argument(name,type=Path)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('keep completed factor checks')
    reports = [json.loads(path.read_text()) for path in (a.compositions,a.base,a.hot)]
    inputs = {}
    for report in reports:
        for name,h in report['input_sha256'].items():
            if name in inputs and inputs[name] != h:
                raise ValueError('inconsistent reference inputs')
            inputs[name] = h
    for path in (a.compositions,a.base,a.hot,Path(__file__)):
        inputs[str(path.resolve())] = digest(path)
    for name,h in inputs.items():
        if digest(Path(name)) != h:
            raise ValueError('changed source input: '+name)
    records,cache = {},{}

    def grouped(root,key):
        if root not in cache:
            cache[root] = source_groups(root)
        b,data,grey=cache[root]
        return uncut_mean(b,data[key],key[0]),grey[key]['rosseland']

    normalized_checks=[]
    for row in reports[0]['records']:
        key=row['temperature_keV'],row['density_atomic_g_cm3']
        uncut,native=grouped(row['group_source'],key)
        factor=row['group_rosseland_atomic_cm2_g']/uncut
        records[row['X'],row['Z'],*key]={'factor':factor,'source':row['group_source']}
        normalized_checks.append({'X':row['X'],'Z':row['Z'],
            'temperature_keV':key[0],'density_atomic_g_cm3':key[1],
            'group_factor':factor,'native_uncut_mean':native,
            'source_normalized_mean':native*factor,
            'relative_to_log_spectrum_reference':native*factor/row['reference_rosseland_atomic_cm2_g']-1})
    wanted={(r['temperature_keV'],r['density_atomic_g_cm3']) for r in reports[0]['records']}
    for row in reports[1]['records']:
        key=row['temperature_keV'],row['density_atomic_g_cm3']
        if key in wanted:
            records[0,.02,*key]={'factor':row['rosseland_atomic_cm2_g']/row['uncut_recombined_rosseland'],
                                'source':row['source']}
    for row in reports[2]['records']:
        key=row['temperature_keV'],row['density_atomic_g_cm3']
        if key in wanted and (0,.02,*key) not in records:
            uncut,_=grouped(row['source'],key)
            records[0,.02,*key]={'factor':row['rosseland_component_estimate']/uncut,'source':row['source']}
    checks=[]
    for (x,z,t,r),row in records.items():
        if x in (.025,.1) and z==.02:
            lo,hi=records[0,z,t,r]['factor'],records[.2,z,t,r]['factor']
            f=x/.2;value=math.exp((1-f)*math.log(lo)+f*math.log(hi))
            checks.append({'axis':'X','X':x,'Z':z,'temperature_keV':t,
                'density_atomic_g_cm3':r,'bracket':[0,.2],
                'direct_factor':row['factor'],'interpolated_factor':value,
                'relative_error':value/row['factor']-1})
        if x in (0,.2) and z==.02 and (x,.01,t,r) in records and (x,.03,t,r) in records:
            value=math.sqrt(records[x,.01,t,r]['factor']*records[x,.03,t,r]['factor'])
            checks.append({'axis':'Z','X':x,'Z':z,'temperature_keV':t,
                'density_atomic_g_cm3':r,'bracket':[.01,.03],
                'direct_factor':row['factor'],'interpolated_factor':value,
                'relative_error':value/row['factor']-1})
    if len(checks)!=24 or len(normalized_checks)!=48:
        raise ValueError('missing composition or normalization comparison')
    maximum=max(abs(r['relative_error']) for r in checks)
    for name,h in inputs.items():
        if digest(Path(name))!=h:
            raise ValueError('input changed during factor comparison')
    result={'scope':__doc__,'accepted_for_stellar_opacity':False,
        'composition_comparisons':len(checks),'relative_interpolation_criterion':.005,
        'maximum_relative_interpolation_error':maximum,
        'sampled_composition_interpolation_passed':maximum<=.005,
        'maximum_normalized_mean_relative_to_log_spectrum':max(abs(r['relative_to_log_spectrum_reference']) for r in normalized_checks),
        'source_normalization_is_constructed_not_independent_check':True,
        'original_full_rectangle_uncut_recovery_check_passed':reports[1]['uncut_recovery_check_passed'],
        'checks':checks,'normalized_mean_checks':normalized_checks,
        'factor_records':[{'X':k[0],'Z':k[1],'temperature_keV':k[2],
            'density_atomic_g_cm3':k[3],**v} for k,v in sorted(records.items())],
        'input_sha256':inputs,
        'limitations':['The ratio is tested at six thermal-density pairs, not across a complete factor table.',
            'Only hydrogen fractions inside 0–0.2 have independent factor-interpolation checks here.',
            'Density/temperature interpolation and a complete stellar opacity family remain to be checked.',
            'The source averaging convention is retained by construction. This does not establish convergence of the underlying atomic frequency grid.',
            'Current stellar opacity inputs are unchanged.']}
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('checks','normalized_mean_checks','factor_records','input_sha256','limitations')}),flush=True)


if __name__=='__main__':
    main()
