#!/usr/bin/env python3
"""Test composition interpolation at fixed source temperatures and densities.

Interpolate log refractive ratio linearly in X and Z, matching Ember's
composition rule. Candidate and independent mixtures use identical photon,
temperature and density coordinates. This separates composition error from
density interpolation; it does not accept a complete stellar opacity family.
"""
import argparse
import bisect
import json
import math
from pathlib import Path

from fetch_tops_composition import digest
from reduce_tops_group_factors import add_inputs, verify


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('plan', type=Path)
    p.add_argument('report', type=Path)
    a = p.parse_args()
    if a.report.exists():
        raise FileExistsError('preserve completed composition comparisons')
    plan = json.loads(a.plan.read_text())
    inputs = {str(path.resolve()):digest(path) for path in
              (a.plan, Path(__file__), Path('src/opacity_table.cpp'), Path('src/opacity_mixture.cpp'))}
    candidate, controls = {}, []
    for item in plan['reductions']:
        path = Path(item['report'])
        r = json.loads(path.read_text())
        add_inputs(inputs, {str(path.resolve()):digest(path)})
        add_inputs(inputs, r['input_sha256'])
        if (r['X'],r['Z']) != (item['X'],item['Z']) or r['states'] != 2 or r['photon_boundaries'] != 993:
            raise ValueError('source reduction differs from its composition plan')
        for row in r['records']:
            if row['temperature_keV'] != item['temperature_keV']:
                raise ValueError('source temperature differs')
            record = {'X':r['X'], 'Z':r['Z'], **row}
            key = r['X'],r['Z'],row['temperature_keV'],row['density_atomic_g_cm3']
            if item['role'] == 'candidate':
                if key in candidate:
                    raise ValueError('duplicate candidate source')
                candidate[key] = record
            elif item['role'] == 'control':
                if r['X'] in plan['candidate_X'] or r['Z'] in plan['candidate_Z']:
                    raise ValueError('control composition was inserted into the candidate grid')
                controls.append(record)
            else:
                raise ValueError('unknown source role')
    conditions={(r['temperature_keV'],r['density_atomic_g_cm3']) for r in controls}
    expected={(x,z,t,rho) for x in plan['candidate_X'] for z in plan['candidate_Z'] for t,rho in conditions}
    if set(candidate) != expected or len(controls) != 24 or len(conditions) != 8:
        raise ValueError('incomplete composition comparison')
    verify(inputs)
    results = []
    for ref in controls:
        x,z,t,rho = (ref[k] for k in ('X','Z','temperature_keV','density_atomic_g_cm3'))
        ix,iz = bisect.bisect_left(plan['candidate_X'],x),bisect.bisect_left(plan['candidate_Z'],z)
        if not 0 < ix < len(plan['candidate_X']) or not 0 < iz < len(plan['candidate_Z']):
            raise ValueError('comparison would extrapolate')
        x0,x1 = plan['candidate_X'][ix-1:ix+1]
        z0,z1 = plan['candidate_Z'][iz-1:iz+1]
        fx,fz = (x-x0)/(x1-x0),(z-z0)/(z1-z0)
        corner = [candidate[xv,zv,t,rho] for xv,zv in ((x0,z0),(x1,z0),(x0,z1),(x1,z1))]
        weights = [(1-fx)*(1-fz),fx*(1-fz),(1-fx)*fz,fx*fz]
        if abs(sum(weights)-1) > 1e-14 or any(w < 0 for w in weights):
            raise ValueError('invalid composition interpolation weights')
        def interpolate(field):
            return math.exp(sum(w*math.log(r[field]) for r,w in zip(corner,weights,strict=True)))
        factor = interpolate('refractive_ratio')
        grey = interpolate('uncut_source_rosseland')
        results.append({'X':x,'Z':z,'temperature_keV':t,'density_atomic_g_cm3':rho,
                        'reference_source':ref['source'],
                        'candidate_sources':[r['source'] for r in corner],
                        'candidate_X_bracket':[x0,x1],'candidate_Z_bracket':[z0,z1],
                        'composition_weights':weights,
                        'direct_ratio':ref['refractive_ratio'],'interpolated_ratio':factor,
                        'ratio_relative_error':factor/ref['refractive_ratio']-1,
                        'grey_relative_error_on_pilot_composition_grid':grey/ref['uncut_source_rosseland']-1,
                        'product_relative_error_on_pilot_composition_grid':grey*factor/ref['native_uncut_times_ratio_atomic_cm2_g']-1})
    criterion = plan['relative_criterion']
    summaries = {}
    for label,rows in [('all',results),*[(f'X_{x:g}',[r for r in results if r['X']==x]) for x in plan['control_X']]]:
        summaries[label] = {'comparisons':len(rows), **{field:{
            'maximum_absolute_relative_error':max(abs(r[field]) for r in rows),
            'failures_above_criterion':sum(abs(r[field])>criterion for r in rows)} for field in
            ('ratio_relative_error','grey_relative_error_on_pilot_composition_grid',
             'product_relative_error_on_pilot_composition_grid')}}
    verify(inputs)
    result = {'scope':__doc__,'analysis_complete':True,'accepted_for_stellar_opacity':False,
              'sampled_composition_check_passed':summaries['all']['ratio_relative_error']['failures_above_criterion']==0,
              'relative_criterion':criterion,'summary':summaries,'records':results,'input_sha256':inputs,
              'limitations':[
                  'These are native-temperature, fixed-density source comparisons, not a resampled family or new trajectory.',
                  'The selected source temperatures are below the conduction activation threshold.',
                  'The pilot tests a wider grey composition spacing than the existing native grey family; its grey errors do not measure that family.',
                  'Frequency quadrature convergence is distinct from full-spectrum accuracy.',
                  'Passing this pilot does not establish interpolation over unsampled temperatures, densities or compositions.']}
    a.report.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'passed':result['sampled_composition_check_passed'],'summary':summaries}),flush=True)


if __name__ == '__main__':
    main()
