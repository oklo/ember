#!/usr/bin/env python3
"""Compare grey density interpolation against retained finer TOPS source values.

No sources or runtime queries are regenerated. Alternative interpolants use
only the coarse native nodes. The finer nodes and separate reference densities
remain validation data. A composition-correction experiment withholds X=0.3.
"""
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path

import numpy as np
from scipy.interpolate import Akima1DInterpolator, CubicSpline, PchipInterpolator
from audit_refractive_opacity_family import read_table
from fetch_tops_composition import digest
from reduce_tops_group_factors import add_inputs, verify


def local_polynomial(x,y,q,order):
    i=int(np.searchsorted(x,q,side='right'))-1
    start=max(0,min(i-(order//2-1),len(x)-order))
    xx=x[start:start+order];yy=y[start:start+order]
    terms=[]
    for j in range(order):
        w=math.prod((q-xx[k])/(xx[j]-xx[k]) for k in range(order) if k!=j)
        terms.append(w*yy[j])
    return math.fsum(terms)


def interpolants(x,y):
    return {'local_cubic':lambda q:local_polynomial(x,y,q,4),
            'local_quintic':lambda q:local_polynomial(x,y,q,6),
            'pchip':PchipInterpolator(x,y,extrapolate=False),
            'akima':Akima1DInterpolator(x,y,extrapolate=False),
            'makima':Akima1DInterpolator(x,y,method='makima',extrapolate=False),
            'cubic_spline':CubicSpline(x,y,extrapolate=False)}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    assembly_path=Path('docs/results/tops_bridge_v2.json');assembly=json.loads(assembly_path.read_text())
    inputs=dict(assembly['input_sha256']);add_inputs(inputs,assembly['output_sha256'])
    for path in (Path(__file__),assembly_path,Path('scripts/audit_refractive_opacity_family.py')):
        add_inputs(inputs,{str(path.resolve()):digest(path)})
    source_path=next(Path(p) for p in assembly['output_sha256'] if Path(p).name=='assembled_states.jsonl')
    groups=defaultdict(list)
    for line in source_path.read_text().splitlines():
        r=json.loads(line)
        if r['X'] in (.1,.3,.7) and r['Z']==.02 and .008<=r['temperature_keV']<=.02:
            groups[r['X'],r['temperature_keV']].append(r)
    coarse={};curves={};fine={};fine_curves={}
    historical={}
    for x in (.1,.3,.7):
        path=Path(f'docs/results/tops_cool_grey_refined_x{round(100*x):03d}-z020_v1.json')
        r=json.loads(path.read_text())
        add_inputs(inputs,{str(path.resolve()):digest(path)})
        add_inputs(historical,r['input_sha256']);add_inputs(inputs,r['output_sha256'])
        table=read_table(Path(next(iter(r['output_sha256']))))
        for i,t in enumerate(r['temperatures_keV']):
            if (x,t) not in groups:continue
            rows=sorted(groups[x,t],key=lambda r:r['density_atomic_g_cm3'])
            xx=np.log10([r['density_atomic_g_cm3'] for r in rows])
            yy=np.log10([r['uncut_rosseland'] for r in rows])
            coarse[x,t]=(xx,yy);curves[x,t]=interpolants(xx,yy)
            fxx=np.array(table['logrho']);fyy=np.array(table['planes'][x][i])
            fine[x,t]=(fxx,fyy)
            fine_curves[x,t]=interpolants(fxx,fyy)
            for r0 in rows:
                rho=r0['density_atomic_g_cm3']
                if rho not in r['densities_atomic_g_cm3']:continue
                j=r['densities_atomic_g_cm3'].index(rho)
                if abs(10**fyy[j]/r0['uncut_rosseland']-1)>1e-12:raise ValueError('coarse/fine source overlap changed')
    for relocation in assembly['source_relocations']:
        key=relocation['original']
        if key in historical:
            if historical[key]!=relocation['sha256']:raise ValueError('historical source versions disagree')
            del historical[key];historical[relocation['retained']]=relocation['sha256']
    add_inputs(inputs,historical);verify(inputs)
    extra=[]
    for (x,t),(fxx,fyy) in sorted(fine.items()):
        coarse_r={round(v,12) for v in coarse[x,t][0]}
        for q,truth in zip(fxx,fyy):
            if round(q,12) in coarse_r:continue
            predictions={name:10**(float(fn(q))-truth)-1 for name,fn in curves[x,t].items()}
            if not all(math.isfinite(v) for v in predictions.values()):raise ValueError('invalid interpolation')
            extra.append(dict(X=x,T_keV=t,rho=10**q,relative_errors=predictions))
    previous_path=Path('docs/results/tops_bridge_runtime_v1.json')
    ratios_path=Path('docs/results/opacity_bridge_existing_ratios_v1.json')
    for path in (previous_path,ratios_path):add_inputs(inputs,{str(path.resolve()):digest(path)})
    prior=json.loads(previous_path.read_text());ratios=json.loads(ratios_path.read_text())
    ratio_map={(r['X'],r['temperature_keV'],r['density_atomic_g_cm3']):r['interpolated_ratio'] for r in ratios['comparisons']}
    controls=[]
    for r in prior['independent_comparisons']:
        key=r['X'],r['temperature_keV']
        if key not in curves or not .01<=r['density_atomic_g_cm3']<=1:continue
        q=math.log10(r['density_atomic_g_cm3']);truth=math.log10(r['native_uncut'])
        predictions={name:10**(float(fn(q))-truth)-1 for name,fn in curves[key].items()}
        finer=float(fine_curves[key]['local_cubic'](q))
        row=dict(X=key[0],T_keV=key[1],rho=r['density_atomic_g_cm3'],source=r['source'],
            coarse_grey_relative_errors=predictions,fine_grey_relative_error=10**(finer-truth)-1)
        factor=ratio_map.get((*key,r['density_atomic_g_cm3']))
        if factor is not None:row['fine_grey_times_retained_ratio_relative_error']=10**finer*factor/r['normalized_mean']-1
        controls.append(row)
    # Test whether sharing the density correction between distant hydrogen
    # abundances would be justified. X=0.3 fine nodes are withheld completely.
    correction=[]
    for (x,t),(fxx,fyy) in sorted(fine.items()):
        if x!=.3:continue
        frac=(.3-.1)/(.7-.1)
        for q,truth in zip(fxx,fyy):
            if not -2<=q<=0:continue
            residuals=[float(fine_curves[xx,t]['local_cubic'](q))-float(curves[xx,t]['local_cubic'](q)) for xx in (.1,.7)]
            prediction=float(curves[x,t]['local_cubic'](q))+(1-frac)*residuals[0]+frac*residuals[1]
            correction.append(dict(X=x,T_keV=t,rho=10**q,relative_error=10**(prediction-truth)-1))
    summaries=[]
    for label,select in [('full_fine_grid',lambda r:True),('nearby_rho_001_1',lambda r:.01-1e-13<=r['rho']<=1+1e-13)]:
        rows=[r for r in extra if select(r)]
        for name in next(iter(curves.values())):
            worst=max(rows,key=lambda r:abs(r['relative_errors'][name]))
            summaries.append(dict(domain=label,method=name,comparisons=len(rows),
                maximum_relative_error=abs(worst['relative_errors'][name]),worst=worst,
                passing_005=int(sum(abs(r['relative_errors'][name])<=.005 for r in rows))))
    product=[r for r in controls if 'fine_grey_times_retained_ratio_relative_error' in r]
    result=dict(scope=__doc__,outcome='completed_interpolation_diagnosis',accepted_for_stellar_opacity=False,
        new_source_requests=0,new_runtime_queries=0,extra_native_density_checks=extra,method_summaries=summaries,
        independent_controls=controls,maximum_fine_grey_control_relative_error=max(abs(r['fine_grey_relative_error']) for r in controls),
        fine_grey_times_ratio_controls=len(product),maximum_fine_grey_times_ratio_relative_error=max(abs(r['fine_grey_times_retained_ratio_relative_error']) for r in product),
        withheld_X03_density_correction=correction,maximum_withheld_X03_density_correction_error=max(abs(r['relative_error']) for r in correction),
        limitations=['Different interpolants are diagnostics; none changes the C++ implementation or selected inputs.',
            'Fine grey tables cover X=0.1,0.3,0.7 at Z=0.02; they do not form the required full composition family.',
            'The correction experiment tests one withheld hydrogen fraction, not arbitrary metallicity or all hydrogen fractions.',
            'All comparisons are at native source temperatures; thermal interpolation and physical source joins remain separate.'],input_sha256=inputs)
    verify(inputs);args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(methods=[{k:v for k,v in r.items() if k!='worst'} for r in summaries],
        maximum_fine_grey_control_relative_error=result['maximum_fine_grey_control_relative_error'],
        product_controls=len(product),maximum_fine_grey_times_ratio_relative_error=result['maximum_fine_grey_times_ratio_relative_error'],
        maximum_withheld_X03_density_correction_error=result['maximum_withheld_X03_density_correction_error'])),flush=True)


if __name__=='__main__':main()
