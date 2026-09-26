#!/usr/bin/env python3
"""Extend only missing strong repulsive-ion collisions for the saved star."""
import hashlib,json,math,re
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from scipy.interpolate import CubicSpline
from ion_collision_integrals import E2,KB
from yukawa_scattering_control_v2 import collision_integrals
from ion_collision_table import ExtendedIonTable


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    work=Path('/tmp/ember-ion-table-extension-run-v1')
    old_path=Path('docs/results/yukawa_collision_table_v1.json')
    heat_path=Path('docs/results/stellar_electron_heat_3890gyr_v1.json')
    header=Path('include/ember/gs98_mixture.hpp')
    old=json.loads(old_path.read_text());heat=json.loads(heat_path.read_text())
    inputs={str(p):sha(p) for p in [Path(__file__),Path('scripts/ion_collision_table.py'),
        Path('scripts/yukawa_scattering_control_v2.py'),Path('scripts/yukawa_scattering_control.py'),
        Path('scripts/ion_collision_integrals.py'),old_path,heat_path,header]}
    for source in [old,heat]:
        for key in ['input_sha256','inputs_sha256']:
            for p,h in source.get(key,{}).items():
                assert sha(p)==h,p
                inputs[p]=h
    base_path=Path(old['table']);assert sha(base_path)==old['table_sha256']
    base=json.loads(base_path.read_text());x=np.array(base['log10_strength']);y=np.array(base['dimensionless_integrals'])
    assert old['status']=='pass' and base['numerical_status']=='pass'
    inputs[str(base_path)]=sha(base_path)
    metals=[tuple(map(float,s.split(','))) for s in re.findall(r'^\s*\{([^{}]+)\}, //',header.read_text(),re.M)]
    charges=np.array([1.,2.,2.]+[m[0] for m in metals])
    bounds=[]
    for row in heat['records']:
        for case in row['cases']:
            g=np.outer(charges,charges)*E2/(KB*row['temperature_K']*case['screening_length_cm'])
            bounds.append(dict(zone=row['zone'],screening=case['screening'],minimum=float(g.min()),maximum=float(g.max())))
    maximum=max(v['maximum'] for v in bounds)
    step=x[-1]-x[-2]
    count=math.ceil((math.log10(maximum)-x[-1])/step)
    if count<=0 or count>16:raise ValueError('unexpected extension scope')
    xn=x[-1]+step*np.arange(1,count+1)
    pairs=[tuple(v) for v in base['moment_pairs']]
    cache={}
    def values(logg,order=32,angle=128,radius=256,tail=8.):
        key=(float(logg),order,angle,radius,tail)
        if key not in cache:
            result=collision_integrals(10**logg,thermal_order=order,angle_order=angle,radius_order=radius,thermal_max=tail)
            cache[key]=[result[k] for k in pairs]
        return np.array(cache[key])
    yn=np.array([values(v) for v in xn])
    spec=dict(scope=__doc__,base_table=str(base_path),base_table_sha256=sha(base_path),
        additional_log10_strength=xn.tolist(),additional_dimensionless_integrals=yn.tolist(),
        interpolation='Original cubic log-moment interpolant through its existing upper limit, then C2 quintic Hermite extension; no extrapolation.',
        collision_strength_min=float(10**x[0]),collision_strength_max=float(10**xn[-1]),
        numerical_status='pending',accepted_for_stellar_evolution=False)
    table_path=work/'extended_collision_integrals.json'
    with table_path.open('x') as f:json.dump(spec,f,indent=2,allow_nan=False);f.write('\n')
    table=ExtendedIonTable(table_path)
    old_spline=CubicSpline(x,np.log(y),axis=0,extrapolate=False)
    queries=np.r_[x,[v['log10_strength'] for v in old['interpolation_controls']]]
    exact_reuse=np.array_equal(table.log_moments(queries),old_spline(queries))
    joins=[]
    for order in range(3):
        left=table.base(table.join,order);right=table.extension(table.join,order)
        joins.append(float(np.max(abs(left-right)/(1+abs(left)))))
    mids=(np.r_[x[-1],xn][:-1]+xn)/2
    controls=[]
    for logg in mids:
        actual=values(logg);predicted=table.moments(10**logg)
        controls.append(dict(log10_strength=float(logg),reference=actual.tolist(),
                             predicted=predicted.tolist(),relative_error=(predicted/actual-1).tolist()))
    quadrature=[]
    for logg in [mids[0],xn[-1]]:
        coarse=values(logg);fine=values(logg,64,192,384);tail=values(logg,32,128,256,10.)
        quadrature.append(dict(log10_strength=float(logg),coarse=coarse.tolist(),fine=fine.tolist(),
                              relative_quadrature_change=(coarse/fine-1).tolist(),relative_tail_change=(tail/coarse-1).tolist()))
    qe=max(abs(e) for c in quadrature for e in c['relative_quadrature_change'])
    te=max(abs(e) for c in quadrature for e in c['relative_tail_change'])
    ie=max(abs(e) for c in controls for e in c['relative_error'])
    unsupported_rejected=bool(np.all(np.isnan(table.log_moments([x[0]-.01,xn[-1]+.01]))))
    passed=exact_reuse and max(joins)<1e-10 and qe<2e-5 and te<2e-5 and ie<1e-3 and unsupported_rejected
    spec['numerical_status']='pass' if passed else 'fail'
    # This source specification is not published/hashed until all checks finish.
    table_path.write_text(json.dumps(spec,indent=2,allow_nan=False)+'\n')
    sources=[dict(log10_strength=k[0],thermal_order=k[1],angle_order=k[2],radius_order=k[3],thermal_max=k[4],moments=v) for k,v in cache.items()]
    for p,h in inputs.items():assert sha(p)==h,p
    report=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),
        status='pass' if passed else 'fail',accepted_for_stellar_evolution=False,
        table=str(table_path),table_sha256=sha(table_path),base_nodes_reused=len(x),additional_nodes=len(xn),
        source_calculations=len(cache),saved_layer_screen_pairs=len(bounds),
        required_collision_strength_max=maximum,table_collision_strength_max=float(10**xn[-1]),
        old_queries_byte_identical=exact_reuse,old_query_checks=len(queries),
        continuity_relative_errors_through_second_derivative=joins,
        maximum_quadrature_change=qe,maximum_tail_change=te,maximum_interpolation_error=ie,
        unsupported_queries_rejected=unsupported_rejected,quadrature_controls=quadrature,
        interpolation_controls=controls,source_calculations_retained=sources,
        input_sha256=inputs,limitations=['Numerical extension of the same classical repulsive Yukawa potential only.',
        'No electron scattering, quantum ion correction or new interacting-mixture prescription.',
        'Metal charges remain prescribed as fully stripped in the support scan.'])
    with Path('docs/results/ion_collision_table_extension_v1.json').open('x') as f:
        json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:report[k] for k in ['status','additional_nodes','source_calculations',
        'required_collision_strength_max','table_collision_strength_max','old_queries_byte_identical',
        'maximum_quadrature_change','maximum_tail_change','maximum_interpolation_error']}))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
