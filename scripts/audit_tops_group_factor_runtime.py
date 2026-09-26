#!/usr/bin/env python3
"""Check a dimensionless refractive-ratio table in Ember's actual interpolator.

This checks conversion, source-node values and derivatives on one rectangular
composition plane. The ratio table is not a standalone stellar opacity.
"""
import argparse
import json
import math
from pathlib import Path

from audit_tops_spectral_means import digest
from audit_tops_spectral_runtime import evaluate
from audit_tops_electron_dispersion import KEV, KB


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('reduction','probe','probe_reference','scratch','output'):
        p.add_argument(name,type=Path)
    a=p.parse_args()
    if a.scratch.exists() or a.output.exists():
        raise FileExistsError('use new factor check outputs')
    old=json.loads(a.reduction.read_text())
    build=json.loads(a.probe_reference.read_text())
    if old['X']!=0 or old['Z']!=.02:
        raise ValueError('this first runtime check expects X=0, Z=.02')
    inputs=dict(old['input_sha256'])
    # The existing compiled-probe report pins both the executable and sources.
    for name,h in build['input_sha256'].items():
        if name in inputs and inputs[name]!=h:
            raise ValueError('probe and reduction dependencies disagree')
        inputs[name]=h
    if inputs.get(str(a.probe.resolve()))!=digest(a.probe):
        raise ValueError('unidentified runtime executable')
    for path in (a.reduction,a.probe_reference,Path(__file__),Path(__file__).with_name('audit_tops_spectral_runtime.py')):
        inputs[str(path.resolve())]=digest(path)
    for name,h in inputs.items():
        if digest(Path(name))!=h:
            raise ValueError('changed check dependency: '+name)
    tt,rr=old['temperatures_keV'],old['densities_atomic_g_cm3']
    if any(n!=len(rr) for n in old['density_prefix_sizes']):
        raise ValueError('this rectangle check does not test variable density prefixes')
    factors={(r['temperature_keV'],r['density_atomic_g_cm3']):
             r['rosseland_atomic_cm2_g']/r['uncut_recombined_rosseland'] for r in old['records']}
    if set(factors)!={(t,r) for t in tt for r in rr} or any(
            not math.isfinite(f) or f<1-1e-14 for f in factors.values()):
        raise ValueError('invalid refractive-ratio rectangle')
    a.scratch.mkdir(parents=True)
    path=a.scratch/'factor.dat'
    lines=[f'1 {len(tt)} {len(rr)} DIMENSIONLESS refractive ratio; not a standalone opacity',
           ' '.join(format(math.log10(r),'.17g') for r in rr),
           ' '.join(format(math.log10(t*KEV/KB),'.17g') for t in tt),'0 0.02']
    lines.extend(' '.join(format(math.log10(factors[t,r]),'.17g') for r in rr) for t in tt)
    path.write_text('\n'.join(lines)+'\n')
    points=[(t,r) for t in tt for r in rr]
    nodes=evaluate(a.probe,path,points,'nodes')
    node_error=max(abs(v[0]/factors[k]-1) for k,v in zip(points,nodes,strict=True))
    checks=[];queries=[];h=1e-4
    centers=[(math.sqrt(tt[it]*tt[it+1]),math.sqrt(rr[ir]*rr[ir+1]))
             for it in (0,5,10,15,20,25,30,34) for ir in (5,25,45,60)]
    for t,r in centers:
        queries.extend([(t,r),(t*math.exp(-h),r),(t*math.exp(h),r),
                        (t,r*math.exp(-h)),(t,r*math.exp(h))])
    values=evaluate(a.probe,path,queries,'derivatives')
    for i,(t,r) in enumerate(centers):
        v=values[5*i:5*i+5]
        for axis,j,k,d in [('temperature',1,2,1),('density',3,4,2)]:
            finite=(math.log(v[k][0])-math.log(v[j][0]))/(2*h)
            scaled=abs(finite-v[0][d])/max(1.,abs(finite),abs(v[0][d]))
            checks.append({'temperature_keV':t,'density_atomic_g_cm3':r,
                          'axis':axis,'runtime_derivative':v[0][d],
                          'finite_difference':finite,'scaled_error':scaled})
    maximum=max(r['scaled_error'] for r in checks)
    for name,hsh in inputs.items():
        if digest(Path(name))!=hsh:
            raise ValueError('input changed during runtime check')
    result={'scope':__doc__,'accepted_for_stellar_opacity':False,
        'source_node_queries':len(points),'maximum_source_node_relative_error':node_error,
        'derivative_comparisons':len(checks),'maximum_scaled_derivative_error':maximum,
        'runtime_numerical_checks_passed':node_error<1e-12 and maximum<1e-7,
        'factor_table':str(path.resolve()),'factor_table_sha256':digest(path),
        'source_uncut_recovery_check_passed':old['uncut_recovery_check_passed'],
        'input_sha256':inputs,'derivative_checks':checks,
        'scratch_sha256':{p.name:digest(p) for p in a.scratch.iterdir() if p.is_file()}}
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('input_sha256','scratch_sha256','derivative_checks')}),flush=True)
    if not result['runtime_numerical_checks_passed']:
        raise ValueError('runtime ratio-table check failed')


if __name__=='__main__':
    main()
