#!/usr/bin/env python3
"""Measure the small-temperature-gradient assumption on the saved star.

Use actual gravity Gm/r^2 and EOS pressure to compare the temperature scale
with the pressure scale height: epsilon_T = |P/(rho g) dlnT/dr|. Reuse cached
physical EOS replies and the existing fully-ionized material classification.
No diffusion velocity, thermal coefficient or accepted omission is inferred.
"""
import argparse
import csv
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import re
import numpy as np


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def intervals(radius,mass,rho,temperature,pressure,G,stride):
    result=[]
    for i in range(0,len(radius)-stride,stride):
        j=i+stride;r=.5*(radius[i]+radius[j]);m=.5*(mass[i]+mass[j])
        density=np.sqrt(rho[i]*rho[j]);p=np.sqrt(pressure[i]*pressure[j]);g=G*m/r**2
        gradient=np.log1p((temperature[j]-temperature[i])/temperature[i])/(radius[j]-radius[i])
        hp=p/(density*g)
        result.append(dict(left=i,right=j,mass_interval_g=float(mass[j]-mass[i]),
                           radius_cm=float(r),gravity_cm_s2=float(g),pressure_scale_height_cm=float(hp),
                           dlnT_dr_cm_inverse=float(gradient),epsilon_T=float(abs(hp*gradient))))
    return result


def summarize(rows,total_mass):
    if not rows:return dict(intervals=0)
    values=np.array([r['epsilon_T'] for r in rows]);weight=np.array([r['mass_interval_g'] for r in rows])
    order=np.argsort(values);cdf=np.cumsum(weight[order])/sum(weight)
    return dict(intervals=len(rows),covered_mass_fraction=float(sum(weight)/total_mass),
                minimum=float(min(values)),maximum=float(max(values)),
                mass_weighted_mean=float(np.dot(values,weight)/sum(weight)),
                mass_weighted_quantiles={str(q):float(values[order[np.searchsorted(cdf,q)]]) for q in (.1,.5,.9)},
                covered_mass_fraction_above={str(t):float(sum(weight[values>t])/total_mass) for t in (.001,.01,.03,.1)})


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('report',type=Path);a=parser.parse_args()
    profile=Path('docs/reports/2026-09-11/evolution_latest_profile.csv')
    regimes=Path('docs/results/diffusion_material_regime_3890gyr_v1.json')
    comparison=Path('docs/results/smooth_eos_refined_retained_sources_v1.json')
    work=Path('/tmp/ember-smooth-eos-refined-physical-run-v1/retained')
    inp,out=work/'retained.input',work/'retained.output'
    source=json.loads(comparison.read_text())
    for q in (inp,out):assert source['artifacts_sha256'][str(q)]==sha(q)
    records=json.loads(regimes.read_text());material={r['zone']:r for r in records['records']}
    data=list(csv.DictReader(profile.open()));n=len(data)
    vectors={k:np.array([float(row[k]) for row in data]) for k in data[0]}
    cache={}
    for q,r in zip(inp.read_text().splitlines(),out.read_text().splitlines()):
        coords=list(map(float,q.split()));reply=json.loads(r)
        if reply['ok']:cache[tuple(coords[:4])]=reply['values'][0]
    pressure=np.array([cache[(float(r['X']),float(r['Y3']),float(r['temperature_K']),float(r['density_g_cm3']))] for r in data])
    G=float(re.search(r'\bG\s*=\s*([\d.eE+-]+)',Path('include/ember/constants.hpp').read_text()).group(1))
    # An analytic hydrostatic power-law family has dlnT/dlnP = b/a.
    # These tests check the scale-height definition and midpoint quadrature.
    controls=[]
    for exponent in (-1.,-4.):
        for b in (0.,-.05,-.4):
            errors=[]
            for points in (33,65,129):
                r=np.geomspace(1.,2.,points);m=r**3;p=r**exponent;rho=-exponent*p/(G*r**2)
                values=intervals(r,m,rho,r**b,p,G,1)
                error=max(abs(v['epsilon_T']-abs(b/exponent)) for v in values)
                errors.append(float(error))
            controls.append(dict(pressure_exponent=exponent,temperature_exponent=b,
                                 expected_epsilon_T=abs(b/exponent),absolute_errors=errors,
                                 passed=bool(errors[-1]<2e-5 and (b==0 or errors[0]>3.9*errors[1]>3.9**2*errors[2]))))
    samples={};summary={};total=vectors['mass_g'][-1]
    for stride in (1,2):
        rows=intervals(vectors['radius_cm'],vectors['mass_g'],vectors['density_g_cm3'],
                       vectors['temperature_K'],pressure,G,stride)
        for r in rows:
            i,j=r['left'],r['right'];hot=all(k in material for k in range(i,j+1))
            r['all_hot']=hot;r['all_nonconvective']=hot and all(not material[k]['convective'] for k in range(i,j+1))
            r['all_T_over_TF_below_0p3']=hot and all(material[k]['temperature_over_fermi_temperature']<.3 for k in range(i,j+1))
        samples[str(stride)]=rows
        summary[str(stride)]={name:summarize([r for r in rows if select(r)],total) for name,select in [
            ('hot',lambda r:r['all_hot']),('hot_nonconvective',lambda r:r['all_nonconvective']),
            ('hot_nonconvective_T_over_TF_below_0p3',lambda r:r['all_nonconvective'] and r['all_T_over_TF_below_0p3'])]}
    passed=all(c['passed'] for c in controls)
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),scope=__doc__,outcome='completed_diagnostic' if passed else 'failed_analytic_checks',
                accepted_thermal_diffusion_omission=False,new_EOS_queries=0,new_FreeEOS_queries=0,
                gravity_definition='G times actual enclosed baryonic mass divided by radius squared; no ideal-pressure gravity inference',
                pressure_definition='Checked refined EOS total pressure, including radiation; sampled on the unchanged saved structure',
                mass_weight_definition='Mass differences between saved mass coordinates within each retained radial interval; boundary-crossing intervals omitted from subgroup summaries',
                analytic_controls=controls,summaries_by_stride=summary,intervals_by_stride=samples,
                inputs_sha256={str(q):sha(q) for q in [Path(__file__),profile,regimes,comparison,inp,out,Path('include/ember/constants.hpp')]},
                references=['https://arxiv.org/html/1710.08424v2#S3.SS2.SSS2'],
                limitations=['epsilon_T is a scale-height diagnostic, not a bound on thermal diffusion or heat-flow effects.',
                             'Electron degeneracy does not itself imply a negligible temperature gradient.',
                             'Stride-two summaries measure sensitivity to radial sampling, not an independently evolved mesh refinement.',
                             'Physical thermal coefficients, near-cancellation of driving forces and reciprocal heat transport remain required.'])
    a.report.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(outcome=report['outcome'],summaries=summary)))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
