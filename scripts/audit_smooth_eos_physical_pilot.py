#!/usr/bin/env python3
"""Compare the native smooth EOS with six retained independent source targets.

Also query the complete saved stellar profile and finite-difference its
thermodynamic responses. No FreeEOS query is repeated. This is a physical
pilot, not acceptance of the full composition/temperature/density domain.
"""
import argparse
import csv
from datetime import datetime,timezone
import json
from pathlib import Path
import numpy as np
from composition_potential_v2 import RGAS,isotope_value
from metal_eos_composition import mixture
from audit_diffusion_eos_two_compositions import normalized_error
from audit_smooth_eos_runtime import query,sha,ARAD


def source_state(raw,x,y,T,rho):
    v=np.array(list(map(float,raw['stdout'].split())))
    scale=mixture(x,y)['source_mass_scale'];pr=ARAD*T**4/3
    P=v[4]+pr;E=v[5]*scale+3*pr/rho;cv=v[10]*scale/T+12*pr/(rho*T)
    chiT=(v[4]*v[8]+4*pr)/P;chiR=v[4]*v[7]/P
    delta=chiT/chiR;cp=cv+P/(rho*T)*chiT*delta;ad=P/(rho*T)*delta/cp
    S=v[6]*scale-isotope_value(x,y)+4*pr/(rho*T)
    return np.array([P,E,S,cv,cp,chiT,chiR,delta,ad,chiR*cp/cv])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('probe',type=Path)
    p.add_argument('family',type=Path);p.add_argument('work',type=Path);p.add_argument('report',type=Path);a=p.parse_args()
    a.work.mkdir()
    reference_path=Path('docs/results/diffusion_eos_two_compositions_v2.json')
    response_path=Path('docs/results/diffusion_eos_composition_responses_v1.json')
    reference=json.loads(reference_path.read_text());responses=json.loads(response_path.read_text())
    direct_path=Path('/tmp/ember-diffusion-eos-two-compositions-run-v1/direct_source.jsonl')
    # Verify the immutable direct queries against the comparison that accepted
    # their quadrature and composition-spacing accuracy.
    if sha(direct_path)!=reference['inputs_sha256'][str(direct_path)]:raise ValueError('changed direct source')
    direct={tuple(v['key']):v for v in map(json.loads,direct_path.read_text().splitlines())}
    controls=[]
    for c in reference['comparisons']:
        m=reference['material_coordinates'][c['material']]
        controls.append([c['X'],c['Y3'],m['T'],m['rho'],1])
    profile_path=Path('docs/reports/2026-09-11/evolution_latest_profile.csv')
    with profile_path.open() as f:profile=list(csv.DictReader(f))
    h=1e-5;rows=controls.copy()
    for r in profile:
        x,y,T,rho=[float(r[k]) for k in ['X','Y3','temperature_K','density_g_cm3']]
        for dt,dr in [(0,0),(h,0),(-h,0),(0,h),(0,-h)]:
            rows.append([x,y,T*np.exp(dt),rho*np.exp(dr),0])
    actual=query(a.probe,a.family,rows,a.work,'physical')
    comparisons=[]
    for i,(c,q) in enumerate(zip(reference['comparisons'],controls)):
        r=actual[i]
        if not r['ok']:
            comparisons.append(dict(query=q,passed=False,error=r['error']));continue
        x,y,T,rho=q[:4];v=np.array(r['values']);raw=direct[1,c['material'],x,y]
        target=source_state(raw,x,y,T,rho);e=v[:10]/target-1
        href=np.array(c['source_Hessian_over_Rgas'])*RGAS
        herror=normalized_error(v[24:28].reshape(2,2)-href,href)
        gradient=float(np.max(np.abs(v[22:24]/RGAS-c['source_gradient_over_Rgas'])))
        dref=np.array(responses['comparisons'][i]['source_derivatives_rows_P_E_columns_X_Y3'])
        derror=v[17:21].reshape(2,2)/dref-1
        passed=abs(e[0])<.001 and abs(e[1])<.002 and max(abs(e[3:]))<.003 and herror<.005 and gradient<1e-4 and max(abs(derror.flatten()))<.005
        comparisons.append(dict(query=q,relative_state_errors=e.tolist(),normalized_Hessian_error=herror,
                           gradient_error_over_Rgas=gradient,relative_composition_response_errors=derror.tolist(),passed=bool(passed)))
    failures=[];maximum=dict(first_law=0.,energy_temperature=0.,heat_capacity_response=0.,delta_response=0.,adiabatic_gradient_response=0.)
    for i,row in enumerate(profile):
        records=actual[len(controls)+5*i:len(controls)+5*i+5]
        if not all(r['ok'] for r in records):
            failures.append(dict(zone=i,responses=records));continue
        base,tp,tm,rp,rm=[np.array(r['values']) for r in records]
        T=float(row['temperature_K']);rho=float(row['density_g_cm3'])
        dt=(tp-tm)/(2*h);dr=(rp-rm)/(2*h)
        defects=dict(first_law=max(abs(dr[1]-base[0]/rho*(1-base[5]))/(T*base[3]),
                                  abs(dr[2]+base[0]*base[5]/(rho*T))/base[3]),
                     energy_temperature=abs(dt[1]/(T*base[3])-1),
                     heat_capacity_response=max(abs(dt[4]-base[11]),abs(dr[4]-base[12]))/base[4],
                     delta_response=max(abs(dt[7]-base[13]),abs(dr[7]-base[14]))/base[7],
                     adiabatic_gradient_response=max(abs(dt[8]-base[15]),abs(dr[8]-base[16]))/base[8])
        for k,v in defects.items():maximum[k]=max(maximum[k],float(v))
        if max(defects.values())>1e-5:failures.append(dict(zone=i,defects=defects))
    passed=all(c['passed'] for c in comparisons) and not failures
    report=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed' if passed else 'failed',
                accepted_for_evolution=False,new_FreeEOS_queries=0,independent_comparisons=comparisons,
                profile_zones=len(profile),profile_failed_checks=failures,profile_maximum_errors=maximum,
                criteria=dict(pressure=.001,energy=.002,thermal=.003,Hessian=.005,gradient_over_Rgas=1e-4,
                              composition_response=.005,profile_identities_and_responses=1e-5),
                inputs_sha256={str(p.resolve()):sha(p) for p in [Path(__file__),a.probe,a.family,reference_path,response_path,direct_path,profile_path]},
                artifacts={str(p):sha(p) for p in a.work.rglob('*') if p.is_file()},
                remaining_work=['Broader independent source coverage including ionization and dense extension',
                                'Source-mask and composition-boundary physical comparisons','Atmosphere checks',
                                'Runtime cost and memory','Acceptance and coupling to diffusion/evolution'])
    with a.report.open('x') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'outcome':report['outcome'],'profile_failures':len(failures),
                      'source_passes':sum(c['passed'] for c in comparisons),'profile_maximum_errors':maximum}))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
