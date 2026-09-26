#!/usr/bin/env python3
"""Check the pressure and energy composition derivatives for the smooth EOS pilot.

No source or runtime query is repeated. Composition derivatives of pressure
and energy give the density/composition and temperature/composition responses
of Phi=F/T through Phi_(lnrho,X)=P_X/(rho*T) and Phi_(lnT,X)=-E_X/T.
This tests these responses at six saved comparison coordinates; full runtime
EOS, heat-capacity derivatives, domain boundaries and transport remain pending.
"""
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root=Path(__file__).resolve().parents[1]
    report=root/'docs/results/diffusion_eos_two_compositions_v2.json'
    ref=json.loads(report.read_text())
    oldpath=root/'docs/results/diffusion_eos_two_compositions_v1.json';old=json.loads(oldpath.read_text())
    for data in [old,ref]:
        for p,h in data['inputs_sha256'].items():assert digest(p)==h
    for p,h in old['output_sha256'].items():assert digest(p)==h
    work=Path('/tmp/ember-diffusion-eos-two-compositions-run-v1')
    source={tuple(v['key']):v for v in map(json.loads,(work/'direct_source.jsonl').read_text().splitlines())}
    cache={tuple(v[:4]):v for name in [str(work/'runtime_stdout.txt'),'/tmp/ember-diffusion-eos-spline-run-v1/stdout.txt'] for v in np.loadtxt(name)}
    family=root/'data/eos/numerical_electron_base_v1/freeeos300_gs98_z020.dat'
    xs=np.array([float(v) for v in family.read_text().splitlines()[1].split()[2:]]);ys=np.array([0.,.06,.12])
    bx=CubicSpline(xs,np.eye(len(xs)),axis=0,extrapolate=False)
    by=CubicSpline(ys,np.eye(len(ys)),axis=0,extrapolate=False)
    grids={};arad=4*5.670374419e-5/2.99792458e10
    for mi,m in enumerate(ref['material_coordinates']):
        grid=np.empty((len(xs),len(ys),2));t,rho=m['T'],m['rho'];pr=arad*t**4/3
        for i,x in enumerate(xs):
            for j,y in enumerate(ys):
                if y==.06:
                    v=source[0,mi,x,y];grid[i,j]=[v['P']+pr,v['E']+3*pr/rho]
                else:grid[i,j]=cache[x,y,t,rho][5:7]
        grids[mi]=grid
    cases=[]
    for c in ref['comparisons']:
        mi,x,y,hx,hy=c['material'],c['X'],c['Y3'],c['hx'],c['hy']
        m=ref['material_coordinates'][mi];t,rho=m['T'],m['rho']
        grid=grids[mi];w=bx(x);wx=bx(x,1);v=by(y);vy=by(y,1)
        value=np.array([w@grid[:,:,k]@v for k in range(2)])
        grad=np.array([[wx@grid[:,:,k]@v,w@grid[:,:,k]@vy] for k in range(2)])
        direct=[];fine=[]
        for bi in range(2):
            def d(s):
                rows=[]
                for k in ['P','E']:
                    dx=(source[bi,mi,x+s*hx,y][k]-source[bi,mi,x-s*hx,y][k])/(2*s*hx)
                    dy=(source[bi,mi,x,y+s*hy][k]-source[bi,mi,x,y-s*hy][k])/(2*s*hy)
                    rows.append([dx,dy])
                return np.array(rows)
            a,b=d(1),d(.5);direct.append((4*b-a)/3);fine.append(b)
        assert np.all(np.abs(direct[1])>0)
        precision=float(np.max(np.abs(direct[0]/direct[1]-1)))
        spacing=float(np.max(np.abs(fine[1]/direct[1]-1)))
        assert precision<1e-4 and spacing<1e-4,(c,precision,spacing)
        raw=source[1,mi,x,y];pr=arad*t**4/3
        target=np.array([raw['P']+pr,raw['E']+3*pr/rho])
        valerr=value/target-1;derr=grad/direct[1]-1
        cases.append({'material':c['material_name'],'X':x,'Y3':y,'T':t,'rho':rho,
                      'relative_value_errors_P_E':valerr.tolist(),
                      'relative_derivative_errors_rows_P_E_columns_X_Y3':derr.tolist(),
                      'source_derivatives_rows_P_E_columns_X_Y3':direct[1].tolist(),
                      'candidate_derivatives_rows_P_E_columns_X_Y3':grad.tolist(),
                      'candidate_Phi_lnrho_composition':(grad[0]/(rho*t)).tolist(),
                      'candidate_Phi_lnT_composition':(-grad[1]/t).tolist(),
                      'source_precision_error':precision,'source_spacing_error':spacing,
                      'passed':bool(np.max(np.abs(valerr))<.003 and np.max(np.abs(derr))<.005)})
    output=root/'docs/results/diffusion_eos_composition_responses_v1.json'
    result={'scope':__doc__,'created_utc':datetime.now(timezone.utc).isoformat(),
            'outcome':'completed_response_pilot','all_comparisons_pass':all(c['passed'] for c in cases),
            'accepted_for_evolution':False,'new_source_queries':0,'new_runtime_queries':0,
            'criteria':{'relative_value':.003,'relative_composition_derivative':.005,
                        'source_precision_and_spacing':1e-4},'comparisons':cases,
            'inputs_sha256':{str(p):digest(p) for p in [Path(__file__),report,oldpath,family,
                work/'direct_source.jsonl',work/'runtime_stdout.txt']}}
    with output.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({'all_pass':result['all_comparisons_pass'],
                      'max_value_error':max(abs(v) for c in cases for v in c['relative_value_errors_P_E']),
                      'max_derivative_error':max(abs(v) for c in cases for row in c['relative_derivative_errors_rows_P_E_columns_X_Y3'] for v in row)}))


if __name__=='__main__':main()
