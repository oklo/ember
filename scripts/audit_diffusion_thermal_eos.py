#!/usr/bin/env python3
"""Check material exchange enthalpies with native EOS and retained FreeEOS data.

One persistent native process supplies only previously unqueried neighboring
states for independent fixed-pressure finite differences. Saved full replies
and physical composition derivatives are reused. This validates sampled EOS
identities, not kinetic coefficients or stellar diffusion velocities.
"""
import argparse
from datetime import datetime, timezone
import hashlib, json, subprocess
from pathlib import Path
import numpy as np
from audit_smooth_eos_runtime import ARAD
from composition_potential_v2 import RGAS
from diffusion_thermal_transport import exchange_enthalpy


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['probe','family','work','report']:p.add_argument(name,type=Path)
    a=p.parse_args();a.work.mkdir()
    checked=[];cache={};reused=0;new=0
    reference_path=Path('docs/results/diffusion_eos_two_compositions_v2.json')
    response_path=Path('docs/results/diffusion_eos_composition_responses_v1.json')
    reference=json.loads(reference_path.read_text());responses=json.loads(response_path.read_text())
    direct_path=Path('/tmp/ember-diffusion-eos-two-compositions-run-v1/direct_source.jsonl')
    assert sha(direct_path)==reference['inputs_sha256'][str(direct_path)]
    assert sha(direct_path)==responses['inputs_sha256'][str(direct_path)]
    direct={tuple(v['key']):v for v in map(json.loads,direct_path.read_text().splitlines())}
    reports=[Path('docs/results/smooth_eos_refined_chemical_domain_v1.json'),
             Path('docs/results/smooth_eos_refined_physical_profile_v1.json')]
    for report in reports:
        d=json.loads(report.read_text());assert d['outcome']=='passed'
        artifacts=d.get('artifacts_sha256',d.get('artifacts',{}))
        for path in artifacts:
            if not path.endswith('.input'):continue
            inp=Path(path);out=inp.with_suffix('.output')
            assert sha(inp)==artifacts[str(inp)] and sha(out)==artifacts[str(out)]
            checked.extend([inp,out])
            inputs=inp.read_text().splitlines();outputs=out.read_text().splitlines()
            assert len(inputs)==len(outputs)
            for row,line in zip(inputs,outputs):
                key=tuple(map(float,row.split()));v=json.loads(line)
                if v['ok']:cache[key]=np.array(v['values'])
    # The active-API build was separately checked to retain all full replies.
    regression_path=Path('docs/results/smooth_eos_active_regression_v2.json')
    regression=json.loads(regression_path.read_text());assert regression['outcome']=='passed'
    assert sha(a.probe)==regression['inputs_sha256'][str(a.probe)]
    expected=json.loads(reports[-1].read_text())['inputs_sha256']
    assert sha(a.family)==expected[str(a.family.resolve())]
    log=(a.work/'native.stderr').open('x')
    input_log=(a.work/'new.input').open('x');output_log=(a.work/'new.output').open('x')
    proc=subprocess.Popen([str(a.probe),str(a.family),'smooth'],stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,stderr=log,text=True,bufsize=1)
    def query(x,y,T,rho,chemical=0):
        nonlocal reused,new
        key=(float(x),float(y),float(T),float(rho),float(chemical))
        # A full chemical reply also includes every thermal-only response.
        full=key[:4]+(1.,)
        if key in cache:reused+=1;return cache[key]
        if not chemical and full in cache:reused+=1;return cache[full]
        text=' '.join(format(v,'.17g') for v in key)+'\n'
        proc.stdin.write(text);proc.stdin.flush();input_log.write(text);input_log.flush()
        line=proc.stdout.readline()
        if not line:raise RuntimeError('native EOS exited before replying')
        output_log.write(line);output_log.flush();new+=1
        reply=json.loads(line)
        if not reply['ok']:raise ValueError('Unsupported EOS query '+text.strip()+': '+reply['error'])
        v=np.array(reply['values']);assert np.all(np.isfinite(v));cache[key]=v
        return v
    pressure_residuals=[]
    def at_pressure(x,y,T,target,rho_guess):
        logrho=np.log(rho_guess)
        for iteration in range(15):
            rho=np.exp(logrho);v=query(x,y,T,rho);pm=v[0]-ARAD*T**4/3
            residual=pm/target-1
            if abs(residual)<2e-13:
                pressure_residuals.append(abs(residual));return rho,v
            correction=(pm-target)/(v[0]*v[6])
            if not np.isfinite(correction) or abs(correction)>.1:
                raise ValueError('constant-pressure density correction outside local check')
            logrho-=correction
        raise ValueError('constant-pressure EOS solve failed')
    def material_enthalpy(v,T,rho):
        pr=ARAD*T**4/3
        return v[1]-3*pr/rho+(v[0]-pr)/rho
    cases=[];checks=[]
    def check(name,actual,expected,scales,tolerance,**meta):
        error=float(np.max(abs(np.asarray(actual)-expected)/scales))
        checks.append(dict(name=name,normalized_error=error,tolerance=tolerance,
                           passed=bool(np.isfinite(error) and error<=tolerance),**meta))
    try:
        for index,(c,response) in enumerate(zip(reference['comparisons'],responses['comparisons'])):
            material=reference['material_coordinates'][c['material']]
            x,y,T,rho=c['X'],c['Y3'],material['T'],material['rho']
            base=query(x,y,T,rho,1);prad=ARAD*T**4/3;pm=base[0]-prad
            delta=(base[0]*base[5]-4*prad)/(base[0]*base[6])
            h=exchange_enthalpy(temperature=T,gradient_phi_lnT=base[28:30],
                               gradient_phi_lnrho=base[30:32],material_delta=delta)
            scales=RGAS*T+abs(h);meta=dict(control=index,material=material['name'],X=x,Y3=y,T=T,rho=rho)
            check('native pressure-energy derivative identity',base[19:21]+delta*base[17:19]/rho,h,scales,2e-11,**meta)
            source=direct[1,c['material'],x,y]
            raw=np.array(list(map(float,source['stdout'].split())))
            source_delta=raw[8]/raw[7]
            source_derivatives=np.array(response['source_derivatives_rows_P_E_columns_X_Y3'])
            source_h=source_derivatives[1]+source_delta*source_derivatives[0]/rho
            check('independent FreeEOS exchange enthalpy',h,source_h,abs(source_h),.005,**meta)
            spacing=[]
            for factor in [1.,.5]:
                measured=[]
                for direction,fraction in enumerate([x,y]):
                    step=factor*min(1e-4,.1*fraction,.1*(.98-x-y))
                    values=[]
                    for offset in [-2,-1,1,2]:
                        xx,yy=[x,y];
                        if direction==0:xx+=offset*step
                        else:yy+=offset*step
                        rr,v=at_pressure(xx,yy,T,pm,rho)
                        values.append(material_enthalpy(v,T,rr))
                    measured.append((values[0]-8*values[1]+8*values[2]-values[3])/(12*step))
                spacing.append(np.array(measured))
            check('independent constant-pressure bulk enthalpy',spacing[-1],h,scales,2e-6,**meta)
            check('composition finite-difference spacing',spacing[0],spacing[1],scales,2e-6,**meta)
            thermal=[]
            for step in [1e-4,5e-5]:
                values=[]
                for sign in [1,-1]:
                    tt=T*np.exp(sign*step);rr,v=at_pressure(x,y,tt,pm,rho)
                    chemical=query(x,y,tt,rr,1);values.append(tt*chemical[22:24])
                thermal.append(T*base[22:24]-(values[0]-values[1])/(2*step))
            check('independent constant-pressure chemical derivative',thermal[-1],h,scales,2e-6,**meta)
            check('temperature finite-difference spacing',thermal[0],thermal[1],scales,2e-6,**meta)
            cases.append(dict(**meta,material_delta=float(delta),exchange_enthalpies_erg_g=h.tolist(),
                              source_enthalpies_erg_g=source_h.tolist(),normalization_scales_erg_g=scales.tolist()))
    finally:
        proc.stdin.close()
        try:code=proc.wait(timeout=10)
        except subprocess.TimeoutExpired:proc.kill();proc.wait();raise
        log.close();input_log.close();output_log.close()
    assert code==0
    passed=all(c['passed'] for c in checks)
    inputs=[Path(__file__),Path('scripts/diffusion_thermal_transport.py'),a.probe,a.family,reference_path,
            response_path,direct_path,regression_path,*reports,*checked]
    report=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed' if passed else 'failed',
                accepted_for_stellar_evolution=False,new_FreeEOS_queries=0,new_native_queries=new,reused_native_replies=reused,
                cases=cases,checks=checks,constant_pressure_solutions=len(pressure_residuals),
                maximum_relative_pressure_residual=max(pressure_residuals),
                inputs_sha256={str(p.resolve()):sha(p) for p in inputs},
                artifacts_sha256={str(p):sha(p) for p in a.work.iterdir() if p.is_file()},
                limitations=['Six sampled physical controls, not full EOS-domain acceptance',
                             'Kinetic thermal coefficients are not supplied',
                             'No stellar velocity, flux or abundance evolution'])
    with a.report.open('x') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(outcome=report['outcome'],controls=len(cases),checks=len(checks),
                         new_native_queries=new,reused_native_replies=reused,failures=[c for c in checks if not c['passed']])))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
