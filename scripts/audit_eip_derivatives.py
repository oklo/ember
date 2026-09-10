#!/usr/bin/env python3
"""Compare EOSFI22 derivatives with independent finite differences of its potential.

The forced-phase source probe excludes ideal electrons/radiation and evaluates
the exact requested density. The audit does not select a stable phase or install
an EOS. Normalization uses ideal-ion energy/pressure, not small cancellation
residuals. All points and finite-difference steps are recorded.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess

from prepare_nongrey_sources import digest

KB=1.380649e-16
MU=1.66053906660e-24
FIELDS=['F','U','S','P','cv','dP_dlnT','dP_dlnrho','electron_rs','Gamma','Tp_over_T']


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source_manifest','output']:p.add_argument(name,type=Path)
    a=p.parse_args();source=json.loads(a.source_manifest.read_text())
    if source['probe_kind']!='forced-phase interaction terms':
        raise ValueError('requires the explicit phase probe, not the native density inversion')
    exe=Path(source['executable'])
    if digest(exe)!=source['executable_sha256']:raise ValueError('source probe changed')
    # Both branches are evaluated as formulae. Their formal thermodynamic
    # identities must hold even though branch stability is not established.
    states=[(phase,z,mass,rho,t) for phase in [0,1]
            for z,mass in [(2.,4.),(6.,12.)]
            for rho in [1e3,1e5]
            for t in [1e5,1e6]]
    steps=[1e-3,3e-4,1e-4]
    queries=[]
    for state in states:
        queries.append(state)
        phase,z,mass,rho,t=state
        for h in steps:
            queries.extend([(phase,z,mass,rho,t*math.exp(h)),
                            (phase,z,mass,rho,t*math.exp(-h)),
                            (phase,z,mass,rho*math.exp(h),t),
                            (phase,z,mass,rho*math.exp(-h),t)])
    text=''.join(' '.join(format(x,'.17g') for x in q)+'\n' for q in queries)
    result=subprocess.run([str(exe)],input=text,text=True,capture_output=True,check=True)
    values=[dict(zip(FIELDS,map(float,line.split()),strict=True)) for line in result.stdout.splitlines()]
    if len(values)!=len(queries) or any(not math.isfinite(x) for r in values for x in r.values()):
        raise ValueError('invalid phase-probe output')
    rows=[];i=0
    for state in states:
        phase,z,mass,rho,t=state;central=values[i];i+=1
        rgas=KB/(mass*MU);pideal=rho*rgas*t;derivatives=[]
        for h in steps:
            tp,tm,rp,rm=values[i:i+4];i+=4
            fd={'S':-(tp['F']-tm['F'])/(2*h*t),
                'P':rho*(rp['F']-rm['F'])/(2*h),
                'cv':(tp['U']-tm['U'])/(2*h*t),
                'dP_dlnT':(tp['P']-tm['P'])/(2*h),
                'dP_dlnrho':(rp['P']-rm['P'])/(2*h)}
            norms={'S':rgas,'P':pideal,'cv':rgas,'dP_dlnT':pideal,'dP_dlnrho':pideal}
            derivatives.append({'logarithmic_step':h,'finite_difference':fd,
                                'defect_in_ideal_ion_units':{k:(central[k]-v)/norms[k] for k,v in fd.items()}})
        coarse,fine=derivatives[-2:]
        ratio=(coarse['logarithmic_step']/fine['logarithmic_step'])**2
        extrapolated={k:(ratio*fine['finite_difference'][k]-coarse['finite_difference'][k])/(ratio-1)
                      for k in fine['finite_difference']}
        rows.append({'phase':phase,'Z':z,'A':mass,'rho':rho,'T':t,'source':central,
                     'F_minus_U_plus_TS_in_ideal_ion_units':(central['F']-central['U']+t*central['S'])/(rgas*t),
                     'richardson_finite_difference':extrapolated,
                     'richardson_defect_in_ideal_ion_units':{k:(central[k]-v)/norms[k] for k,v in extrapolated.items()},
                     'differences':derivatives})
    if digest(exe)!=source['executable_sha256']:raise ValueError('probe changed during audit')
    maxima={k:max(abs(row['differences'][-1]['defect_in_ideal_ion_units'][k]) for row in rows)
            for k in ['S','P','cv','dP_dlnT','dP_dlnrho']}
    report={'scope':__doc__,'source_manifest':source,'source_manifest_sha256':digest(a.source_manifest),
            'script_sha256':digest(__file__),'states':rows,'finest_step_maximum_defects':maxima,
            'richardson_maximum_defects':{k:max(abs(row['richardson_defect_in_ideal_ion_units'][k]) for row in rows)
                                         for k in maxima}}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'finest_step':maxima,'richardson':report['richardson_maximum_defects']}))


if __name__=='__main__':main()
