#!/usr/bin/env python3
"""Repair only the under-resolved occupation reference in electron heat audit v1.

All original reports and sources remain immutable. The original quadrature
placed no subdivision just below the Fermi surface at eta=1e4 and missed its
narrow occupation deficit. Explicit subdivisions on both sides, spacing
changes and an independent Sommerfeld expansion resolve that reference.
Other verified comparisons are reused; their numerical criteria are unchanged.
"""
import argparse,hashlib,json,math
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import expit
from electron_ion_born import fermi_half
from electron_ion_heat import current_basis


def occupation_energy(eta,edge_width):
    upper=math.sqrt(max(eta,0)+edge_width)
    points=[math.sqrt(eta+offset) for offset in [-edge_width,0] if 0<eta+offset<upper**2]
    return quad(lambda t:2*t**4*expit(eta-t*t),0,upper,points=points,
                epsabs=1e-100,epsrel=2e-12,limit=200)[0]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('report',type=Path);a=p.parse_args()
    original=Path('docs/results/electron_ion_heat_v1.json');old=json.loads(original.read_text())
    assert old['outcome']=='failed'
    for path,h in old['inputs_sha256'].items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==h,path
    failed=[i for i,c in enumerate(old['checks']) if not c['passed']]
    assert len(failed)==1
    index=failed[0];bad=old['checks'][index]
    assert bad['name']=='independent ideal electron enthalpy' and bad['eta']==1e4
    checks=[];data=[]
    for eta in [100.,1e3,1e4]:
        values=[occupation_energy(eta,width) for width in [40.,60.,80.]]
        sommerfeld=.4*eta**2.5+math.pi**2/4*math.sqrt(eta)-7*math.pi**4/(960*eta**1.5)
        for label,value,reference,tolerance in [('occupation integration subdivisions',values[0],values[-1],2e-10),
                ('occupation integration tail',values[1],values[-1],2e-10),
                ('independent Sommerfeld occupation energy',values[-1],sommerfeld,2e-10)]:
            error=abs(value/reference-1)
            checks.append(dict(name=label,eta=eta,normalized_error=error,tolerance=tolerance,passed=bool(error<=tolerance)))
        basis=current_basis(eta);enthalpy=(5/3)*values[1]/fermi_half(eta)
        error=abs((eta+basis.centered_mean)/enthalpy-1)
        data.append(dict(eta=eta,energy_integrals=values,sommerfeld=sommerfeld,
                         enthalpy_from_current_basis=eta+basis.centered_mean,enthalpy_from_occupation=enthalpy,
                         relative_enthalpy_difference=error))
        if eta==bad['eta']:
            corrected={**bad,'normalized_error':error,'passed':bool(error<=bad['tolerance'])}
    result={**old,'created_utc':datetime.now(timezone.utc).isoformat(),
            'scope':__doc__+'\nThe retained collision calculations have the scope recorded in v1.',
            'checks':[corrected if i==index else c for i,c in enumerate(old['checks'])]+checks,
            'reference_correction':dict(previous_report=str(original),previous_report_sha256=hashlib.sha256(original.read_bytes()).hexdigest(),
                                        original_failed_check=bad,corrected_check=corrected,controls=data,
                                        repeated_collision_calculations=0,changed_criteria=False),
            'inputs_sha256':{**old['inputs_sha256'],str(original):hashlib.sha256(original.read_bytes()).hexdigest(),
                             str(Path(__file__)):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}}
    passed=all(c['passed'] for c in result['checks']);result['outcome']='passed' if passed else 'failed'
    with a.report.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(outcome=result['outcome'],checks=len(result['checks']),corrected_check=corrected,
                         repeated_collision_calculations=0,failures=[c for c in result['checks'] if not c['passed']])),flush=True)
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
