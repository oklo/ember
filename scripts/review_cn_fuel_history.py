#!/usr/bin/env python3
"""Review Fable's sealed exposure estimate without changing its stellar history."""
import argparse
import hashlib
import json
from pathlib import Path
import textwrap
from datetime import datetime, timezone

import numpy as np


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('output',type=Path)
    args=ap.parse_args();assert not args.output.exists()
    root=Path(__file__).resolve().parents[1]
    base=root/'docs/research/fable';batch=base/'results/fable-cn-fuel-v1'
    ready=json.loads((batch/'READY.json').read_text());identities={}
    for group in ('files_sha256','input_sha256'):
        for name,h in ready[group].items():
            p=((batch if group=='files_sha256' else base)/name).resolve()
            assert hashlib.sha256(p.read_bytes()).hexdigest()==h,p
            identities[str(p)]=h
    for p in (batch/'READY.json',Path(__file__)):
        identities[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    # Reuse only the frozen input/rate preparation. The estimates and their
    # comparisons below are separate from the worker's result construction.
    source=(batch/'cn_fuel_history.py').read_text()
    header,body=source.split('\ndef main():',1)
    header=header.replace('open("cno_rate_v2.py")',f'open({str(batch/"cno_rate_v2.py")!r})')
    setup=textwrap.dedent(body.split('    results = {}',1)[0])
    setup=setup.replace('root = "../../../../"',f'root = {str(root)+"/"!r}')
    ns={};exec(header,ns);exec(setup,ns)
    cum_c,cum_n,dm,rc,rn=(ns[k] for k in ('cum_c','cum_n','dm','ratio_c','ratio_n'))
    iend=ns['i_end'];yc,yn=ns['Y_C'],ns['Y_N'];cases=[]
    for label,fc,fn in [('mass_weighted',ns['homog_c'],ns['homog_n']),('central',1.,1.),('none',0.,0.)]:
        ec=np.minimum(cum_c,cum_c[iend])[:,None]*fc+np.maximum(cum_c-cum_c[iend],0.)[:,None]*rc
        en=np.minimum(cum_n,cum_n[iend])[:,None]*fn+np.maximum(cum_n-cum_n[iend],0.)[:,None]*rn
        dc=np.diff(ec,axis=0);dn=np.diff(en,axis=0)
        phi=np.ones_like(dc);np.divide(-np.expm1(-dc),dc,out=phi,where=dc>0)
        cycles=((yn+yc)*dn-yc*np.exp(-ec[:-1])*dn*phi).sum(axis=0)
        final_shortcut=(yn+yc*(-np.expm1(-ec[-1])))*en[-1]
        assert np.min(cycles)>=0 and np.all(cycles<=final_shortcut+1e-20)
        conversion=2*yc*(-np.expm1(-ec[-1]))
        cases.append(dict(case=label,
            original_cycle_hydrogen_Msun=float(dm@(4*final_shortcut)/ns['MSUN']),
            time_integrated_cycle_hydrogen_Msun=float(dm@(4*cycles)/ns['MSUN']),
            shortcut_cycle_relative_overestimate=float(dm@(final_shortcut-cycles)/(dm@cycles)),
            time_integrated_total_hydrogen_Msun=float(dm@(conversion+4*cycles)/ns['MSUN']),
            central_dX=float(conversion[0]+4*cycles[0])))
    exponent_before=ns['homog_c']*cum_c[iend]
    correct_after=exponent_before+cum_c[iend+1]-cum_c[iend]
    shortcut_after=cum_c[iend+1]
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='needs_changes',
        batch='fable-cn-fuel-v1',input_sha256=identities,cases=cases,
        proton_number_correction_at_fixed_screening=ns['A_H'],
        convective_transition=dict(age_before_years=float(ns['age'][iend]),
            age_after_years=float(ns['age'][iend+1]),
            carbon_fraction_before=float(-np.expm1(-exponent_before)),
            continuous_carbon_fraction_after=float(-np.expm1(-correct_after)),
            worker_carbon_fraction_after=float(-np.expm1(-shortcut_after))),
        accepted=['Historical carbon conversion and CN cycling require a consistent fuel treatment before predicting core exhaustion.',
                  'The worker preserved the pp-only stellar trajectory and did not subtract this estimate from its central hydrogen.'],
        required_changes=['Use the evolving N14 inventory in the cycle integral, rather than its final inventory throughout.',
            'Carry the mixed carbon exposure continuously across the end of full convection.',
            'Use one abundance convention throughout: the declared baryon-basis hydrogen needs proton molality X, not X/A_H.',
            'Keep the present-profile scaling and effective burning mass as conditional approximations, not error bounds.',
            'Choose a practical fuel-conserving restart approximation with quantified fixed-mixture EOS and opacity errors (P111).'],
        limitations=['This audit retains the worker\'s rate history, screening and spatial approximation.',
            'The revised cycle integral is exact for piecewise-constant capture frequencies with the same integrated exposures; it is not a stellar evolution calculation.',
            'The quoted proton correction holds screening fixed; the entire abundance/screening convention still needs a consistent comparison.'])
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='input_sha256'},indent=2))


if __name__=='__main__':
    main()
