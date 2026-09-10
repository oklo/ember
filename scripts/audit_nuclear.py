#!/usr/bin/env python3
"""Audit the nuclear prescriptions on a saved ember-evolve final profile."""
import argparse
import json
import subprocess
from pathlib import Path

COLUMNS = ['T_K','rho_g_cm3','bare_pp_sfii_over_legacy','bare_33_sfii_over_legacy','bare_34_sfii_over_legacy',
           'legacy_pp_H','legacy_33_H','sfii_legacy_eps','debye_pp_H','debye_33_H','sfii_debye_eps',
           'svh_pp_H','svh_33_H','sfii_svh_eps','electron_eta','electron_susceptibility','gamma_e','pp_zeta',
           'eps_nuclear_neutrinos','reduced_ppII_over_pp']

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('probe',type=Path);p.add_argument('profile',type=Path);p.add_argument('output',type=Path)
    args=p.parse_args()
    data=json.loads(args.profile.read_text());profile=data['profile']
    weights=[profile[0][0]]+[0.]*(len(profile)-1)
    for i in range(1,len(profile)):
        dm=.5*(profile[i][0]-profile[i-1][0]);weights[i-1]+=dm;weights[i]+=dm
    # The nuclear model returns exactly zero below this declared cutoff.
    active=[(w,row) for w,row in zip(weights,profile) if row[3]>=1e5]
    rows=''.join(f'{r[3]:.17g} {r[2]:.17g} {r[5]:.17g} {r[6]:.17g}\n' for _,r in active)
    metal=data.get('metal_inventory','').lower().startswith('gs98')
    result=subprocess.run([str(args.probe.resolve())]+(['--gs98'] if metal else []),input=rows,text=True,capture_output=True,check=True)
    values=[dict(zip(COLUMNS,map(float,line.split()),strict=True)) for line in result.stdout.splitlines()]
    assert len(values)==len(active)
    powers={key:sum(w*v[key] for (w,_),v in zip(active,values)) for key in ['sfii_legacy_eps','sfii_debye_eps','sfii_svh_eps','eps_nuclear_neutrinos']}
    audit={
        'profile':str(args.profile), 'nuclear_model':data['nuclear_model'], 'points':data['points'],
        'age_yr':data['history'][-1][0], 'central':values[0], 'integrated_power_erg_s':powers,
        'max_pp_zeta_where_T_ge_1e5':max(v['pp_zeta'] for v in values),
        'max_reduced_ppII_over_pp':max(v['reduced_ppII_over_pp'] for v in values),
        'metal_inventory':data.get('metal_inventory','carried_isotopes'),
        'interpretation':'Prescription comparisons at identical final thermal/composition states; not independently relaxed stars or calibrated uncertainty bounds. Electron gas fully ionized; Z=.02 with the selected metal inventory. ppII is a reduced steady branch.'}
    args.output.write_text(json.dumps(audit,indent=2)+'\n')

if __name__=='__main__':
    main()
