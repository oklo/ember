#!/usr/bin/env python3
"""Report bracketed exhaustion and cooling milestones without extrapolation.

All recrossings are retained. A threshold crossing is not a claim of permanent
extinction, a white-dwarf classification, or a converged endpoint age.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def crossings(rows, quantity, threshold):
    result=[]
    for lo,hi in zip(rows,rows[1:]):
        a,b=quantity(lo),quantity(hi)
        if (a<=threshold)==(b<=threshold):continue
        result.append({'direction':'down' if b<=threshold else 'up',
                       'age_bracket_yr':[lo['age_yr'],hi['age_yr']],
                       'value_bracket':[a,b],
                       'hydrogen_mass_bracket_Msun':[lo.get('hydrogen_mass_Msun'),hi.get('hydrogen_mass_Msun')],
                       'central_X_bracket':[lo['central_X'],hi['central_X']]})
    return {'threshold':threshold,'already_below_at_segment_start':quantity(rows[0])<=threshold,
            'crossings':result,'below_at_segment_end':quantity(rows[-1])<=threshold}


def summarize(track):
    columns=track['columns']
    if len(columns)!=len(set(columns)) or not track['history']:
        raise ValueError('duplicate columns or empty history')
    rows=[dict(zip(columns,row,strict=True)) for row in track['history']]
    required=['age_yr','central_X','Teff_K','L_Lsun']
    for row in rows:
        if any(k not in row or not isinstance(row[k],(int,float)) or not math.isfinite(row[k]) for k in required):
            raise ValueError('missing or nonfinite endpoint diagnostic')
        if row['age_yr']<0 or not 0<=row['central_X']<=1 or row['Teff_K']<=0 or row['L_Lsun']<=0:
            raise ValueError('invalid physical endpoint diagnostic')
    if any(b['age_yr']<=a['age_yr'] for a,b in zip(rows,rows[1:])):
        raise ValueError('history ages must be strictly increasing')
    nuclear=all('nuclear_deposited_Lsun' in r and isinstance(r['nuclear_deposited_Lsun'],(int,float))
                and math.isfinite(r['nuclear_deposited_Lsun']) and r['nuclear_deposited_Lsun']>=0 for r in rows)
    return {'scope':__doc__,'age_range_yr':[rows[0]['age_yr'],rows[-1]['age_yr']],
            'accepted_states':len(rows),'run_reached_requested_age':track.get('converged',False),
            'core_hydrogen':[crossings(rows,lambda r:r['central_X'],v) for v in [1e-3,1e-4,1e-5]],
            'nuclear_deposited_over_photon_luminosity':
                [crossings(rows,lambda r:r['nuclear_deposited_Lsun']/r['L_Lsun'],v) for v in [.1,.01,.001]] if nuclear else None,
            'nuclear_diagnostic_scope':'Deposited luminosity of the selected network; pp-only in the present driver. Missing legacy data are not interpreted as zero burning.',
            'temperature':[crossings(rows,lambda r:r['Teff_K'],v) for v in [1000,500,100]],
            'temperature_scope':'Downward crossings are cooling-temperature milestones only. Establish a degenerate remnant, residual-burning regime and source validity separately.',
            'final':{k:rows[-1].get(k) for k in required+['hydrogen_mass_Msun','nuclear_deposited_Lsun']}}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('track',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();raw=a.track.read_bytes();result=summarize(json.loads(raw))
    result['input_sha256']=hashlib.sha256(raw).hexdigest()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result['final']))


if __name__=='__main__':main()
