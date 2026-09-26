#!/usr/bin/env python3
"""Compare source grey means, frequency groups and independent full spectra.

Distinguish frequency averaging differences from the refractive correction;
do not change completed reports or their criteria.
"""
import argparse
import json
import math
from pathlib import Path

from scipy.integrate import quad

from audit_tops_spectral_means import digest, source, Integrals
from audit_tops_electron_dispersion import mean, electron_moments, KEV, KB, RW


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('spectral_plan', 'reduction_report', 'output'):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('keep completed comparisons')
    plan = json.loads(a.spectral_plan.read_text())
    old = json.loads(a.reduction_report.read_text())
    inputs = dict(old['input_sha256'])
    for path in [a.spectral_plan, a.reduction_report, Path(old['table']), Path(__file__)]:
        inputs[str(path.resolve())] = digest(path)
    if digest(Path(old['table'])) != old['table_sha256']:
        raise ValueError('candidate table changed')
    for name,h in inputs.items():
        if digest(Path(name)) != h:
            raise ValueError('changed input: '+name)
    refs = {(r['temperature_keV'],r['density_atomic_g_cm3']):r for r in old['records']}
    rows = []
    for job in plan['requests']:
        root = Path(job['work'])
        spectra, grey = source(root)
        receipt = json.loads((root/'receipt.json').read_text())
        if (receipt['X'],receipt['Z']) != (old['X'],old['Z']):
            raise ValueError('spectral composition differs')
        for name in ('source.txt','request.json','receipt.json','recipe.json'):
            path = root/name; inputs[str(path.resolve())] = digest(path)
        for key,data in spectra.items():
            row = refs[key]
            e = electron_moments(key[0]*KEV/KB,row['free_electron_density_cm3'])
            minimum_u = float(data[0,0]/key[0])
            cutoff_covered = row['cutoff_u'] >= minimum_u
            calculations=[]
            for order in (8,16):
                integrals = Integrals(data,key[0],order)
                uncut = float(RW/integrals.tail[0,1])
                corrected = (float(mean(data,key[0],row['cutoff_u'],e['vstar_squared'],order))
                             if cutoff_covered else None)
                calculations.append({'order':order,'uncut':uncut,'refractive':corrected})
            low,high = calculations
            error=max(abs(high[k]/low[k]-1) for k in ('uncut','refractive') if high[k] is not None)
            if error>1e-7:
                raise ValueError('full-spectrum quadrature did not converge')
            rows.append({'temperature_keV':key[0],'density_atomic_g_cm3':key[1],
                'group_mean':row['uncut_recombined_rosseland'],
                'source_uncut_grey_mean':row['uncut_source_rosseland'],
                'full_spectrum_uncut_mean':high['uncut'],
                'full_spectrum_refractive_mean':high['refractive'],
                'group_vs_full_relative_error':row['uncut_recombined_rosseland']/high['uncut']-1,
                'source_grey_vs_full_relative_error':row['uncut_source_rosseland']/high['uncut']-1,
                'full_spectrum_refractive_relative_change':(high['refractive']/high['uncut']-1
                                                           if cutoff_covered else None),
                'group_refractive_relative_change':row['rosseland_atomic_cm2_g']/row['uncut_recombined_rosseland']-1,
                'source_minimum_u':minimum_u, 'source_covers_plasma_cutoff':cutoff_covered,
                'thermal_weight_below_source':quad(lambda u:u**4*math.exp(-u)/(-math.expm1(-u))**2,
                                                   0,minimum_u)[0]/RW,
                'quadrature_relative_change':error,'full_spectrum_source':str(root)})
    for name,h in inputs.items():
        if digest(Path(name)) != h:
            raise ValueError('input changed during comparison')
    maximum=max(abs(r['group_vs_full_relative_error']) for r in rows)
    result={'scope':__doc__,'accepted_for_stellar_opacity':False,
        'X':old['X'],'Z':old['Z'],'states':len(rows),
        'maximum_group_vs_full_relative_error':maximum,
        'group_representation_relative_criterion':.005,
        'group_representation_check_passed':maximum<=.005,
        'original_uncut_recovery_check_passed':old['uncut_recovery_check_passed'],
        'comparison_quantity':'Uncut mean integrated over the available monochromatic source. No missing low-frequency opacity is supplied; refractive means are omitted when the plasma cutoff is outside the source frequency range.',
        'records':rows,'input_sha256':inputs}
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='input_sha256'}),flush=True)


if __name__ == '__main__':
    main()
