#!/usr/bin/env python3
"""Check cold frequency groups against full spectra at identical source states.

Reuse verified group quadrature, independently integrate full component spectra,
and compare both the refractive ratio and the mean opacity. These checks do
not establish interpolation accuracy or bound-electron dielectric physics.
"""
import argparse
import json
import math
from pathlib import Path

from audit_tops_spectral_means import source, Integrals
from audit_tops_electron_dispersion import electron_moments, mean, controls, KEV, KB, HBAR, ME, RW
from fetch_tops_composition import digest
from reduce_tops_group_factors import add_inputs, verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('plan',type=Path);p.add_argument('report',type=Path)
    a=p.parse_args()
    if a.report.exists():raise FileExistsError('preserve completed comparisons')
    plan=json.loads(a.plan.read_text());inputs={str(a.plan.resolve()):digest(a.plan)}
    for name in ('audit_tops_cool_full_spectra.py','audit_tops_spectral_means.py',
                 'audit_tops_electron_dispersion.py','fetch_tops_composition.py',
                 'reduce_tops_group_factors.py','import_tops_composition.py'):
        path=Path(__file__).with_name(name);inputs[str(path.resolve())]=digest(path)
    cached={}
    for j in plan['reference_groups']:
        path=Path(j['group_reduction_report'])
        if str(path) not in cached:
            r=json.loads(path.read_text());cached[str(path)]=r
            add_inputs(inputs,r['input_sha256']);inputs[str(path.resolve())]=digest(path)
    verify(inputs);limit_checks=controls();rows=[]
    for i,job in enumerate(plan['requests']):
        root=Path(job['work']);data,grey=source(root)
        receipt=json.loads((root/'receipt.json').read_text())
        requested={(t,rho) for t in job['temperatures_keV'] for rho in job['densities_g_cm3']}
        if set(data)!=requested:raise ValueError('full spectral coordinates differ from plan')
        for n in ('receipt.json','request.json','source.txt','recipe.json'):
            path=root/n;inputs[str(path.resolve())]=digest(path)
        part=next(v for v in plan['reference_groups'] if v['request_index']==i)
        reference=cached[part['group_reduction_report']]
        if (receipt['X'],receipt['Z'])!=(job['X'],job['Z']) or (reference['X'],reference['Z'])!=(job['X'],job['Z']):
            raise ValueError('source compositions differ')
        comp=(root/'source.txt').read_text().split('No. Fraction Mass Fraction  At. No.  Chem. Sym.  Mat ID.\n')[1].split('Temperature grid')[0]
        helium=[v.split() for v in comp.splitlines() if len(v.split())==5 and v.split()[3]=='He']
        if len(helium)!=1:raise ValueError('missing source helium normalization')
        number,mass=map(float,helium[0][:2]);atomic_mass=4.002602*number/mass
        for (t,rho),spectrum in data.items():
            group=next(r for r in reference['records'] if (r['temperature_keV'],r['density_atomic_g_cm3'])==(t,rho))
            group_receipt=json.loads((Path(group['source'])/'receipt.json').read_text())
            if any(receipt[k]!=group_receipt[k] for k in ('X','Z','metals')):
                raise ValueError('source metal mixtures differ')
            ne=rho/(atomic_mass*1.66053906660e-24)*grey[t,rho]['free_electrons_per_ion']
            if ne!=group['free_electron_density_cm3']:raise ValueError('source electron counts differ')
            if max(group['quadrature_relative_change'],group['uncut_quadrature_relative_change'])>1e-6:
                raise ValueError('group quadrature was not converged')
            e=electron_moments(t*KEV/KB,ne)
            cutoff=HBAR*math.sqrt(4*math.pi*4.80320471257e-10**2*ne/ME)/(t*KEV)*math.sqrt(e['plasma_frequency_squared_ratio'])
            trials=[];previous=None
            for order in (8,16,32,64):
                full=float(mean(spectrum,t,cutoff,e['vstar_squared'],order))
                uncut=float(RW/Integrals(spectrum,t,order).tail[0][1])
                if not all(math.isfinite(v) and v>0 for v in (full,uncut)):
                    raise ValueError('invalid full-spectrum integral')
                change=None if previous is None else max(abs(v/w-1) for v,w in zip((full,uncut),previous))
                trials.append({'order':order,'transport_mean':full,'uncut_mean':uncut,'relative_change':change})
                if change is not None and change<=1e-7:break
                previous=full,uncut
            else:raise ValueError('full-spectrum quadrature did not converge')
            rows.append({'X':job['X'],'Z':job['Z'],'temperature_keV':t,'density_atomic_g_cm3':rho,
                         'spectral_points':len(spectrum),'full_spectrum_source':str(root),'group_source':group['source'],
                         'free_electron_density_cm3':ne,'cutoff_u':cutoff,'full_quadrature_trials':trials,
                         'reference_rosseland_atomic_cm2_g':full,'reference_uncut_rosseland_atomic_cm2_g':uncut,
                         'reference_refractive_ratio':full/uncut,'group_refractive_ratio':group['refractive_ratio'],
                         'component_relative_error':group['rosseland_atomic_cm2_g']/full-1,
                         'refractive_ratio_relative_error':group['refractive_ratio']/(full/uncut)-1,
                         'native_grey_relative_error':group['uncut_source_rosseland']/uncut-1,
                         'normalized_group_relative_error':group['native_uncut_times_ratio_atomic_cm2_g']/full-1})
    verify(inputs)
    fields=('component_relative_error','refractive_ratio_relative_error','native_grey_relative_error','normalized_group_relative_error')
    summary={k:{'maximum_absolute_difference':max(abs(r[k]) for r in rows),
                'count_above_0p5_percent':sum(abs(r[k])>.005 for r in rows)} for k in fields}
    result={'scope':__doc__,'accepted_for_stellar_opacity':False,'states':len(rows),'limit_checks':limit_checks,
            'relative_comparison_criterion':.005,'summaries':summary,'records':rows,'input_sha256':inputs,
            'group_representation_check_passed':all(summary[k]['count_above_0p5_percent']==0 for k in fields)}
    a.report.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'states':len(rows),'summaries':summary,'group_representation_check_passed':result['group_representation_check_passed']}))
    if not result['group_representation_check_passed']:raise SystemExit(1)


if __name__=='__main__':main()
