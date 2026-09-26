#!/usr/bin/env python3
"""Fuel-budget and homology comparisons; no new stellar evolution or causal fit."""
import csv
import hashlib
import json
import math
from pathlib import Path
from datetime import datetime,timezone


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open() as f:
        return [{k:v if k=='run' else float(v) if v else None for k,v in r.items()}
                for r in csv.DictReader(f)]


def crossing(history, target):
    for a,b in zip(history,history[1:]):
        if a['central_X'] > target >= b['central_X']:
            f=(a['central_X']-target)/(a['central_X']-b['central_X'])
            return dict(central_X=target,age_bracket_yr=[a['age_yr'],b['age_yr']],
                        interpolated_age_yr=a['age_yr']+f*(b['age_yr']-a['age_yr']))
    return None


def main():
    directory=Path('docs/reports/2026-09-11')
    ep=directory/'evolution_latest_provenance.json';fp=directory/'f77_history_provenance.json'
    eprov,fprov=json.loads(ep.read_text()),json.loads(fp.read_text())
    eh=directory/eprov['history_csv'];fh=directory/'f77_history.csv'
    assert sha(eh)==eprov['history_csv_sha256'] and sha(fh)==fprov['csv_sha256']
    published_path=Path('docs/results/f77_lba97_comparison_v1.json')
    published=json.loads(published_path.read_text())['published']
    lba_pdf=Path('/Users/greglaughlin/Projects/low_mass_stars/papers/lba97.pdf')
    assert sha(lba_pdf)==published['source_pdf_sha256']
    cb_pdf=Path('/tmp/ember-lba97-lifetime-analytic-v1/chabrier_baraffe_1997.pdf')
    assert sha(cb_pdf)=='e393536343ac3ef02be7154f95de974fa0dfd268f67addc57e87ee8ef0562fa1'
    ember,f77=rows(eh),rows(fh)
    sigma,lsun,rsun=5.670374419e-5,3.828e33,6.957e10
    G,msun,year=6.67430e-8,1.98847e33,31557600.
    initial=published['main_sequence_start']
    T0,L0=initial['Teff_K'],10**initial['log10_L_Lsun']
    R0=math.sqrt(L0*lsun/(4*math.pi*sigma*T0**4))/rsun
    models=[dict(name='LBA97',age_yr=initial['age_yr'],Teff_K=T0,L_Lsun=L0,
                 R_Rsun=R0,radius_inferred_from_L_and_T=True,
                 peak_Y3_age_yr=published['peak_helium3']['age_yr'],
                 peak_Y3=published['peak_helium3']['central_Y3'])]
    e0=ember[0];peak=max(ember,key=lambda r:r['central_Y3'])
    models.append(dict(name='Ember',age_yr=e0['age_yr'],Teff_K=e0['Teff_K'],
                       L_Lsun=e0['L_Lsun'],R_Rsun=e0['R_Rsun'],
                       peak_Y3_age_yr=peak['age_yr'],peak_Y3=peak['central_Y3']))
    for name in ['ajr-20000','ferguson-20000']:
        sequence=[r for r in f77 if r['run']==name]
        start=min((r for r in sequence if r['age_yr']<1e10),key=lambda r:r['L_Lsun'])
        peak=max(sequence,key=lambda r:r['central_Y3'])
        models.append(dict(name='F77 '+name,age_yr=start['age_yr'],Teff_K=start['Teff_K'],
                           L_Lsun=start['L_Lsun'],R_Rsun=start['R_Rsun'],
                           starting_stage='early luminosity minimum; not exact nuclear/photospheric equality',
                           peak_Y3_age_yr=peak['age_yr'],peak_Y3=peak['central_Y3']))
    models.append(dict(name='Chabrier & Baraffe 1997 Table 2',age_yr=1e9,
                       Teff_K=2811.,L_Lsun=10**(-3.070),R_Rsun=.863e10/rsun,
                       central_T_K=10**6.658,central_density_g_cm3=10**2.607,
                       metallicity_M_H=0.,helium_mass_fraction=.275))
    for model in models:
        t=model['Teff_K']/T0
        model['luminosity_over_LBA97_start']=model['L_Lsun']/L0
        model['radius_over_LBA97_start']=model['R_Rsun']/R0
        model['temperature_over_LBA97_start']=t
        model['kelvin_helmholtz_GM2_over_RL_yr']=G*(.1*msun)**2/(model['R_Rsun']*rsun*model['L_Lsun']*lsun)/year
        model['homology_luminosity_ratio_for_pp_exponents_4_and_6']=[t**(4*(nu+3)/(nu+5)) for nu in [4.,6.]]
        if 'peak_Y3_age_yr' in model:
            model['LBA97_over_model_helium3_peak_age']=models[0]['peak_Y3_age_yr']/model['peak_Y3_age_yr']
    sequences={'Ember':ember,**{n:[r for r in f77 if r['run']==n] for n in ['ajr-20000','ferguson-20000']}}
    milestones={n:[crossing(s,x) for x in [.2034,.001]] for n,s in sequences.items()}
    photon=sum(.5*(a['L_Lsun']+b['L_Lsun'])*(b['age_yr']-a['age_yr'])*year*lsun for a,b in zip(ember,ember[1:]))
    q=6e18
    paths=[Path(__file__),ep,fp,eh,fh,published_path,lba_pdf,cb_pdf]
    report=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),models=models,
        milestones=milestones,homology=dict(assumptions='Same mass and composition; homologous nearly ideal-gas interior; epsilon_pp proportional to rho T^nu with nu between 4 and 6; L_nuclear equals L_surface.',
            radius_scaling='R proportional to Teff^[-4/(nu+5)]',
            luminosity_scaling='L proportional to Teff^[4(nu+3)/(nu+5)]'),
        fuel_budget=dict(assumed_deposited_energy_per_gram_H_erg_g=q,
            initial_H_mass_Msun=.07,initial_fuel_erg=q*.07*msun,
            initial_gravitational_scale_erg=G*(.1*msun)**2/(R0*rsun),
            ember_photon_energy_integral_erg=photon,
            ember_mean_photon_luminosity_Lsun=photon/(ember[-1]['age_yr']*year*lsun),
            ember_H_consumed_Msun=ember[0]['hydrogen_mass_Msun']-ember[-1]['hydrogen_mass_Msun']),
        limitations=['A cooler initial guess with unchanged physics relaxes on a thermal timescale; it does not set the trillion-year lifetime.',
            'Initial luminosity ratios need not persist with composition. This is not a controlled attribution to opacity, atmosphere or EOS.',
            'Helium-3 peak abundances differ, particularly in F77, so their ages are only a rough clock comparison.',
            'Core hydrogen thresholds are not equal total-fuel states once the star is incompletely mixed.',
            'The energy-per-hydrogen estimate neglects transient He3 storage and detailed neutrino losses.',
            'No full Ember nuclear lifetime or WD cooling age has been computed.'],
        input_sha256={str(p):sha(p) for p in paths})
    with Path('docs/results/lba97_lifetime_scale_v1.json').open('x') as f:
        json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(models=models,fuel_budget=report['fuel_budget'],milestones=milestones)))


if __name__=='__main__':main()
