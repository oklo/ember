#!/usr/bin/env python3
"""Diagnostic latent-enthalpy scale from pinned FastChem equilibrium constants.

This is not a replacement EOS or a grain-enthalpy acceptance rule. Formation
enthalpies follow the van't Hoff derivative of log(Kp). The optional ideal
neutral-atom reference omits atomic excitation; the native gas EOS and its
partition functions remain a separate approximation. Compare multiple finite
difference intervals across phase boundaries before interpreting derivatives.
"""
import argparse
import bisect
import gzip
import json
from pathlib import Path
import subprocess

from generate_nongrey_grid import composition, SOURCE_HMASS
from prepare_nongrey_sources import SOURCES, digest

KB=1.380649e-16


def coefficients(path):
    lines=Path(path).read_text().splitlines();i=3;result={}
    while i<len(lines):
        if not lines[i].strip():i+=1;continue
        header=lines[i];symbol=header.split()[0]
        tokens=header.split(':',1)[1].split('#',1)[0].split()
        stoichiometry={tokens[j]:int(tokens[j+1]) for j in range(0,len(tokens),2)}
        limits=list(map(float,lines[i+2].split()))
        fits=[list(map(float,row.split())) for row in lines[i+3:i+3+len(limits)]]
        if not limits or any(len(row)!=5 for row in fits):raise ValueError('unrecognized thermochemical fit')
        result[symbol]={'limits':limits,'fits':fits,'stoichiometry':stoichiometry,'source':header}
        i+=3+len(limits)
    return result


def formation_enthalpy(t,record):
    if t>record['limits'][-1]:raise ValueError('enthalpy query exceeds source temperature support')
    a=record['fits'][bisect.bisect_left(record['limits'],t)]
    # Kp is dimensionless with partial pressures relative to 1 bar.
    return KB*(-a[0]+a[1]*t+a[3]*t*t+2*a[4]*t**3)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('model',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--temperature',type=float,required=True);p.add_argument('--pressure-bar',type=float,required=True)
    p.add_argument('--gas-cp',type=float,help='native gas heat capacity at exactly this T/P, for a diagnostic ratio')
    a=p.parse_args();provenance=json.loads((a.model/'provenance.json').read_text())
    fc=provenance['prepared']['depletion']['fastchem'];source=Path(fc['source'])
    data=source/'input/logK/logK_condensates.dat';fits=coefficients(data)
    if digest(fc['executable'])!=fc['executable_sha256']:raise ValueError('chemistry probe changed')
    abundance,mass=composition(provenance['XH'],provenance['X3'],provenance['specification']['metals'])
    symbols=json.loads((SOURCES/'synple-elements.json').read_text())['symbol']
    mass={s.capitalize():v*SOURCE_HMASS for s,v in zip(symbols,mass,strict=True)}
    steps=[.005,.002,.001,.0005,.0002]
    ts=[a.temperature]+[a.temperature*(1+sign*h) for h in steps for sign in [-1,1]]
    inputs=''.join(f'{t:.17g} {a.pressure_bar:.17g}\n' for t in ts)
    result=subprocess.run([fc['executable'],str(source),str((a.model/'condensate-abundances.dat').resolve()),'equilibrium'],
                          input=inputs,text=True,capture_output=True,check=True)
    chemistry=json.loads(result.stdout)
    if chemistry['status'] or len(chemistry['rows'])!=len(ts):raise ValueError('incomplete equilibrium scan')
    mean_mass=sum(n*mass[s] for s,n in zip(chemistry['elements'],chemistry['abundances'],strict=True) if s!='e-')/sum(chemistry['abundances'])
    values=[]
    for row in chemistry['rows']:
        if row['flag'] or not all(row['element_conserved']):raise ValueError('equilibrium conservation failure')
        formation=0.;reference=0.;condensed_mass=0.;species=[]
        for entry,n in zip(chemistry['condensates'],row['condensed'],strict=True):
            if n/row['total_element_density']<=1e-14:continue
            record=fits[entry['symbol']];n_per_mass=n/(row['total_element_density']*mean_mass)
            h=formation_enthalpy(row['T_K'],record)
            ideal=2.5*KB*row['T_K']*sum(record['stoichiometry'].values())
            fraction=n_per_mass*sum(mass[s]*v for s,v in record['stoichiometry'].items())
            formation+=n_per_mass*h;reference+=n_per_mass*(h+ideal);condensed_mass+=fraction
            species.append({'symbol':entry['symbol'],'mass_fraction':fraction,'formation_enthalpy_erg_per_particle':h})
        values.append({'T_K':row['T_K'],'formation_enthalpy_erg_g':formation,
                       'ideal_atom_referenced_enthalpy_erg_g':reference,'condensed_mass_fraction':condensed_mass,'species':species})
    derivatives=[]
    for i,h in enumerate(steps):
        lo,hi=values[1+2*i:3+2*i];dt=hi['T_K']-lo['T_K']
        d={k:(hi[k]-lo[k])/dt for k in ['formation_enthalpy_erg_g','ideal_atom_referenced_enthalpy_erg_g','condensed_mass_fraction']}
        if a.gas_cp:d['ideal_atom_referenced_cp_over_native_gas_cp']=d['ideal_atom_referenced_enthalpy_erg_g']/a.gas_cp
        derivatives.append({'relative_temperature_step':h,'derivatives_per_K':d})
    report={'description':__doc__,'model':str(a.model),'temperature_K':a.temperature,'pressure_bar':a.pressure_bar,
            'native_gas_cp_erg_g_K':a.gas_cp,'values':values,'finite_differences':derivatives,
            'provenance':{'probe_sha256':digest(fc['executable']),'thermochemical_data_sha256':digest(data),
                          'model_provenance_sha256':digest(a.model/'provenance.json')},
            'limitations':'Diagnostic formation-energy estimate. Not a complete gas-plus-grain EOS, pseudoadiabat, flux bound, or permission to relax the grain-enthalpy guard.'}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    a.output.with_suffix('.source.json.gz').write_bytes(gzip.compress(result.stdout.encode(),mtime=0))
    a.output.with_suffix('.input.dat').write_text(inputs)
    a.output.with_suffix('.log').write_text(result.stderr)
    print(json.dumps({'central':values[0],'derivatives':derivatives}))


if __name__=='__main__':main()
