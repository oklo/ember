#!/usr/bin/env python3
"""Offline FastChem equilibrium/rainout audit of validated TLUSTY profiles.

Requires numpy and the separately built GPL FastChem adapter. This diagnoses
condensation on fixed gas-atmosphere structures; it does not feed cloud or
depletion opacity back into radiative transfer.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
import math
import math
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from generate_nongrey_grid import composition
from import_nongrey_grid import import_grid, read_text

ROOT=Path(__file__).resolve().parents[1]
def sha(data):return hashlib.sha256(data).hexdigest()
def pack(value):return (json.dumps(value,separators=(',',':'),allow_nan=False)+'\n').encode()


def closure(result, profile, numbers, masses):
    if result.get('mode') not in ['equilibrium','rainout']:
        raise ValueError('unrecognized chemistry mode')
    elements=result['elements'];h=elements.index('H');electron=elements.index('e-')
    rows=result['rows'];n=len(rows)
    if result['status'] or n!=len(profile) or any(r['flag'] or not all(r['element_conserved']) for r in rows):
        raise ValueError('chemistry source failed convergence or conservation')
    if any(not math.isclose(r['T_K'],p[3],rel_tol=1e-12) or
           not math.isclose(r['P_bar']*1e6,p[6],rel_tol=1e-12)
           for r,p in zip(rows,profile,strict=True)):
        raise ValueError('chemistry state does not match atmosphere profile')
    gas=np.array([r['gas'] for r in rows]);cond=np.array([r['condensed'] for r in rows])
    if not np.isfinite(gas).all() or not np.isfinite(cond).all() or gas.min()<0 or cond.min()<0:
        raise ValueError('invalid species number density')
    atoms_gas=gas@np.array([s['stoichiometry'] for s in result['gas_species']])
    atoms_cond=cond@np.array([s['stoichiometry'] for s in result['condensates']])
    atoms=atoms_gas+atoms_cond
    total_density=np.array([r['total_element_density'] for r in rows])
    if (not np.isfinite(atoms).all() or np.any(atoms[:,h]<=0) or
            not np.isfinite(total_density).all() or np.any(total_density<=0)):
        raise ValueError('invalid conserved nuclei density')
    selected=np.arange(len(elements))!=electron
    initial=np.array([numbers.get(s,0) for s in elements]);initial/=initial[h]
    reservoir=np.tile(initial,(n,1))
    if result['mode']=='rainout':reservoir[1:]=atoms_gas[:-1]/atoms_gas[:-1,h,None]
    ratio=atoms/atoms[:,h,None]
    element_error=float(np.max(np.abs(ratio[:,selected]/reservoir[:,selected]-1)))
    pressure_error=float(np.max(np.abs(gas.sum(axis=1)*1.380649e-16*
        np.array([r['T_K'] for r in rows])/np.array([r['P_bar'] for r in rows])/1e6-1)))
    nuclei_error=float(np.max(np.abs(atoms[:,selected].sum(axis=1)/total_density-1)))
    if not all(math.isfinite(v) for v in [element_error,pressure_error,nuclei_error]) or max(element_error,pressure_error,nuclei_error)>2e-7:
        raise ValueError(f'independent chemistry closure failed: {element_error}, {pressure_error}, {nuclei_error}')
    weights=np.array([masses.get(s,0) for s in elements])
    condensed_mass=(atoms_cond@weights)/(atoms@weights)
    depletion=1-atoms_gas[:,selected]/atoms_gas[:,h,None]/initial[selected]
    tau=np.array([r[2] for r in profile])
    thresholds={}
    for threshold in [0,.001,.01,.1,1,100]:
        mask=tau>=threshold
        thresholds[str(threshold)]={'max_condensed_mass_fraction':float(condensed_mass[mask].max()) if mask.any() else None,
            'max_element_gas_depletion':float(depletion[mask].max()) if mask.any() else None}
    present=[]
    for j,s in enumerate(result['condensates']):
        fraction=cond[:,j]/np.array([r['total_element_density'] for r in rows])
        if fraction.max()>1e-14:
            active=fraction>1e-14
            present.append({'species':s['symbol'],'max_number_per_nucleus':float(fraction.max()),
                            'largest_tau':float(tau[active].max())})
    return {'element_closure':element_error,'pressure_closure':pressure_error,'nuclei_closure':nuclei_error,
            'thresholds':thresholds,'condensates':present}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('manifest',type=Path);p.add_argument('source_receipt',type=Path)
    p.add_argument('probe',type=Path);p.add_argument('work',type=Path)
    p.add_argument('--jobs',type=int,default=4);p.add_argument('--limit',type=int)
    a=p.parse_args()
    if not 1<=a.jobs<=8:raise ValueError('invalid jobs')
    spec=json.loads(a.manifest.read_text());source=json.loads(a.source_receipt.read_text())
    with tempfile.TemporaryDirectory() as temp:import_grid(a.manifest,Path(temp)/'verified.dat')
    metadata=json.loads((ROOT/'data/atmosphere/sources/synple-elements.json').read_text())
    a.work.mkdir(parents=True,exist_ok=True)
    provenance={**source,'probe_sha256':sha(a.probe.read_bytes()),
        'adapter_sha256':sha(Path(__file__).with_name('fastchem_probe.cpp').read_bytes()),
        'manifest_sha256':sha(a.manifest.read_bytes()),'accuracy':1e-8}
    for file in ['input/logK/logK.dat','input/logK/logK_condensates.dat']:
        provenance[file+'_sha256']=sha((Path(source['source'])/file).read_bytes())
    (a.work/'provenance.json').write_bytes(pack(provenance))

    def run(job):
        i,model,mode=job
        d=a.work/f'model-{i:03d}'/mode;d.mkdir(parents=True,exist_ok=True)
        logpath=a.manifest.parent/model['log']
        if sha(logpath.read_bytes())!=model['log_sha256']:raise ValueError('atmosphere checksum mismatch')
        log=read_text(logpath);profile=[]
        for line in log.rsplit('FINAL MODEL ATMOSPHERE',1)[1].splitlines():
            t=line.replace('D','E').split()
            if len(t)==11 and t[0].isdigit():profile.append(list(map(float,t)))
        profile.reverse()
        abundance,weights=composition(model['XH'],model['X3'],spec['metals'])
        numbers={};masses={}
        for s,n,w in zip(metadata['symbol'],abundance,weights,strict=True):
            if n>1e-90:numbers[s.capitalize()]=n;masses[s.capitalize()]=w*1.67333e-24/1.66053906660e-24
        abundance_text='# Shared GS98 baryonic element numbers\ne- 0\n'+''.join(
            f'{s} {12+math.log10(n):.17g}\n' for s,n in numbers.items())
        inputs=''.join(f'{r[3]:.17g} {r[6]/1e6:.17g}\n' for r in profile)
        (d/'abundances.dat').write_text(abundance_text);(d/'input.dat').write_text(inputs)
        fingerprint=sha(pack(provenance)+abundance_text.encode()+inputs.encode()+mode.encode())
        output=d/'source.json.gz';receipt=d/'receipt.json'
        if output.exists() and receipt.exists():
            saved=json.loads(receipt.read_text())
            if saved['input_sha256']!=fingerprint or saved['output_sha256']!=sha(output.read_bytes()):
                raise ValueError('cached chemistry checksum mismatch')
            result=json.loads(gzip.decompress(output.read_bytes()))
        else:
            r=subprocess.run([str(a.probe.resolve()),source['source'],str((d/'abundances.dat').resolve()),mode],
                input=inputs,text=True,capture_output=True,timeout=600)
            (d/'run.log').write_text(r.stderr)
            if r.returncode:raise ValueError(f'{d}: chemistry failed; see run.log')
            result=json.loads(r.stdout);output.write_bytes(gzip.compress(r.stdout.encode(),mtime=0))
            receipt.write_bytes(pack({'input_sha256':fingerprint,'output_sha256':sha(output.read_bytes())}))
        audit=closure(result,profile,numbers,masses)
        record={'coordinates':[model[k] for k in ['XH','X3','teff_K','log_g']],
                'mode':mode,**audit,'source':str(output.relative_to(a.work)),
                'source_sha256':sha(output.read_bytes())}
        (d/'audit.json').write_bytes(pack(record));print(i,mode,audit['thresholds'],flush=True)
        return record

    models=spec['models'][:a.limit] if a.limit else spec['models']
    jobs=[(i,m,mode) for i,m in enumerate(models) for mode in ['equilibrium','rainout']]
    with ThreadPoolExecutor(a.jobs) as pool:records=list(pool.map(run,jobs))
    (a.work/'audit.json').write_bytes(pack({'provenance':provenance,'models':len(models),'records':records,
        'limitation':'Fixed gas-atmosphere profiles; condensate chemistry and rainout audited without radiative feedback or grain opacity'}))

if __name__=='__main__':main()
