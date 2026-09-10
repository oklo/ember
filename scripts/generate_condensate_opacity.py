#!/usr/bin/env python3
"""Generate the depleted-gas opacity endmember with original source receipts."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import shutil
from generate_nongrey_grid import opacity_inputs,execute,sequence,temperatures,composition,input_fingerprint
from nongrey_opacity import validate_table,merge_isotherms
from prepare_nongrey_sources import digest,SOURCES


def depletion_inputs(directory,prepared,spec,x,y,mode):
    fc=prepared['depletion']['fastchem'];abundance,_=composition(x,y,spec['metals'])
    symbols=json.loads((SOURCES/'synple-elements.json').read_text())['symbol']
    path=directory/'condensate-abundances.dat'
    path.write_text('# Same baryonic element numbers as the atmosphere\ne- 0\n'+''.join(
        f'{s.capitalize()} {12+math.log10(v):.17g}\n' for s,v in zip(symbols,abundance,strict=True) if v>1e-90))
    (directory/'ember-condensates.cfg').write_text(mode+'\n'+fc['source']+'\n'+str(path.resolve())+'\n')
    hashes={}
    for name,expected in fc['chemistry_data_sha256'].items():
        actual=digest(Path(fc['source'])/'input/logK'/name)
        if actual!=expected:raise ValueError('chemistry source identity changed')
        hashes[name]=actual
    (directory/'condensates.sha256').write_text(json.dumps(hashes,sort_keys=True)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('prepared',type=Path);p.add_argument('specification',type=Path);p.add_argument('work',type=Path)
    p.add_argument('--hydrogen',type=float,required=True);p.add_argument('--helium3',type=float,default=0)
    p.add_argument('--mode',choices=['gas','equilibrium'],default='equilibrium');p.add_argument('--jobs',type=int,default=4)
    p.add_argument('--indices',type=int,nargs='+')
    p.add_argument('--reuse-opacity',type=Path,
                   help='reuse independently revalidated identical isotherms from another material grid')
    a=p.parse_args();prepared=json.loads(a.prepared.read_text());spec=json.loads(a.specification.read_text())
    if not 1<=a.jobs<=8:raise ValueError('invalid worker count')
    if digest(prepared['synspec'])!=prepared['executables']['synspec']:raise ValueError('source executable changed')
    a.work.mkdir(parents=True,exist_ok=True);ts=temperatures(spec);rs=sequence(spec['log_density'])
    indices=a.indices if a.indices is not None else list(range(len(ts)))
    if len(set(indices))!=len(indices) or any(i<0 or i>=len(ts) for i in indices):raise ValueError('invalid isotherms')
    provenance={'prepared':prepared,'specification':spec,'XH':a.hydrogen,'X3':a.helium3,'mode':a.mode}
    receipt=a.work/'provenance.json'
    if receipt.exists() and json.loads(receipt.read_text())!=provenance:raise ValueError('work inputs changed')
    receipt.write_text(json.dumps(provenance,indent=2)+'\n')
    reuse={}
    if a.reuse_opacity:
        old=json.loads((a.reuse_opacity/'provenance.json').read_text())
        if any(old.get(k)!=provenance[k] for k in ['prepared','XH','X3','mode']):
            raise ValueError('reused opacity source or composition differs')
        for key in ['metals','log_density','opacity_frequencies','opacity_method','wavelength_A',
                    'line_threshold','synthesis_spacing_A','microturbulence_km_s']:
            if old['specification'].get(key)!=spec.get(key):raise ValueError('reused opacity physical settings differ')
        reuse={t:a.reuse_opacity/f'temperature-{i:03d}' for i,t in enumerate(temperatures(old['specification']))}
    def run(i):
        d=a.work/f'temperature-{i:03d}'
        if ts[i] in reuse:
            old= reuse[ts[i]];saved=json.loads((old/'completed.json').read_text())
            if input_fingerprint(prepared['executables']['synspec'],old)!=saved['input_sha256']:
                raise ValueError('reused opacity input fingerprint changed')
            if not {'fort.63','fort.29','run.log'}.issubset(saved['outputs']):
                raise ValueError('incomplete reused opacity source receipt')
            for name,expected in saved['outputs'].items():
                if digest(old/name)!=expected:raise ValueError('reused opacity output changed')
            abundance,_=composition(a.hydrogen,a.helium3,spec['metals'])
            validate_table(old/'fort.63',abundance,[ts[i]],rs)
            if not d.exists():
                shutil.copytree(old,d,symlinks=True)
                (d/'reused.json').write_text(json.dumps({'source_directory':str(old.resolve()),
                    'source_provenance_sha256':digest(a.reuse_opacity/'provenance.json'),
                    'source_receipt_sha256':digest(old/'completed.json'),
                    'note':'Identical physical isotherm; original input paths and source receipt preserved'},indent=2)+'\n')
            if input_fingerprint(prepared['executables']['synspec'],d)!=saved['input_sha256'] or any(
                    digest(d/name)!=expected for name,expected in saved['outputs'].items()):
                raise ValueError('copied opacity artifacts disagree with original source')
            print(f'reused validated {a.mode} opacity T={ts[i]:.8g}',flush=True)
            return d/'fort.63'
        abundance,_=opacity_inputs(d,prepared,spec,a.hydrogen,a.helium3,ts[i])
        depletion_inputs(d,prepared,spec,a.hydrogen,a.helium3,a.mode)
        execute(Path(prepared['synspec']),d,['fort.63','fort.29'])
        validate_table(d/'fort.63',abundance,[ts[i]],rs)
        print(f'validated {a.mode} opacity T={ts[i]:.8g}',flush=True)
        return d/'fort.63'
    with ThreadPoolExecutor(a.jobs) as pool:paths=list(pool.map(run,indices))
    if indices==list(range(len(ts))):
        merge_isotherms(paths,a.work/'opacity.bin')
        abundance,_=composition(a.hydrogen,a.helium3,spec['metals'])
        validate_table(a.work/'opacity.bin',abundance,ts,rs)
        print('complete',digest(a.work/'opacity.bin'),flush=True)


if __name__=='__main__':main()
