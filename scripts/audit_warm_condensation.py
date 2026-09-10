#!/usr/bin/env python3
"""Check condensates on verified warm gas-atmosphere profiles.

This is a fixed-structure equilibrium diagnostic, not an atmosphere with
grain feedback. The independent chemistry query stops at 6000 K and uses
the existing explicit vaporized-join check for hotter atmosphere layers.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import json
import math
from pathlib import Path
import subprocess
import numpy as np

from assemble_nongrey_grid import load_continuation
from condensate_chemistry_profile import audit_profile, HOT_JOIN_TEMPERATURE
from generate_nongrey_grid import composition
from import_nongrey_grid import read_text
from prepare_nongrey_sources import digest, SOURCES


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['assembly','source','work','output']:p.add_argument(name,type=Path)
    p.add_argument('--temperatures',type=float,nargs='+',default=[3200.,3400.])
    p.add_argument('--hydrogen',type=float,nargs='+',default=[.1,.15,.2])
    p.add_argument('--logg',type=float,default=5.15)
    p.add_argument('--jobs',type=int,default=2)
    a=p.parse_args()
    if not 1<=a.jobs<=8:raise ValueError('invalid worker count')
    a.work.mkdir(parents=True,exist_ok=False)
    source=json.loads(a.source.read_text());assembly=json.loads(a.assembly.read_text())
    executable=Path(source['executable'])
    if digest(executable)!=source['executable_sha256']:raise ValueError('chemistry executable changed')
    metadata=json.loads((SOURCES/'synple-elements.json').read_text())
    selections=[e for e in assembly['extension'] if e['coordinates'][0] in a.hydrogen
                and e['coordinates'][2] in a.temperatures and e['coordinates'][3]==a.logg]
    if len(selections)!=len(a.hydrogen)*len(a.temperatures)*2:
        raise ValueError('missing or duplicate requested warm atmosphere')

    def calculate(entry):
        work=Path(entry['work']);key,state,spec,record,dependencies=load_continuation(work)
        if list(key)!=entry['coordinates'] or digest(work/'final/validated.json')!=entry['validation_sha256']:
            raise ValueError('assembly atmosphere identity differs')
        logfile=work/record['log'];profile=[]
        for line in read_text(logfile).rsplit('FINAL MODEL ATMOSPHERE',1)[1].splitlines():
            words=line.replace('D','E').split()
            if len(words)==11 and words[0].isdigit():profile.append(list(map(float,words)))
        profile.reverse()
        abundance,weights=composition(*key[:2],spec['metals']);numbers={};masses={}
        for symbol,n,mass in zip(metadata['symbol'],abundance,weights,strict=True):
            if n>1e-90:
                numbers[symbol.capitalize()]=n
                masses[symbol.capitalize()]=mass*1.67333e-24/1.66053906660e-24
        label=f'x{key[0]:g}-y{key[1]:g}-t{key[2]:g}-g{key[3]:g}'
        folder=a.work/label;folder.mkdir()
        abundances=folder/'abundances.dat'
        abundances.write_text('# Shared atmosphere element numbers\ne- 0\n'+''.join(
            f'{s} {12+math.log10(n):.17g}\n' for s,n in numbers.items()))
        def query(rows,name):
            inputs=''.join(f'{r[3]:.17g} {r[6]/1e6:.17g}\n' for r in rows)
            (folder/(name+'.input')).write_text(inputs)
            result=subprocess.run([str(executable),source['source'],str(abundances.resolve()),'equilibrium'],
                                  input=inputs,text=True,capture_output=True,timeout=600)
            (folder/(name+'.log')).write_text(result.stderr)
            result.check_returncode()
            (folder/(name+'.json.gz')).write_bytes(gzip.compress(result.stdout.encode(),mtime=0))
            return json.loads(result.stdout)
        cool=[row for row in profile if row[3]<=HOT_JOIN_TEMPERATURE]
        hot=[row.copy() for row in profile if row[3]>HOT_JOIN_TEMPERATURE]
        for row in hot:row[3]=HOT_JOIN_TEMPERATURE
        chemistry=query(cool,'cool');joined=query(hot,'hot-join') if hot else None
        audit,condensed=audit_profile(chemistry,profile,numbers,masses,joined)
        result={'coordinates':list(key),'atmosphere':str(work),'profile_layers':len(profile),
                'layers_with_condensates':int(np.count_nonzero(condensed)),'chemistry':audit,
                'atmosphere_input_sha256':dependencies,
                'diagnostic_files_sha256':{str(f):digest(f) for f in folder.iterdir() if f.is_file()}}
        (folder/'audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print(json.dumps({'coordinates':list(key),'condensed_layers':result['layers_with_condensates'],
                         'maximum_condensed_mass_fraction':audit['thresholds']['0']['max_condensed_mass_fraction'],
                         'condensates':audit['condensates']}),flush=True)
        return result

    with ThreadPoolExecutor(a.jobs) as pool:results=list(pool.map(calculate,selections))
    report={'scope':__doc__,'chemistry_source':source,'input_sha256':{str(f):digest(f) for f in [a.assembly,a.source,Path(__file__)]},
            'helper_sha256':{name:digest(Path(__file__).with_name(name)) for name in ['audit_condensates.py','condensate_chemistry_profile.py']},
            'records':results}
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
