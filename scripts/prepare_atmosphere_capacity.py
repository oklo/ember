#!/usr/bin/env python3
"""Build an isolated TLUSTY executable with adjusted array capacities.

Increase MTABT (opacity temperatures) or MDEPTH (atmospheric depths). An
optional MFREQ limit reduces unused storage while retaining the requested
frequency resolution. The equations and numerical settings retain their
original source. A replay must check the executable before using its results.
"""
import argparse
import difflib
import json
from pathlib import Path
import re
import shutil

from prepare_nongrey_sources import digest,run


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('prepared',type=Path);p.add_argument('work',type=Path)
    capacity=p.add_mutually_exclusive_group(required=True)
    capacity.add_argument('--temperatures',type=int)
    capacity.add_argument('--depths',type=int)
    p.add_argument('--frequency-capacity',type=int,
                   help='optional smaller MFREQ capacity; does not change model frequency sampling')
    p.add_argument('--depletion-object',type=Path);p.add_argument('--fastchem-library',type=Path)
    a=p.parse_args();base=json.loads(a.prepared.read_text());source=Path(base['tlusty']).parent
    if 'initialization_only' in base or 'array_capacity' in base:raise ValueError('requires the original canonical source')
    if digest(base['tlusty'])!=base['executables']['tlusty']:raise ValueError('canonical executable changed')
    parameter='MTABT' if a.temperatures is not None else 'MDEPTH'
    value=a.temperatures if a.temperatures is not None else a.depths
    original=(source/'BASICS.FOR').read_text();pattern=r'(\b'+parameter+r'\s*=\s*)(\d+)'
    matches=list(re.finditer(pattern,original))
    if len(matches)!=1:raise ValueError('unrecognized array capacity parameter')
    old=int(matches[0][2])
    if not old<value:raise ValueError('array capacity must increase')
    if parameter=='MTABT' and value>100:raise ValueError('temperature reader limit is 100')
    additional={}
    if a.frequency_capacity is not None:
        fmatches=list(re.finditer(r'(\bMFREQ\s*=\s*)(\d+)',original))
        if len(fmatches)!=1:raise ValueError('unrecognized frequency capacity')
        fold=int(fmatches[0][2])
        if not 2<=a.frequency_capacity<fold:raise ValueError('frequency capacity must decrease and remain >=2')
        additional['MFREQ']={'original':fold,'value':a.frequency_capacity}
    if a.work.exists():raise FileExistsError('use a fresh source directory')
    extra=[]
    if 'depletion' in base:
        if not a.depletion_object or not a.fastchem_library:raise ValueError('requires original depletion object and library')
        if digest(a.fastchem_library)!=base['depletion']['library_sha256']:raise ValueError('FastChem library changed')
        extra=[str(a.depletion_object.resolve()),str(a.fastchem_library.resolve()),'-lc++']
    elif a.depletion_object or a.fastchem_library:raise ValueError('gas source has no depletion objects')
    a.work.mkdir(parents=True);target=a.work.resolve()/'tlusty';shutil.copytree(source,target)
    modified=re.sub(pattern,lambda m:m[1]+str(value),original)
    for name,change in additional.items():
        modified=re.sub(r'(\b'+name+r'\s*=\s*)(\d+)',
                        lambda m:m[1]+str(change['value']),modified)
    (target/'BASICS.FOR').write_text(modified)
    patch=a.work/'array-capacity.patch'
    patch.write_text(''.join(difflib.unified_diff(original.splitlines(True),modified.splitlines(True),
                    fromfile='a/tlusty/BASICS.FOR',tofile='b/tlusty/BASICS.FOR')))
    files=[f for f in source.iterdir() if f.suffix.lower() in ['.f','.for']]
    unchanged={f.name:digest(f) for f in files if f.name!='BASICS.FOR'}
    if any(digest(target/name)!=value for name,value in unchanged.items()):raise ValueError('equation source changed')
    exe=target/Path(base['tlusty']).name
    command=['gfortran','-O2','-g','-fno-automatic','-std=legacy','-fallow-argument-mismatch',
             '-fcheck=bounds','-fbacktrace','-o',str(exe),'tlusty208.f',*extra]
    run(command,target,a.work/'build.log')
    variant={'canonical_prepared':base,'parameter':parameter,'original':old,'value':value,
             'unchanged_source_sha256':unchanged,'original_parameter_file_sha256':digest(source/'BASICS.FOR'),
             'parameter_file_sha256':digest(target/'BASICS.FOR'),'patch_sha256':digest(patch),'command':command}
    if additional:variant['additional_capacities']=additional
    if a.depletion_object:variant['depletion_object_sha256']=digest(a.depletion_object)
    receipt={**base,'tlusty':str(exe),'executables':{**base['executables'],'tlusty':digest(exe)},'array_capacity':variant}
    (a.work/'prepared.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(a.work/'prepared.json')


if __name__=='__main__':main()
