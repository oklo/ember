#!/usr/bin/env python3
"""Build an isolated TLUSTY executable with a larger opacity temperature array.

Only MTABT changes; the atmosphere equations and numerical settings retain
their canonical source. A same-table replay must check the new executable
before interpreting a refined-table comparison.
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
    p.add_argument('--temperatures',type=int,required=True)
    p.add_argument('--depletion-object',type=Path);p.add_argument('--fastchem-library',type=Path)
    a=p.parse_args();base=json.loads(a.prepared.read_text());source=Path(base['tlusty']).parent
    if 'initialization_only' in base or 'array_capacity' in base:raise ValueError('requires the original canonical source')
    if digest(base['tlusty'])!=base['executables']['tlusty']:raise ValueError('canonical executable changed')
    original=(source/'BASICS.FOR').read_text();pattern=r'(\bMTABT\s*=\s*)(\d+)'
    matches=list(re.finditer(pattern,original))
    if len(matches)!=1:raise ValueError('unrecognized temperature array parameter')
    old=int(matches[0][2])
    if not old<a.temperatures<=100:raise ValueError('capacity must increase within the original reader limit of 100')
    if a.work.exists():raise FileExistsError('use a fresh source directory')
    extra=[]
    if 'depletion' in base:
        if not a.depletion_object or not a.fastchem_library:raise ValueError('requires original depletion object and library')
        if digest(a.fastchem_library)!=base['depletion']['library_sha256']:raise ValueError('FastChem library changed')
        extra=[str(a.depletion_object.resolve()),str(a.fastchem_library.resolve()),'-lc++']
    elif a.depletion_object or a.fastchem_library:raise ValueError('gas source has no depletion objects')
    a.work.mkdir(parents=True);target=a.work.resolve()/'tlusty';shutil.copytree(source,target)
    modified=re.sub(pattern,lambda m:m[1]+str(a.temperatures),original)
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
    variant={'canonical_prepared':base,'parameter':'MTABT','original':old,'value':a.temperatures,
             'unchanged_source_sha256':unchanged,'original_parameter_file_sha256':digest(source/'BASICS.FOR'),
             'parameter_file_sha256':digest(target/'BASICS.FOR'),'patch_sha256':digest(patch),'command':command}
    if a.depletion_object:variant['depletion_object_sha256']=digest(a.depletion_object)
    receipt={**base,'tlusty':str(exe),'executables':{**base['executables'],'tlusty':digest(exe)},'array_capacity':variant}
    (a.work/'prepared.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(a.work/'prepared.json')


if __name__=='__main__':main()
