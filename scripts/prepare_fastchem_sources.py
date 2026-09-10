#!/usr/bin/env python3
"""Build the pinned offline FastChem condensate audit adapter (GPL).

The adapter is a separate source-generation tool, never linked into Ember.
No network or Fortran/C++ source build is required by ordinary Ember builds.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import urllib.request

COMMIT='ae67cbd559bc64a3233a1cee6030b8e6b50520de'
ARCHIVE_SHA='6e7c04d8e9d7e8f79b929e6865f5d06d5f5fc4077b5fae7e6ab39675a75706de'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('work',type=Path);p.add_argument('--offline',action='store_true')
    p.add_argument('--cpu',default='');p.add_argument('--jobs',type=int,default=4)
    p.add_argument('--compiler',default='c++');a=p.parse_args()
    if not 1<=a.jobs<=8 or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in a.cpu):
        raise ValueError('invalid build configuration')
    root=a.work.resolve();root.mkdir(parents=True,exist_ok=True);archive=root/'source.tar.gz'
    url=f'https://codeload.github.com/NewStrangeWorlds/FastChem/tar.gz/{COMMIT}'
    if not archive.exists():
        if a.offline:raise FileNotFoundError('offline source archive missing')
        with urllib.request.urlopen(url,timeout=120) as response:archive.write_bytes(response.read())
    if digest(archive)!=ARCHIVE_SHA:raise ValueError('unrecognized FastChem source archive')
    source=root/f'FastChem-{COMMIT}'
    # Always extract a clean source copy into a distinct verified directory;
    # upstream CMake writes an executable beside its source files.
    verified=root/'verified-source'
    if verified.exists():shutil.rmtree(verified)
    verified.mkdir()
    with tarfile.open(archive) as tar:tar.extractall(verified,filter='data')
    source=verified/source.name;build=root/'verified-build';build.mkdir(exist_ok=True)
    flags=['-O3','-DNDEBUG']+([f'-mcpu={a.cpu}'] if a.cpu else [])
    commands=[['cmake','-S',str(source),'-B',str(build),'-DCMAKE_BUILD_TYPE=Release',
               '-DCMAKE_CXX_COMPILER='+a.compiler,'-DCMAKE_CXX_FLAGS_RELEASE='+' '.join(flags)],
              ['cmake','--build',str(build),'--target','fastchem_lib','-j',str(a.jobs)],
              [a.compiler,'-std=c++17',*flags,'-I'+str(source/'fastchem_src'),
               str(Path(__file__).with_name('fastchem_probe.cpp').resolve()),str(build/'libfastchem_lib.a'),
               '-o',str(root/'verified-probe')]]
    with (root/'verified-build.log').open('w') as log:
        for command in commands:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    receipt={'repository':'https://github.com/NewStrangeWorlds/FastChem','commit':COMMIT,
        'archive_sha256':ARCHIVE_SHA,'source':str(source),'executable':str(root/'verified-probe'),
        'executable_sha256':digest(root/'verified-probe'),'commands':commands,
        'adapter_sha256':digest(Path(__file__).with_name('fastchem_probe.cpp')),
        'chemistry_data_sha256':{name:digest(source/'input/logK'/name) for name in ['logK.dat','logK_condensates.dat']}}
    (root/'prepared.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(root/'prepared.json')

if __name__=='__main__':main()
