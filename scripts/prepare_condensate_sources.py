#!/usr/bin/env python3
"""Build separate offline TLUSTY/SYNSPEC settled-grain source experiments.

Requires already verified atmosphere and FastChem sources. Normal Ember
builds neither download nor link FastChem. This experiment combines FastChem
equilibrium elemental removal with the existing gas partition/line data;
grain opacity, grain enthalpy and time-dependent settling are not included.
"""
import argparse
import difflib
import json
from pathlib import Path
import shutil
import subprocess
from prepare_nongrey_sources import digest, run

ROOT=Path(__file__).resolve().parents[1]


def patch_molecules(text,kind):
    begin=text.index('      subroutine moleq(')
    end=text.index('\n      END\n',begin)+len('\n      END\n') if kind=='tlusty' else text.index('\n        END\n',begin)+len('\n        END\n')
    part=text[begin:end]
    # Original bulk numbers must survive calls at perturbed temperatures.
    anchor='      data nmetal/92/'
    if part.count(anchor)!=1:raise ValueError('unrecognized molecular source')
    part=part.replace(anchor,'      dimension cbase(100)\n      integer icond\n'+anchor)
    part=part.replace('ccomp(ia)=abndd(ia,id)','ccomp(ia)=abndd(ia,id)\n'+(' '*9 if kind=='synspec' else ' '*12)+'cbase(ia)=abndd(ia,id)')
    anchor='c---- end of reading atomic and molecular data  ----------------------'
    if part.count(anchor)!=1:raise ValueError('unrecognized molecular initialization')
    hook='''c
c     Equilibrium gas depletion, with grains removed from opacity.
c     Re-evaluate at every trial T/P, including derivative calls.
      call ember_condense(tt,pgas,cbase,ccomp,nmetal,icond)
      if(icond.ne.0) stop 'EMBER DEPLETION FAILED'
      if(ipri.ne.0) then
         ytot(id)=0.d0
         wmy(id)=0.d0
         do i=1,nmetal
            abndd(i,id)=ccomp(i)
            if(iatex(i).ge.0) then
               ytot(id)=ytot(id)+ccomp(i)
               wmy(id)=wmy(id)+ccomp(i)*amas(i)
            end if
         end do
      end if
'''
    part=part.replace(anchor,anchor+'\n'+hook)
    result=text[:begin]+part+text[end:]
    if kind=='tlusty':
        # RUSSEL's current mixture must not come from depth 1 or from the
        # last unperturbed state when numerical T/P derivatives are taken.
        result=result.replace('        HEH=YTOT(1)-UN','        HEH=-UN\n        do i=1,nmetal\n           HEH=HEH+CCOMP(i)/CCOMP(1)\n        end do')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('atmosphere',type=Path);p.add_argument('fastchem',type=Path)
    p.add_argument('library',type=Path);p.add_argument('work',type=Path)
    a=p.parse_args();base=json.loads(a.atmosphere.read_text());fc=json.loads(a.fastchem.read_text())
    if a.work.exists():raise FileExistsError('use a fresh source directory')
    a.work.mkdir(parents=True);a.work=a.work.resolve()
    if digest(fc['executable'])!=fc['executable_sha256']:raise ValueError('FastChem executable identity changed')
    for name,sha in fc['chemistry_data_sha256'].items():
        if digest(Path(fc['source'])/'input/logK'/name)!=sha:raise ValueError('FastChem chemistry changed')
    adapter=ROOT/'scripts/fastchem_depletion.cpp';obj=a.work/'depletion.o'
    run(['c++','-std=c++17','-O3','-mcpu=apple-m4','-I'+str(Path(fc['source'])/'fastchem_src'),
         '-c',str(adapter),'-o',str(obj)],a.work,a.work/'adapter-build.log')
    receipt={**base,'depletion':{'fastchem':fc,'adapter_sha256':digest(adapter),
        'library_sha256':digest(a.library),'patches':{},'source_sha256':{}}}
    for kind in ['tlusty','synspec']:
        oldexe=Path(base[kind]);source=oldexe.parent
        if digest(oldexe)!=base['executables'][kind]:raise ValueError('atmosphere executable changed')
        target=a.work/kind;shutil.copytree(source,target)
        name='tlusty208.f' if kind=='tlusty' else 'synspec54.f'
        original=(source/name).read_text();modified=patch_molecules(original,kind)
        (target/name).write_text(modified)
        patch=a.work/(kind+'-depletion.patch')
        patch.write_text(''.join(difflib.unified_diff(original.splitlines(True),modified.splitlines(True),
            fromfile='a/'+kind+'/'+name,tofile='b/'+kind+'/'+name)))
        exe=target/oldexe.name
        run(['gfortran','-O2','-g','-fno-automatic','-std=legacy','-fallow-argument-mismatch',
             '-fcheck=bounds','-fbacktrace','-o',str(exe),name,str(obj),str(a.library.resolve()),'-lc++'],
             target,a.work/(kind+'-build.log'))
        receipt[kind]=str(exe);receipt['executables'][kind]=digest(exe)
        receipt['depletion']['patches'][kind]=digest(patch)
        receipt['depletion']['source_sha256'][kind]={'base':digest(source/name),'patched':digest(target/name)}
    (a.work/'prepared.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(a.work/'prepared.json')


if __name__=='__main__':main()
