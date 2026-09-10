#!/usr/bin/env python3
"""Build an optional TLUSTY initializer with a common Newton step limiter.

The usual componentwise temperature clipping can change the direction of
the Newton step and disrupt a nearly adiabatic deep gradient. This variant
scales the entire temperature correction by one positive factor. It changes
only the route to a solution; every resulting atmosphere must be replayed
with its original source and pass the independent final diagnostics.
"""
import argparse
import difflib
import json
from pathlib import Path
import shutil
import subprocess
from prepare_nongrey_sources import digest, run


def global_step(text):
    start = text.index('      SUBROUTINE RYBCHN(CHANGT)')
    end = text.index('\n      END\n', start)
    part = text[start:end]
    anchor = '      IF(ITER.EQ.1) WRITE(9,800)\n'
    if part.count(anchor) != 1:
        raise ValueError('unrecognized Rybicki temperature update')
    # CHANGT is left intact: the source convergence diagnostic continues
    # to measure the undamped Newton correction, as in the original code.
    part = part.replace(anchor, anchor + '''C     One common step factor preserves the Newton direction.
      STEPFC=UN
      DO ID=1,ND
         CHT=CHANGT(ID)/TEMP(ID)
         IF(CHT.GT.DPLP) STEPFC=MIN(STEPFC,DPLP/CHT)
         IF(CHT.LT.DPLM) STEPFC=MIN(STEPFC,DPLM/CHT)
      END DO
''')
    if part.count('         CHAN=CHT\n') != 1:
        raise ValueError('unrecognized local correction limiter')
    part = part.replace('         CHAN=CHT\n', '         CHAN=CHT*STEPFC\n')
    return text[:start] + part + text[end:]


def indexed_russel(text):
    """Preserve update order while replacing a repeated element-list search.

    Counts, rather than booleans, preserve the original number of updates
    even if an element occurs more than once in NELEMX. Molecular sums and
    every floating-point expression retain their original ordering.
    """
    start=text.index('      SUBROUTINE RUSSEL(TEM,PG)')
    end=text.index('\n        END\n',start)
    part=text[start:end]
    declaration='     *            UIIDU2(100)\n'
    if part.count(declaration)!=1:raise ValueError('unrecognized Russell declarations')
    part=part.replace(declaration,declaration+'''      INTEGER EMBER_METAL_COUNT(100)
      DO I=1,100
         EMBER_METAL_COUNT(I)=0
      END DO
      DO I=1,NMETAL
         NELEMI=NELEMX(I)
         EMBER_METAL_COUNT(NELEMI)=EMBER_METAL_COUNT(NELEMI)+1
      END DO
''')
    search='''              DO I=1,NMETAL
                 NELEMI=NELEMX(I)
                 IF(NELEMJ.EQ.NELEMI) THEN
                    FX(NELEMI)=FX(NELEMI)+ATOMJ*PMOLJ
                    DFX(NELEMI)=DFX(NELEMI)+ATOMJ**2*
     *                          PMOLJ/P(NELEMI)
                 END IF
              END DO'''
    indexed='''              DO I=1,EMBER_METAL_COUNT(NELEMJ)
                 NELEMI=NELEMJ
                    FX(NELEMI)=FX(NELEMI)+ATOMJ*PMOLJ
                    DFX(NELEMI)=DFX(NELEMI)+ATOMJ**2*
     *                          PMOLJ/P(NELEMI)
              END DO'''
    if part.count(search)!=1:raise ValueError('unrecognized Russell element search')
    return text[:start]+part.replace(search,indexed)+text[end:]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('prepared', type=Path)
    p.add_argument('work', type=Path)
    p.add_argument('--depletion-object', type=Path)
    p.add_argument('--fastchem-library', type=Path)
    p.add_argument('--opacity-step', type=float,
                   help='optional smaller forward temperature step for the initializer opacity derivative (1e-6..1e-2)')
    p.add_argument('--local-step',action='store_true',
                   help='retain the original per-layer limiter; isolate the opacity-derivative change')
    p.add_argument('--convective-refinement',action='store_true',
                   help='only repair the native CONREF initializer bounds; retain the original Newton limiter')
    p.add_argument('--indexed-russel',action='store_true',
                   help='replace the chemistry element-list search with a count lookup, preserving arithmetic order; initializer only')
    a = p.parse_args()
    if a.local_step and a.opacity_step is None and not a.convective_refinement:
        raise ValueError('local-step control requires an opacity-derivative change')
    base = json.loads(a.prepared.read_text())
    source = Path(base['tlusty']).parent
    if digest(base['tlusty']) != base['executables']['tlusty']:
        raise ValueError('source executable identity changed')
    if a.work.exists():
        raise FileExistsError('use a fresh initializer directory')
    extra = []
    if 'depletion' in base:
        if not a.depletion_object or not a.fastchem_library:
            raise ValueError('depletion initializer requires its adapter and library')
        if digest(a.fastchem_library) != base['depletion']['library_sha256']:
            raise ValueError('FastChem library identity changed')
        extra = [str(a.depletion_object.resolve()), str(a.fastchem_library.resolve()), '-lc++']
    elif a.depletion_object or a.fastchem_library:
        raise ValueError('gas initializer does not use depletion objects')
    a.work.mkdir(parents=True)
    target = a.work.resolve() / 'tlusty'
    shutil.copytree(source, target)
    original = (source / 'tlusty208.f').read_text()
    if a.convective_refinement:
        patch_source=Path(__file__).resolve().parents[1]/'data/atmosphere/sources/tlusty208-conref.patch'
        subprocess.run(['patch','--batch','-p1','-i',str(patch_source)],cwd=a.work.resolve(),check=True,
                       capture_output=True,text=True)
        modified=(target/'tlusty208.f').read_text()
    else:
        modified = original if a.local_step else global_step(original)
    if a.indexed_russel:modified=indexed_russel(modified)
    if a.opacity_step is not None:
        if not 1e-6<=a.opacity_step<=1e-2:
            raise ValueError('invalid opacity derivative step')
        start=modified.index('      SUBROUTINE OPACTR(IJ)')
        end=modified.index('\n      END\n',start)
        part=modified[start:end]
        anchor='      PARAMETER (DELT=1.D-2)'
        if part.count(anchor)!=1:
            raise ValueError('unrecognized Rybicki opacity derivative step')
        part=part.replace(anchor,'      PARAMETER (DELT='+format(a.opacity_step,'.12e').replace('e','D')+')')
        modified=modified[:start]+part+modified[end:]
    (target / 'tlusty208.f').write_text(modified)
    patch = a.work / 'global-step.patch'
    patch.write_text(''.join(difflib.unified_diff(original.splitlines(True), modified.splitlines(True),
        fromfile='a/tlusty/tlusty208.f', tofile='b/tlusty/tlusty208.f')))
    exe = target / Path(base['tlusty']).name
    command = ['gfortran', '-O2', '-g', '-fno-automatic', '-std=legacy',
               '-fallow-argument-mismatch', '-fcheck=bounds', '-fbacktrace',
               '-o', str(exe), 'tlusty208.f', *extra]
    run(command, target, a.work / 'build.log')
    receipt = {**base, 'tlusty': str(exe),
               'executables': {**base['executables'], 'tlusty': digest(exe)},
               'initialization_only': {'method': ('native CONREF with explicit bounds repairs'
                    if a.convective_refinement else ('original per-layer limiter, smaller opacity derivative step'
                                                if a.local_step else 'common Newton temperature-step scaling')),
                   'canonical_prepared': base, 'command': command,
                   'base_source_sha256': digest(source / 'tlusty208.f'),
                   'patched_source_sha256': digest(target / 'tlusty208.f'),
                   'patch_sha256': digest(patch)}}
    receipt['initialization_only']['opacity_derivative_relative_step']=a.opacity_step or .01
    receipt['initialization_only']['indexed_russel_lookup']=a.indexed_russel
    if a.depletion_object:
        receipt['initialization_only']['depletion_object_sha256'] = digest(a.depletion_object)
    (a.work / 'prepared.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(a.work / 'prepared.json')


if __name__ == '__main__':
    main()
