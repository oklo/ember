#!/usr/bin/env python3
"""Build a diagnostic atmosphere source with extra outer-layer heat capacity.

This is a sensitivity experiment, not a grain EOS or production atmosphere.
Density derivatives stay fixed; the adiabatic gradient changes consistently
with the added heat capacity. The perturbation equals the requested cP below
2000 K and tapers smoothly to zero at 2400 K. Original source and inputs must
be retained, and output from this executable is never a canonical grid cell.
"""
import argparse
import difflib
import json
from pathlib import Path
import shutil

from prepare_nongrey_sources import digest,run


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('prepared',type=Path);p.add_argument('work',type=Path)
    p.add_argument('--extra-cp',type=float,required=True,help='erg g^-1 K^-1; positive diagnostic perturbation')
    p.add_argument('--depletion-object',type=Path,required=True)
    p.add_argument('--fastchem-library',type=Path,required=True)
    a=p.parse_args()
    if not 0<a.extra_cp<=1e8:raise ValueError('invalid diagnostic heat capacity')
    base=json.loads(a.prepared.read_text());source=Path(base['tlusty']).parent
    if 'initialization_only' in base or 'depletion' not in base:raise ValueError('canonical depleted source required')
    if digest(base['tlusty'])!=base['executables']['tlusty']:raise ValueError('source executable changed')
    if digest(a.fastchem_library)!=base['depletion']['library_sha256']:raise ValueError('chemistry library changed')
    if a.work.exists():raise FileExistsError('use a new diagnostic source directory')
    a.work.mkdir(parents=True);target=a.work.resolve()/'tlusty';shutil.copytree(source,target)
    original=(source/'tlusty208.f').read_text()
    start=original.index('      SUBROUTINE TRMDER(');end=original.index('\n      END\n',start)
    part=original[start:end];anchor='      RETURN'
    if part.count(anchor)!=1:raise ValueError('unrecognized native thermodynamic routine')
    value=format(a.extra_cp,'.17e').replace('e','D')
    part=part.replace(anchor,f'''C     Diagnostic positive thermal reservoir; not production physics.
      EMBWGT=MAX(0.D0,MIN(1.D0,(2400.D0-T)/400.D0))
      EMBWGT=EMBWGT*EMBWGT*(3.D0-2.D0*EMBWGT)
      EMBCP={value}*EMBWGT
      IF(HEATCP.LE.0.D0) STOP 'INVALID THERMAL CONTROL GAS CP'
      GRDADB=GRDADB*HEATCP/(HEATCP+EMBCP)
      HEATCP=HEATCP+EMBCP
'''+anchor)
    modified=original[:start]+part+original[end:]
    (target/'tlusty208.f').write_text(modified)
    patch=a.work/'thermal-control.patch'
    patch.write_text(''.join(difflib.unified_diff(original.splitlines(True),modified.splitlines(True),
                                              fromfile='a/tlusty/tlusty208.f',tofile='b/tlusty/tlusty208.f')))
    executable=target/Path(base['tlusty']).name
    command=['gfortran','-O2','-g','-fno-automatic','-std=legacy','-fallow-argument-mismatch',
             '-fcheck=bounds','-fbacktrace','-o',str(executable),'tlusty208.f',
             str(a.depletion_object.resolve()),str(a.fastchem_library.resolve()),'-lc++']
    run(command,target,a.work/'build.log')
    receipt={**base,'tlusty':str(executable),'executables':{**base['executables'],'tlusty':digest(executable)},
             'initialization_only':{'method':'outer-layer heat-capacity sensitivity control; not an accepted atmosphere',
                'extra_cp_erg_g_K':a.extra_cp,'temperature_taper_K':[2000,2400],
                'canonical_prepared':base,'base_source_sha256':digest(source/'tlusty208.f'),
                'patched_source_sha256':digest(target/'tlusty208.f'),'patch_sha256':digest(patch),
                'depletion_object_sha256':digest(a.depletion_object),'command':command}}
    (a.work/'prepared.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(a.work/'prepared.json')


if __name__=='__main__':main()
