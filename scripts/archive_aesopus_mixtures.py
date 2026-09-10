#!/usr/bin/env python3
"""Preserve the original AESOPUS members used by the extended family, offline."""
import argparse
import hashlib
from pathlib import Path
from zipfile import ZipFile,ZipInfo,ZIP_DEFLATED
from import_aesopus import SOURCE_SHA256

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('archive',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    if hashlib.sha256(a.archive.read_bytes()).hexdigest()!=SOURCE_SHA256:raise ValueError('unrecognized original archive')
    with ZipFile(a.archive) as src,ZipFile(a.output,'w',compression=ZIP_DEFLATED,compresslevel=9) as dst:
        for z in [.01,.02,.03]:
            for x in [0,.1,.2,.35,.5,.7,.8,.9,.95]:
                for region in ['lowT','highT','highR']:
                    name=f'03-GS98/GS98_a0.0_OPALZ_{region}/aesopus2.0_gasbroad_GS98_Z{z:.6f}_X{x:g}.tab'
                    info=ZipInfo(name,(1980,1,1,0,0,0));info.compress_type=ZIP_DEFLATED;info.external_attr=0o100644<<16;info.create_system=3
                    dst.writestr(info,src.read(name),compress_type=ZIP_DEFLATED,compresslevel=9)
    print(hashlib.sha256(a.output.read_bytes()).hexdigest())
if __name__=='__main__':main()
