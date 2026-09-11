#!/usr/bin/env python3
"""Build the offline probe from the checksum-pinned FreeEOS 3.0 archive.

FreeEOS is external GPL-2.0-or-later software. Its source/build remain in
the supplied working directory; ember does not link to it at runtime.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

SHA256='4ab1c15a51385a3eab3b08c6f3f240739c0105d92ec828d635ac95720edefb09'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('archive',type=Path)
    ap.add_argument('work',type=Path)
    ap.add_argument('--compiler',default='gfortran')
    ap.add_argument('--jobs',type=int,choices=range(1,9),default=4)
    ap.add_argument('--electron-quadrature-error',type=float,choices=[1e-9,1e-11],default=1e-9,
                    help='numerical electron integration target; does not change EOS formulas')
    a=ap.parse_args()
    if hashlib.sha256(a.archive.read_bytes()).hexdigest()!=SHA256:
        raise ValueError('unexpected FreeEOS source checksum')
    compiler=shutil.which(a.compiler)
    if compiler is None: raise ValueError('Fortran compiler unavailable')
    work=a.work.resolve();work.mkdir(parents=True,exist_ok=True)
    with tarfile.open(a.archive,'r:gz') as tar:
        tar.extractall(work,filter='data')
    src=work/'free_eos-3.0.0';build=work/'build'
    # On arm64 gfortran, extended and quadruple precision resolve to the
    # same kind. Remove duplicate generic diagnostic overloads in that
    # case. The quadruple routines still handle both; no EOS math changes.
    check=work/'check_kinds.f90'
    check.write_text('program kinds\nprint *, selected_real_kind(18,4900), selected_real_kind(32,4900)\nend program kinds\n')
    subprocess.run([compiler,str(check),'-o',str(work/'check_kinds')],check=True)
    kinds=subprocess.check_output([str(work/'check_kinds')],text=True).split()
    if kinds[0]==kinds[1]:
        f=src/'src/mod_free_eos_debug.f90';s=f.read_text()
        for n in (1,2,3):
            old=f'     module procedure set_snan_ep_fp_kind_{n}d'
            if s.count(old)!=1: raise ValueError('unexpected source diagnostic interface')
            s=s.replace(old,f'     ! ep and qp coincide: qp overload handles both ({n}d).')
        f.write_text(s)
    precision_patch={}
    if a.electron_quadrature_error!=1e-9:
        f=src/'src/fermi_dirac_direct.f90';s=f.read_text()
        original=hashlib.sha256(f.read_bytes()).hexdigest()
        old='real(fp_kind), parameter :: fderr = 1.e-09_fp_kind'
        if s.count(old)!=1:raise ValueError('unexpected numerical electron integration target')
        f.write_text(s.replace(old,'real(fp_kind), parameter :: fderr = 1.e-11_fp_kind'))
        precision_patch={'file':'src/fermi_dirac_direct.f90','original_sha256':original,
                         'modified_sha256':hashlib.sha256(f.read_bytes()).hexdigest()}
    subprocess.run(['cmake','-S',str(src),'-B',str(build),'-G','Ninja',
                    '-DCMAKE_Fortran_COMPILER='+compiler,'-DCMAKE_Fortran_FLAGS=-O3',
                    '-DBUILD_TEST=ON','-DBUILD_DOX_DOC=OFF',
                    '-DCMAKE_INSTALL_PREFIX='+str(work/'install')],check=True)
    subprocess.run(['cmake','--build',str(build),'--target','free_eos','-j',str(a.jobs)],check=True)
    probe=Path(__file__).resolve().with_name('freeeos_probe.f90')
    subprocess.run([compiler,'-O3','-I'+str(build/'src'),str(probe),
                    '-L'+str(build/'src'),'-lfree_eos','-Wl,-rpath,'+str(build/'src'),
                    '-o',str(work/'probe')],check=True)
    library=build/'src/libfree_eos.1.0.0.dylib'
    if not library.exists():library=next((build/'src').glob('libfree_eos.so.*'))
    record={'source_archive_sha256':SHA256,'electron_quadrature_error':a.electron_quadrature_error,
            'precision_patch':precision_patch,'probe_sha256':hashlib.sha256((work/'probe').read_bytes()).hexdigest(),
            'library':str(library),'library_sha256':hashlib.sha256(library.read_bytes()).hexdigest(),
            'builder_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (work/'build_receipt.json').write_text(json.dumps(record,indent=2)+'\n')
    print(work/'probe')


if __name__=='__main__':
    main()
