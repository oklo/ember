#!/usr/bin/env python3
"""Build a pinned, offline dense-plasma EOS reference; no runtime installation."""
import argparse
import difflib
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import urllib.request

URL='https://www.ioffe.ru/astro/EIP/eos22.f'
SOURCE_SHA='f514e1781f4f476c8e23ad03848488445aa539dd7ccbdda049dd9ed91d7f8ae9'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('work',type=Path);p.add_argument('--source',type=Path)
    p.add_argument('--compiler',default='gfortran')
    p.add_argument('--phase-probe',action='store_true',
                   help='audit EOSFI22 at explicitly selected phase and exact density; excludes ideal electrons')
    p.add_argument('--repair-quantum-density-derivative',action='store_true',
                   help='isolated EOSFI22 repair: use PDRQL in its density derivative; preserve original source and patch')
    p.add_argument('--repair-quantum-pressure-derivatives',action='store_true',
                   help='also differentiate the LIQUBC free energy consistently, including density-dependent coefficients')
    a=p.parse_args()
    if a.work.exists():raise FileExistsError('use a fresh source-build directory')
    a.work.mkdir(parents=True)
    source=a.work/'eos22.f'
    if a.source:shutil.copyfile(a.source,source)
    else:
        with urllib.request.urlopen(URL,timeout=60) as response:source.write_bytes(response.read())
    if sha(source)!=SOURCE_SHA:raise ValueError('Ioffe source checksum differs')
    original_sha=sha(source)
    patch=None
    if a.repair_quantum_density_derivative or a.repair_quantum_pressure_derivatives:
        original=source.read_text()
        anchor='         PDRi=PDRii+1.d0+PDTQL'
        if original.count(anchor)!=1:raise ValueError('unrecognized quantum density-derivative assignment')
        modified=original.replace(anchor,'         PDRi=PDRii+1.d0+PDRQL')
        if a.repair_quantum_pressure_derivatives:
            start=modified.index('      subroutine LIQUBC(')
            end=modified.index('      subroutine FSCRliq8(',start)
            part=modified[start:end]
            replacements={
                '           DI=DI1\n':'           DI=DI1\n           DIDR=-DI1*(1.d0-DI1)\n',
                '           DI=0.\n':'           DI=0.\n           DIDR=0.\n',
                '           DI=-(CI1/CI)**2*DI1\n':
                    '           DI=-(CI1/CI)**2*DI1\n           DIDR=DI*(3.d0*DI1-2.d0*DI-1.d0)\n',
                '           PDTQLI=-FQLI\n           PDRQLI=UQLI\n':'',
                '           PDTQLI=.5d0*CVQLI-DI*CVQLI/3.d0 ! BC\'22(38)\n'
                '           PDRQLI=.75d0*UQLI-.25d0*CVQLI-.5d0*DI*(UQLI-CVQLI/3.d0)-\n'
                '     -       RSI/(P3+RSI)*DI/9.d0*UQLI ! BC\'22(39)\n':'',
                '         FQL=FQL+FQLI\n':
                    'c     Differentiate the same f(C_i(RSI)*TPT) potential.\n'
                    'c     AA=d ln(C_i*TPT)/d ln(rho), DIDR=dDI/dln(RSI).\n'
                    '         AA=.5d0-DI/3.d0\n'
                    '         PDTQLI=AA*CVQLI\n'
                    '         PDRQLI=(AA+AA**2+DIDR/9.d0)*UQLI-\n'
                    '     -      AA**2*CVQLI\n'
                    '         FQL=FQL+FQLI\n'}
            for before,after in replacements.items():
                if part.count(before)!=1:raise ValueError('unrecognized LIQUBC derivative source')
                part=part.replace(before,after)
            modified=modified[:start]+part+modified[end:]
        source.rename(a.work/'eos22-original.f')
        source.write_text(modified)
        patch=a.work/'eos22-quantum-derivatives.patch'
        patch.write_text(''.join(difflib.unified_diff(original.splitlines(True),modified.splitlines(True),
                                                    fromfile='a/eos22.f',tofile='b/eos22.f')))
    wrapper=Path(__file__).with_name('eip_phase_probe.f90' if a.phase_probe else 'eip_probe.f90').resolve()
    executable=a.work.resolve()/'eip-probe'
    command=[a.compiler,'-O2','-std=legacy','-fcheck=bounds','-fbacktrace',
             str(wrapper),str(source.resolve()),'-o',str(executable)]
    with (a.work/'build.log').open('w') as log:
        subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    record={'scope':('Locally repaired' if patch else 'Unmodified')+' Ioffe dense-plasma reference, not a runtime EOS',
            'probe_kind':'forced-phase interaction terms' if a.phase_probe else 'native mixture',
            'source_url':URL,'upstream_source_sha256':original_sha,'source_sha256':sha(source),'wrapper_sha256':sha(wrapper),
            'executable':str(executable),'executable_sha256':sha(executable),'command':command,
            'compiler':subprocess.check_output([a.compiler,'--version'],text=True).splitlines()[0],
            'limitations':['Fully ionized matter only; no bound states or positrons.',
                           'Native density inversion is approximate; requested and inferred density are both reported.',
                           'Native phase choice is the classical effective Gamma=175 threshold; quantum melting is not solved.',
                           'Material quantities exclude radiation. Native ionic spin convention is retained.',
                           'No phase separation, diffusion or integration with a stellar model.']}
    if a.phase_probe:
        record['limitations']=['Single fully ionized species; includes ideal ions and ii/ie/ee interactions.',
                              'Excludes ideal electrons and radiation, so its pressure need not be positive.',
                              'Phase is forced by caller; no phase equilibrium, melting or domain acceptance implied.',
                              'Thermodynamic consistency is not an accuracy bound on the underlying physical fits.']
    if patch:record['local_repair']={'patch_sha256':sha(patch),
                                   'repair_LIQUBC':a.repair_quantum_pressure_derivatives,
                                   'reason':'EOSFI22 liquid PDRi uses PDTQL instead of PDRQL. Optional LIQUBC repair differentiates its free energy with the full density coefficient chain rule.',
                                   'validation':'Compare audit_eip_derivatives.py on original and repaired phase probes.'}
    (a.work/'source_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    print(executable)


if __name__=='__main__':main()
