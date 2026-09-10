#!/usr/bin/env python3
"""Build pinned, optional grain optics tools outside the Ember library.

LX-MIE is GPL-3.0 and Optool is MIT licensed. Their original sources and
optical-data references remain in the supplied work directory. No normal
Ember build downloads, compiles, or links either package.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import urllib.request
from prepare_nongrey_sources import digest, extract, run

ROOT = Path(__file__).resolve().parents[1]
SOURCES = {
    'lxmie': ('NewStrangeWorlds/LX-MIE', 'a078eddbd26d89e7706b75fef6fe4ca08bb271db',
              '4fa99171b0334756dfb9b116073e8b149200573adccc29e28b003edae729369c'),
    'optool': ('cdominik/optool', '0b1fa6df95e3a694dc64abdbf8e36dad127165ab',
               'e271e7118461b6ec6af0a00ac807ce291208779903bfcde59e09af52919c45a5'),
}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('work', type=Path)
    p.add_argument('--archives', type=Path, help='read pinned lxmie.tar.gz and optool.tar.gz locally')
    p.add_argument('--cxx', default='c++')
    a = p.parse_args()
    a.work.mkdir(parents=True, exist_ok=False); work = a.work.resolve()
    receipt = {'format': 1, 'scope': 'Offline spherical-grain optics; no atmosphere installation', 'sources': {}}
    for key, (repository, commit, expected) in SOURCES.items():
        url = f'https://codeload.github.com/{repository}/tar.gz/{commit}'
        archive = work / (key + '.tar.gz')
        if a.archives:
            shutil.copyfile(a.archives / archive.name, archive)
        else:
            with urllib.request.urlopen(url, timeout=60) as response, archive.open('wb') as out:
                shutil.copyfileobj(response, out)
        if digest(archive) != expected:
            raise ValueError(f'{key}: source archive checksum differs')
        source = extract(archive, work / (key + '-source'))
        receipt['sources'][key] = {'repository': 'https://github.com/' + repository,
            'commit': commit, 'archive_url': url, 'archive_sha256': expected, 'source': str(source)}
    source = Path(receipt['sources']['lxmie']['source'])
    adapter = ROOT / 'scripts/grain_mie_probe.cpp'; probe = work / 'grain-mie-probe'
    command = [a.cxx, '-std=c++17', '-O2', '-I' + str(source / 'src/mie'),
               str(adapter), str(source / 'src/mie/mie.cpp'), '-o', str(probe)]
    run(command, work, work / 'lxmie-build.log')
    receipt['probe'] = str(probe); receipt['probe_sha256'] = digest(probe)
    receipt['adapter_sha256'] = digest(adapter); receipt['build_command'] = command
    receipt['compiler'] = subprocess.check_output([a.cxx, '--version'], text=True)
    receipt['optical_data'] = {path.name: {'sha256': digest(path), 'reference': path.read_text().splitlines()[0]}
                             for path in sorted((source / 'compilation').glob('*.dat'))}
    source = Path(receipt['sources']['optool']['source'])
    run(['make', '-j2'], source, work / 'optool-build.log')
    receipt['optool'] = str(source / 'optool'); receipt['optool_sha256'] = digest(source / 'optool')
    (work / 'prepared.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(work / 'prepared.json')


if __name__ == '__main__':
    main()
