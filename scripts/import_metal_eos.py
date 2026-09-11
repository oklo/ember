#!/usr/bin/env python3
"""Validate and import the complete baryonic GS98 FreeEOS H/He3 family."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
from fetch_tops_composition import fraction_label
from eos_source_coverage import absent_source_rows, inconsistent_source_rows

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('work',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--jobs',type=int,default=4)
    a=p.parse_args()
    if not 1<=a.jobs<=8:raise ValueError('invalid jobs')
    spec=json.loads((a.work/'specification.json').read_text())
    a.output.mkdir(parents=True,exist_ok=True);archive=a.output/'sources';archive.mkdir(exist_ok=True)
    def run(job):
        i,x,y=job;source=a.work/f'plane-{i:03d}'/'source.json.gz'
        raw=json.loads(gzip.decompress(source.read_bytes()))
        if raw['hydrogen']!=x or raw['helium3']!=y or raw['logT']!=spec['logT'] or raw['logQ']!=spec['logQ']:
            raise ValueError('plane metadata does not match family specification')
        if raw['probe_sha256']!=spec['probe_sha256'] or raw['source_archive_sha256']!=spec['source_archive_sha256']:
            raise ValueError('inconsistent source identity')
        if raw['options']!=spec.get('source_options',[3,1,-2]):
            raise ValueError('source physics differs from family specification')
        if raw.get('source_coverage')!=spec.get('source_coverage'):
            raise ValueError('source coverage differs from family specification')
        if raw.get('precision_fallback')!=spec.get('precision_fallback'):
            raise ValueError('source numerical precision differs from family specification')
        absent=absent_source_rows(raw)
        excluded=inconsistent_source_rows(raw)
        if excluded and spec.get('source_consistency_exclusion_limit')!=1e-7:
            raise ValueError('source consistency exclusions are not declared by the family')
        # Refinement may require sub-per-mille abundances. Rounded labels
        # alias distinct physical planes; preserve the exact request value.
        name=f'freeeos300_gs98_x{fraction_label(x,1000)}_he3{fraction_label(y,1000)}'
        target=a.output/(name+'_potential.dat');saved=archive/(name+'_raw.json.gz')
        if saved.exists() and digest(saved)!=digest(source):raise ValueError('refusing to replace a different raw source')
        shutil.copy2(source,saved)
        subprocess.run([sys.executable,str(Path(__file__).with_name('import_freeeos_potential.py')),str(saved),str(target)],check=True)
        return {'hydrogen':x,'helium3':y,'potential':target.name,'potential_sha256':digest(target),
                'source':str(saved.relative_to(a.output)),'source_sha256':digest(saved),
                'failed_source_states':sum(r is not None and r[0]!=0 for r in raw['data']),
                'absent_source_states':sum(absent),'inconsistent_source_states':len(excluded)}
    jobs=[(i,x,y) for i,(x,y) in enumerate((x,y) for x in spec['hydrogen'] for y in spec['helium3'])]
    with ThreadPoolExecutor(a.jobs) as pool:planes=list(pool.map(run,jobs))
    lines=['EMBER_METAL_HELMHOLTZ 1','hydrogen '+str(len(spec['hydrogen']))+' '+' '.join(map(str,spec['hydrogen'])),
           'helium3 '+str(len(spec['helium3']))+' '+' '.join(map(str,spec['helium3']))]
    lines += [json.dumps(r['potential']) for r in planes]
    family=a.output/'freeeos300_gs98_z020.dat';family.write_text('\n'.join(lines)+'\n')
    provenance={**spec,'planes':planes,'family_sha256':digest(family),
        'scripts':{name:digest(Path(__file__).with_name(name)) for name in
                   ['generate_metal_eos.py','metal_eos_composition.py','import_freeeos_potential.py','import_metal_eos.py','eos_source_coverage.py']}}
    (archive/'freeeos300_gs98_manifest.json').write_text(json.dumps(provenance,indent=2)+'\n')

if __name__=='__main__':main()
