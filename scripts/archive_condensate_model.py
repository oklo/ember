#!/usr/bin/env python3
"""Archive one accepted depleted atmosphere, including its numerical initializer."""
import argparse
import gzip
import json
from pathlib import Path
from prepare_nongrey_sources import digest
from validate_condensate_model import validate,original

MODEL_FILES=['provenance.json','validated.json','completed.json','fort.5','fort.7','fort.8','fort.9',
    'fort.15','run.log','tas','ember-masses.dat','opacity.sha256','physics.json',
    'ember-condensates.cfg','condensate-abundances.dat','condensates.sha256',
    'chemistry-audit.json','chemistry-audit-input.dat','chemistry-audit-source.json','chemistry-audit.log']
INITIALIZER_FILES=[name for name in MODEL_FILES if not name.startswith('chemistry-') and name!='validated.json']+['initializer.json']
HOT_JOIN_FILES=['chemistry-hot-join-input.dat','chemistry-hot-join-source.json','chemistry-hot-join.log']


def archive_model(work,destination,*,diagnostic_grain_transport=False):
    work=Path(work);destination=Path(destination)
    result=validate(work,diagnostic_grain_transport=diagnostic_grain_transport);hashes={}
    def copy(directory,name,relative):
        source=original(directory,name);target=destination/relative/(name+'.gz')
        data=source.read_bytes() if source.suffix=='.gz' else gzip.compress(source.read_bytes(),compresslevel=9,mtime=0)
        target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists() and target.read_bytes()!=data:raise ValueError('archive already contains different source bytes')
        target.write_bytes(data);hashes[str(target.relative_to(destination))]=digest(target)
    for name in MODEL_FILES:copy(work,name,Path())
    if 'vaporized_continuation' in result['chemistry']:
        for name in HOT_JOIN_FILES:copy(work,name,Path())
    if 'initialization_provenance' in result['provenance']:
        initial=work/'initializer'
        if not initial.exists():initial=Path(result['provenance']['initial'])
        for name in INITIALIZER_FILES:copy(initial,name,Path('initializer'))
    # Reimport from compressed originals, with no dependence on working
    # directory paths for the source atmosphere or its initializer.
    reloaded=validate(destination,diagnostic_grain_transport=diagnostic_grain_transport)
    if reloaded!=result:raise ValueError('compressed archive does not reproduce source validation')
    record={'kind':'Original canonical equilibrium-depleted atmosphere; zero grain opacity',
        'coordinates':result['coordinates'],'diagnostics':result['diagnostics'],
        'chemistry':result['chemistry'],'files_sha256':hashes}
    if diagnostic_grain_transport:
        record['kind']='Diagnostic atmosphere with omitted grain enthalpy; not accepted for a production grid'
    path=destination/'archive.json'
    content=json.dumps(record,indent=2)+'\n'
    if path.exists() and path.read_text()!=content:raise ValueError('existing archive receipt differs')
    path.write_text(content)
    return result,hashes


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('work',type=Path);p.add_argument('archive',type=Path)
    a=p.parse_args();result,hashes=archive_model(a.work,a.archive)
    print(a.archive,len(hashes),'verified source artifacts')
