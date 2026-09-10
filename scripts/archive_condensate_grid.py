#!/usr/bin/env python3
"""Archive and import a complete independently validated condensate family."""
import argparse
import gzip
import itertools
import json
from pathlib import Path
import shutil
import tempfile
from import_nongrey_grid import import_grid,CONDENSATE_CALCULATION
from prepare_nongrey_sources import digest
from nongrey_opacity import validate_table
from generate_nongrey_grid import composition,temperatures,sequence
from validate_condensate_model import validate,original,PHYSICAL_KEYS
from archive_condensate_model import archive_model


def archive_family(prepared,specification,plan,models,destination):
    destination=Path(destination).resolve()
    if destination.exists():raise FileExistsError(destination)
    spec=json.loads(Path(specification).read_text());source=json.loads(Path(prepared).read_text())
    plan=json.loads(Path(plan).read_text());models=json.loads(Path(models).read_text())
    expected=set(itertools.product(spec['hydrogen'],spec['helium3'],spec['teff_K'],spec['log_g']))
    supplied=[tuple(r[k] for k in ['XH','X3','teff_K','log_g']) for r in models]
    if len(supplied)!=len(expected) or set(supplied)!=expected:raise ValueError('incomplete or duplicate condensate family')
    opacity=[]
    for x,y in itertools.product(spec['hydrogen'],spec['helium3']):
        planes=[p for p in plan['opacity_planes'] if (p['XH'],p['X3'])==(x,y)]
        if len(planes)!=1:raise ValueError('missing opacity plane')
        work=Path(planes[0]['directory']);provenance=json.loads((work/'provenance.json').read_text())
        if (any(provenance.get(k)!=v for k,v in [('prepared',source),('XH',x),('X3',y),('mode','equilibrium')]) or
                any(provenance['specification'].get(k)!=spec.get(k) for k in PHYSICAL_KEYS)):
            raise ValueError('opacity plane source identity differs')
        abundance,_=composition(x,y,spec['metals'])
        table=validate_table(work/'opacity.bin',abundance,temperatures(spec),sequence(spec['log_density']))
        archive=Path(planes[0]['archive'])
        receipt=json.loads((archive/'archive.json').read_text())
        if receipt['mode']!='equilibrium' or receipt['opacity_sha256']!=digest(work/'opacity.bin'):
            raise ValueError('opacity archive disagrees with model input')
        for name,expected_hash in receipt['files_sha256'].items():
            if digest(archive/name)!=expected_hash:raise ValueError('opacity source archive changed')
        opacity.append({'XH':x,'X3':y,'sha256':receipt['opacity_sha256'],
            'shape_frequency_density_temperature':table['shape'],
            'source_archive':str(archive.resolve()),'source_archive_receipt_sha256':digest(archive/'archive.json')})
    destination.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent) as temp:
        root=Path(temp)/'archive';root.mkdir();hashes={};records=[]
        def copy(path,target):
            target.parent.mkdir(parents=True,exist_ok=True)
            if path.suffix=='.gz':shutil.copyfile(path,target)
            else:target.write_bytes(gzip.compress(path.read_bytes(),compresslevel=9,mtime=0))
            hashes[str(target.relative_to(root))]=digest(target)
        for i,model in enumerate(sorted(models,key=lambda r:tuple(r[k] for k in ['XH','X3','teff_K','log_g']))):
            key=tuple(model[k] for k in ['XH','X3','teff_K','log_g']);work=Path(model['directory'])
            plane=next(p for p in opacity if (p['XH'],p['X3'])==key[:2])
            result=validate(work,spec,source,key,plane['sha256'])
            relative=Path(f'model-{i:03d}');target=root/relative
            _,model_hashes=archive_model(work,target)
            hashes.update({str(relative/name):value for name,value in model_hashes.items()})
            hashes[str(relative/'archive.json')]=digest(target/'archive.json')
            records.append({**dict(zip(['XH','X3','teff_K','log_g'],key)),
                            'condensate_model':str(relative),'diagnostics':result['diagnostics'],'chemistry':result['chemistry']})
        for i,plane in enumerate(opacity):
            copy(Path(plane['source_archive'])/'archive.json',root/f'opacity-{i:03d}-archive.json.gz')
        copy(Path(prepared),root/'prepared.json.gz');copy(Path(specification),root/'specification.json.gz')
        manifest={**spec,'format':1,'calculation':CONDENSATE_CALCULATION,'provenance':source,
            'condensates':{'mode':'equilibrium','grain_opacity':0,
                'grain_enthalpy':'reject condensing layers with convective flux fraction above 1e-8'},
            'opacity':{'temperature_K':temperatures(spec),'density_g_cm3':sequence(spec['log_density'])},
            'opacity_planes':opacity,'models':records,'archive_files_sha256':hashes}
        (root/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        import_grid(root/'manifest.json',root/'atmosphere.dat')
        root.rename(destination)
    return destination/'atmosphere.dat'


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['prepared','specification','plan','models','destination']:p.add_argument(name,type=Path)
    a=p.parse_args();print(archive_family(a.prepared,a.specification,a.plan,a.models,a.destination))
