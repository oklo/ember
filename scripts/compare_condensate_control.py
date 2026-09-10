#!/usr/bin/env python3
"""Archive and compare depleted and gas atmospheres with matching source grids."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

from archive_condensate_model import archive_model, MODEL_FILES
from generate_nongrey_grid import composition, input_fingerprint, temperatures, sequence
from import_nongrey_grid import read_text, source_inputs, source_state
from prepare_nongrey_sources import digest
from nongrey_opacity import validate_table
from validate_condensate_model import original, validate, PHYSICAL_KEYS


def gas_state(work, reference, opacity=None):
    work=Path(work)
    provenance=json.loads(read_text(original(work,'provenance.json')))
    other=reference['provenance']
    if provenance['mode']!='gas' or provenance['prepared']!=other['prepared']:
        raise ValueError('gas control source or mode differs')
    if any(provenance[k]!=other[k] for k in ['XH','X3','teff_K','log_g']):
        raise ValueError('gas control coordinates differ')
    spec=provenance['specification']
    if any(spec.get(k)!=other['specification'].get(k) for k in PHYSICAL_KEYS):
        raise ValueError('gas control material/spectral/boundary settings differ')
    if opacity is None:opacity=json.loads(read_text(original(work,'opacity-control.json')))
    op=opacity['provenance']
    if (opacity['sha256']!=provenance['opacity_sha256'] or op['mode']!='gas'
            or op['prepared']!=provenance['prepared']
            or any(op[k]!=provenance[k] for k in ['XH','X3'])
            or any(op['specification'].get(k)!=spec.get(k) for k in PHYSICAL_KEYS)):
        raise ValueError('gas opacity source physics differs')
    cfg=read_text(original(work,'ember-condensates.cfg')).splitlines()
    if 'gas' not in cfg or 'equilibrium' in cfg:
        raise ValueError('gas source did not disable depletion')
    receipt=json.loads(read_text(original(work,'completed.json')))
    if input_fingerprint(provenance['prepared']['executables']['tlusty'],work,archived=True)!=receipt['input_sha256']:
        raise ValueError('gas input fingerprint differs')
    for name, expected in receipt['outputs'].items():
        path=original(work,name);data=path.read_bytes()
        if path.suffix=='.gz':data=gzip.decompress(data)
        if hashlib.sha256(data).hexdigest()!=expected:
            raise ValueError('gas original output differs')
    if read_text(original(work,'opacity.sha256')).strip()!=provenance['opacity_sha256']:
        raise ValueError('gas opacity identity differs')
    log=read_text(original(work,'run.log'))
    names={'atmosphere_input':'fort.5','parameters':'tas',
           'element_masses':'ember-masses.dat','initial_structure':'fort.8'}
    coords=[provenance[k] for k in ['XH','X3','teff_K','log_g']]
    source_inputs({k:read_text(original(work,v)) for k,v in names.items()},spec,*coords,log)
    result=source_state(log,read_text(original(work,'fort.9')),*coords[2:],
                        {'temperature_K':temperatures(spec),'density_g_cm3':sequence(spec['log_density'])},spec['tau'])
    saved=json.loads(read_text(original(work,'validated.json')))
    if saved['provenance']!=provenance or saved['diagnostics']!=result:
        raise ValueError('gas validation receipt differs')
    return result


def compare(work, gas, archive, output):
    result=validate(work)
    gas=Path(gas);table=(gas/'opacity.bin').resolve(strict=True)
    spec=result['provenance']['specification']
    opacity={'provenance':json.loads((table.parent/'provenance.json').read_text()),'sha256':digest(table)}
    abundance,_=composition(*result['coordinates'][:2],spec['metals'])
    validate_table(table,abundance,temperatures(spec),sequence(spec['log_density']))
    control=gas_state(gas,result,opacity)
    result,hashes=archive_model(work,archive)
    archive=Path(archive);gas=Path(gas)
    for name in [n for n in MODEL_FILES if not n.startswith('chemistry-')]:
        source=original(gas,name);target=archive/'gas-control'/(name+'.gz')
        data=source.read_bytes() if source.suffix=='.gz' else gzip.compress(source.read_bytes(),mtime=0)
        target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists() and target.read_bytes()!=data:raise ValueError('gas archive differs')
        target.write_bytes(data);hashes[str(target.relative_to(archive))]=digest(target)
    target=archive/'gas-control/opacity-control.json.gz'
    data=gzip.compress((json.dumps(opacity,indent=2)+'\n').encode(),mtime=0)
    if target.exists() and target.read_bytes()!=data:raise ValueError('gas opacity archive differs')
    target.write_bytes(data);hashes[str(target.relative_to(archive))]=digest(target)
    if gas_state(archive/'gas-control',result)!=control:
        raise ValueError('gas archive cannot reproduce control')
    report={'status':'One matched-grid depletion/gas comparison; not a stellar lifetime or full atmosphere family.',
            'coordinates':result['coordinates'],'archive':str(archive),
            'diagnostics':result['diagnostics'],'chemistry':result['chemistry'],
            'same_material_grid_gas_control':control,
            'isolated_condensation_relative_boundary_change':{
                k:result['diagnostics'][k]/control[k]-1 for k in ['T','Pgas','source_density']},
            'comparison_note':'Identical source version, composition, material/spectral grid and atmosphere settings; equilibrium depletion enabled versus disabled. Grain opacity and enthalpy omitted.',
            'files_sha256':hashes}
    Path(output).write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['work','gas','archive','output']:p.add_argument(name,type=Path)
    a=p.parse_args()
    print(json.dumps(compare(a.work,a.gas,a.archive,a.output)['isolated_condensation_relative_boundary_change'],indent=2))
