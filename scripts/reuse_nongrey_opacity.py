"""Reuse a retained opacity table only with its accepted source provenance."""
import json
from pathlib import Path

from prepare_nongrey_sources import digest


def reuse_accepted_table(donor, destination, prepared, spec, hydrogen, helium3):
    """Validate an existing single-composition table, then link it without synthesis.

    Effective temperature, gravity and atmosphere resolution may change. The
    opacity composition, material data, sampling and table axes must not change.
    Missing or incompatible evidence is an error, never a synthesis fallback.
    """
    # Import lazily: the source reader also imports the atmosphere generator.
    from assemble_helium_fraction_atmospheres import load_plane
    from generate_nongrey_grid import temperatures, sequence, composition
    from nongrey_opacity import validate_table

    donor=Path(donor).resolve(strict=True);destination=Path(destination)
    old=json.loads((donor/'specification.json').read_text())
    provenance=json.loads((donor/'provenance.json').read_text())
    fields=['hydrogen','helium3','metals','log_temperature','log_density',
            'wavelength_A','line_threshold','microturbulence_km_s',
            'opacity_method','opacity_frequencies','synthesis_spacing_A']
    for key in fields:
        if old[key]!=spec[key]:
            raise ValueError('retained opacity setting differs: '+key)
    if temperatures(old)!=temperatures(spec):
        raise ValueError('retained opacity temperature values differ')
    if spec['hydrogen']!=[hydrogen] or spec['helium3']!=[helium3]:
        raise ValueError('retained opacity requires the exact single composition')
    for key in ['data_sha256','line_list_sha256','opacity_method']:
        if provenance[key]!=prepared[key]:
            raise ValueError('retained opacity source identity differs: '+key)
    if provenance['executables']['synspec']!=prepared['executables']['synspec']:
        raise ValueError('retained opacity synthesis executable differs')
    fraction=helium3/(1-hydrogen-sum(spec['metals']))
    _,rows,_,files,_=load_plane(donor,hydrogen,fraction)
    if not rows:
        raise ValueError('retained opacity has no independently verifiable accepted source')
    source=donor/'plane-000/opacity/fort.63'
    abundance,_=composition(hydrogen,helium3,spec['metals'])
    table=validate_table(source,abundance,temperatures(spec),sequence(spec['log_density']))
    checksum=digest(source)
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.symlink_to(source)
    if digest(destination)!=checksum:
        raise ValueError('retained opacity link checksum differs')
    receipt=dict(donor=str(donor),table=str(source),table_sha256=checksum,
                 target=str(destination),shape=table['shape'],
                 accepted_source_columns=len(rows),input_sha256=files,
                 helper_sha256=digest(__file__),
                 synthesis_performed=False,opacity_settings_compared=fields)
    destination.with_name('merged-reuse.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt
