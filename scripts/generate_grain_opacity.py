#!/usr/bin/env python3
"""Generate a separate grain spectrum from an explicit material/size recipe.

Writes wavelength, absorption, scattering and scattering asymmetry. Opacities
are per gram of that grain material, not per gram of stellar atmosphere.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from grain_opacity import prepared, spectrum
from prepare_nongrey_sources import digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path); p.add_argument('specification', type=Path)
    p.add_argument('work', type=Path)
    a = p.parse_args(); receipt = prepared(a.source); spec = json.loads(a.specification.read_text())
    if spec.get('format') != 1 or spec.get('shape') != 'independent compact spheres':
        raise ValueError('unsupported grain specification')
    if not isinstance(spec['wavelength_count'], int) or not 2 <= spec['wavelength_count'] <= 100000:
        raise ValueError('invalid wavelength count')
    low, high = spec['wavelength_range_um']
    if not 0 < low < high:
        raise ValueError('invalid wavelength interval')
    for name in ['density_reference', 'material_assumption', 'size_distribution_assumption']:
        if not isinstance(spec.get(name), str) or not spec[name].strip():
            raise ValueError('explicit physical assumption/reference required: ' + name)
    a.work.mkdir(parents=True, exist_ok=False)
    wave = np.geomspace(low, high, spec['wavelength_count'])
    result = spectrum(receipt, spec['optical_data_file'], spec['density_g_cm3'],
                      spec['radii_um'], spec['number_weights'], wave)
    target = a.work / 'grain-opacity.csv'
    np.savetxt(target, np.column_stack(list(result.values())), delimiter=',',
               header=','.join(result), comments='', fmt='%.17g')
    record = {'format': 1, 'specification': spec, 'specification_sha256': digest(a.specification),
        'source_receipt_sha256': digest(a.source), 'probe_sha256': receipt['probe_sha256'],
        'optical_source': receipt['optical_data'][spec['optical_data_file']],
        'opacity': target.name, 'opacity_sha256': digest(target),
        'scope': 'Absorption and scattering per dust mass. No chemistry, settling, or atmosphere coupling is implied.',
        'interpolation': 'Linear n and logarithmic positive k in log wavelength; linear k on intervals touching zero; no extrapolation'}
    (a.work / 'receipt.json').write_text(json.dumps(record, indent=2) + '\n')
    print(a.work / 'receipt.json')


if __name__ == '__main__':
    main()
