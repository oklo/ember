#!/usr/bin/env python3
"""Make the warmer OP pilot's missing envelope density coverage explicit."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from fetch_tops_composition import digest
from op_native_spectra_v2 import read_spectra
from reduce_tops_group_factors import verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    pilot_path=Path('docs/results/op_hot_pilot_v1.json');pilot=json.loads(pilot_path.read_text())
    verify(pilot['input_sha256']);verify(pilot['output_sha256'])
    warm_path=Path('docs/results/op_warm_envelope_v1.json');composition=json.loads(warm_path.read_text())['composition']
    profile_path=Path('docs/results/opacity_connections_v1.json');saved=json.loads(profile_path.read_text())
    rows=sorted(saved['rows'],key=lambda r:r['query'][4])
    locus=PchipInterpolator(np.log10([r['query'][4] for r in rows]),np.log10([r['query'][5] for r in rows]),extrapolate=False)
    config_path=Path('/tmp/ember-op-hot-pilot-v1/configuration.json');config=json.loads(config_path.read_text())
    records=[]
    for entry in pilot['source_records']:
        source=json.loads(Path(entry['path']).read_text())
        if (source['X'],source['Z'])!=(composition['atomic_X'],composition['atomic_Z']):continue
        it=source['temperature_index'];logt=.025*it
        logr=float(locus(logt))+math.log10(composition['atomic_mass_scale'])-3*logt+18
        states=source['native_states'];last=states[-1]
        native=read_spectra(config['source_paths'][str(it)]['1'])
        at_source_end=last['electron_index']==native['states'][-1]['electron_index']
        covered=states[0]['logR']<=logr<=last['logR']
        records.append(dict(temperature_index=it,temperature_K=10**logt,target_logR=logr,
            source_logR_min=states[0]['logR'],source_logR_max=last['logR'],covered=covered,
            last_mixture_electron_index=last['electron_index'],last_source_electron_index=native['states'][-1]['electron_index'],
            reaches_native_source_density_end=at_source_end,source_masks=source['source_masks'],
            failed_native_states=len(source['failures']),
            exclusion_reason=None if covered else 'Saved-envelope density exceeds the available OP native electron-density grid.'))
    records.sort(key=lambda r:r['temperature_index'])
    if any(not r['covered'] and not r['reaches_native_source_density_end'] for r in records):
        raise ValueError('unexplained missing pilot density support')
    result=dict(scope=__doc__,outcome='incomplete_native_density_coverage',accepted_for_stellar_opacity=False,
        sampled_current_composition_isotherms=len(records),supported_isotherms=sum(r['covered'] for r in records),
        unsupported_isotherms=sum(not r['covered'] for r in records),records=records,
        new_source_or_mixture_calculations=0,
        conclusion='The hotter OP spectra do not supply a continuous source connection along the current envelope; expanding composition nodes cannot supply missing density data.',
        input_sha256={str(q.resolve()):digest(q) for q in (Path(__file__),pilot_path,warm_path,profile_path,config_path,Path('scripts/op_native_spectra_v2.py'))})
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('outcome','sampled_current_composition_isotherms','supported_isotherms','unsupported_isotherms')}),flush=True)


if __name__=='__main__':main()
