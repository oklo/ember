#!/usr/bin/env python3
"""Measure C-to-N composition effects using the retained atmosphere source.

The comparison holds final hydrogen, helium-3, temperature and density fixed.
It contrasts the fixed GS98 lookup with a physical C12+2H->N14 inventory.
Only the lookup replaces captured hydrogen by helium-4; physical helium-4
is lower by the captured-proton mass. No stellar state is modified.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path
import shutil

import numpy as np
import generate_nongrey_grid as gen
from nongrey_opacity import validate_table
from prepare_nongrey_sources import digest,data_digest


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('output',type=Path)
    ap.add_argument('--scratch',type=Path,required=True)
    args=ap.parse_args();assert not args.output.exists() and not args.scratch.exists()
    root=args.scratch.resolve();root.mkdir()
    prepared_path=Path('/tmp/ember-primary-nongrey-restore-v1/prepared.json')
    base=Path('/tmp/ember-nongrey-x020-opacity-13k-v1')
    prepared=json.loads(prepared_path.read_text());spec=json.loads((base/'specification.json').read_text())
    assert prepared['opacity_method']==spec['opacity_method']=='sampling'
    assert spec['hydrogen']==[.2] and spec['helium3'][0]==0.
    for name in ('synspec','tlusty'):assert digest(prepared[name])==prepared['executables'][name]
    for path,h in zip(prepared['line_lists'],prepared['line_list_sha256'],strict=True):assert digest(path)==h
    assert data_digest(Path(prepared['synple'])/'data')==prepared['data_sha256']
    identities={str(p.resolve()):digest(p) for p in [Path(__file__),Path(gen.__file__),prepared_path,base/'specification.json']}
    original=gen.composition;physical_compositions={};tasks=[]
    temperatures=gen.temperatures(spec);densities=gen.sequence(spec['log_density'])
    for fraction,label in [(0.,'GS98'),(.5,'half_carbon'),(1.,'converted_carbon')]:
        def composition(x,y,metals,f=fraction):
            assert x==.2 and y==0.
            abundance,weights=original(x,y,metals)
            original_carbon=abundance[5]
            converted=original_carbon*f
            # N_i/N_H: remove C, add the same number of N nuclei, remove
            # half as many He4 nuclei from the fixed-Z lookup composition.
            abundance[5]=max(original_carbon*(1-f),1e-99)
            abundance[6]+=converted
            abundance[1]-=.5*converted
            assert min(abundance)>0
            mass=sum(a*w for a,w in zip(abundance,weights))*x*gen.SOURCE_HMASS/gen.MU
            assert abs(mass-1)<1e-12
            old,_=original(x,y,metals)
            # Fully ionized electron count is identical between this physical
            # inventory and the proxy at the same final XH.
            electron_difference=sum((a-b)*(i+1)*x for i,(a,b) in enumerate(zip(abundance,old)))
            assert abs(electron_difference)<1e-14
            physical_compositions[label]=dict(converted_carbon_fraction=f,
                hydrogen=x,helium3=y,physical_Z=sum(metals)+2*converted*x,
                physical_helium4=1-x-sum(metals)-2*converted*x,
                captured_proton_baryon_fraction=2*converted*x,
                fully_ionized_electron_moment_difference=electron_difference,
                ion_moment_difference=sum((a-b)*x for a,b in zip(abundance,old)),
                number_abundances_relative_to_H=abundance,source_mass_weights=weights)
            return abundance,weights
        gen.composition=composition
        for index in [3,5,7,9]:
            directory=root/label/f'temperature-{index:03d}';T=temperatures[index]
            gen.opacity_inputs(directory,prepared,spec,.2,0.,T)
            abundance,_=composition(.2,0.,spec['metals'])
            if fraction==0:
                old=base/'plane-000/opacity'/f'temperature-{index:03d}'
                receipt=json.loads((old/'completed.json').read_text())
                assert receipt['input_sha256']==gen.input_fingerprint(prepared['executables']['synspec'],directory)
                for filename in ['fort.63','fort.29','run.log','completed.json']:
                    if filename!='completed.json':assert digest(old/filename)==receipt['outputs'][filename]
                    shutil.copy2(old/filename,directory/filename)
                    identities[str((old/filename).resolve())]=digest(old/filename)
            tasks.append((label,T,directory,abundance))
    gen.composition=original
    inputs=dict(created_utc=datetime.now(timezone.utc).isoformat(),prepared=prepared,specification=spec,
                physical_compositions=physical_compositions,input_sha256=identities,
                task_count=len(tasks),new_source_isotherms=8,reused_source_isotherms=4)
    (root/'inputs.json').write_text(json.dumps(inputs,indent=2)+'\n')

    def calculate(task):
        label,T,directory,abundance=task
        gen.execute(prepared['synspec'],directory,['fort.63','fort.29'])
        table=validate_table(directory/'fort.63',abundance,[T],densities)
        nu=np.asarray(table['frequency'])[::-1]
        logk=np.asarray(table['log_opacity'],dtype=float).reshape(table['shape'])[::-1,:,0]
        x=6.6256e-27*nu/(1.38054e-16*T)
        w=x**4*np.exp(-x)/(-np.expm1(-x))**2
        norm=np.trapezoid(w,x)
        means=norm/np.trapezoid(w[:,None]*np.exp(-logk),x,axis=0)
        assert np.isfinite(means).all() and np.min(means)>0
        print(json.dumps(dict(case=label,temperature_K=T,completed=True)),flush=True)
        return dict(case=label,temperature_K=T,density_g_cm3=densities,absorption_mean_cm2_g=means.tolist(),
                    rosseland_weight_fraction=float(norm/(4*math.pi**4/15)),
                    source_seconds=json.loads((directory/'completed.json').read_text())['seconds'],
                    table_sha256=digest(directory/'fort.63'),table_path=str(directory/'fort.63'))
    with ThreadPoolExecutor(max_workers=2) as pool:records=list(pool.map(calculate,tasks))
    for row in records:
        baseline=next(r for r in records if r['case']=='GS98' and r['temperature_K']==row['temperature_K'])
        row['ratio_to_GS98']=(np.array(row['absorption_mean_cm2_g'])/baseline['absorption_mean_cm2_g']).tolist()
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed',
        accepted_for_stellar_evolution=False,physical_compositions=physical_compositions,records=records,
        source_job_count=8,reused_source_isotherms=4,input_sha256=identities,
        source_mean='Rosseland absorption mean over the source frequency interval; scattering and grains excluded.',
        limitations=['This is a fixed material-state sensitivity, not a coupled atmosphere or stellar response.',
            'The samples use XH=.20 and X3=0; retained restart models have nonzero X3.',
            'The fully ionized particle-count identity does not bound cool molecular-opacity errors.',
            'No historical carbon conversion, mixing or nuclear fuel is changed by this calculation.'])
    report['input_sha256'][str(root/'inputs.json')]=digest(root/'inputs.json')
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(outcome='passed',source_job_count=8)),flush=True)


if __name__=='__main__':
    main()
