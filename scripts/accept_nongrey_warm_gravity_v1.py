#!/usr/bin/env python3
"""Verify and install the completed warm atmosphere gravity extension.

The 48 added sources cover X=0.15/0.2, helium-3 fractions 0/0.12,
5400..6400 K, and log g=5.7/6.0. Existing values and missing-state masks
remain exact. This accepts a gas boundary in its supported domain, not a
cold white-dwarf atmosphere or a new stellar evolutionary trajectory.
"""
import json
from datetime import datetime, timezone
from pathlib import Path
from prepare_nongrey_sources import digest


def main():
    output=Path('docs/results/nongrey_warm_gravity_acceptance_v1.json')
    if output.exists():raise FileExistsError('preserve completed acceptance')
    inputs={str(Path(__file__).resolve()):digest(__file__)};reports={}
    def load(path):
        path=Path(path);r=json.loads(path.read_text());reports[str(path)]=digest(path)
        for name in ('input_files_sha256','input_sha256'):
            for key,value in r.get(name,{}).items():
                key=str(Path(key).resolve())
                if key in inputs and inputs[key]!=value:raise ValueError('inconsistent input identity')
                inputs[key]=value
        return r
    previous=load('docs/results/nongrey_t6400_acceptance_v1.json')
    eos_acceptance=load('docs/results/numerical_electron_hot_dense_exhaustion_acceptance_v1.json')
    if not previous['accepted_for_gas_trajectory'] or not eos_acceptance['accepted_for_gas_trajectory']:
        raise ValueError('previous atmosphere or EOS is unaccepted')
    runtime=load('docs/results/nongrey_warm_gravity_interpolation_v1.json')
    if runtime['source_models']!=300 or runtime['runtime_supported_source_nodes']!=300 or runtime['unsupported_EOS_nodes']:
        raise ValueError('incomplete atmosphere/EOS support')
    if runtime['script_sha256']!=digest('scripts/audit_nongrey_extension.py'):
        raise ValueError('runtime audit source changed')
    if len(runtime['heldout_comparisons'])!=6:raise ValueError('missing independent atmospheres')
    maximum=max(abs(v) for r in runtime['heldout_comparisons'] for v in r['relative_difference'].values())
    if maximum>.01 or max(runtime['retained_max_relative_difference'].values())>1e-12:
        raise ValueError('independent interpolation or retention check fails')
    retention=load('docs/results/nongrey_warm_gravity_retained_v1.json')
    if not retention['passed'] or retention['retained_source_states']!=252 or retention['candidate_source_states']!=300:
        raise ValueError('old atmosphere states were not retained')
    planpath=Path('data/atmosphere/sources/nongrey_warm_gravity_assembly_v1_specification.json');plan=load(planpath)
    if digest(plan['previous_acceptance'])!=plan['previous_acceptance_sha256']:
        raise ValueError('previous acceptance changed')
    manifestpath=Path('/tmp/ember-nongrey-warm-gravity-candidate-v1.manifest.json');manifest=load(manifestpath)
    candidate=Path('/tmp/ember-nongrey-warm-gravity-candidate-v1.dat')
    if digest(candidate)!=manifest['table_sha256'] or manifest['accepted_states']!=300:
        raise ValueError('candidate differs from source assembly')
    depths=[]
    paths=[r['report'] for r in plan['completed_repair_depth_checks']]
    paths += [str(p) for p in sorted(Path('docs/results').glob('nongrey_*_warm_column_v1.json'))]
    paths += ['docs/results/nongrey_x150_y000_t6400_g570_warm_column_v2.json',
              'docs/results/nongrey_x150_y000_t6000_g570_depth400_v1.json',
              'docs/results/nongrey_x200_y000_t6000_g570_depth400_v1.json',
              'docs/results/nongrey_x150_y000_t6000_g600_column_v1.json']
    if len(paths)!=17 or len(set(paths))!=17:raise ValueError('unexpected depth comparison set')
    for path in paths:
        r=load(path)
        if not r['passed'] or r['relative_matching_state_tolerance']!=1e-4 or r['maximum_matching_state_relative_difference']>1e-4:
            raise ValueError('depth comparison failed')
        depths.append(r)
    chemistry=[]
    for g in ('570','600'):
        r=load(f'docs/results/nongrey_warm_gravity_condensation_g{g}_v1.json')
        for name,value in r['helper_sha256'].items():inputs[str((Path('scripts')/name).resolve())]=value
        for row in r['records']:
            if row['layers_with_condensates']!=0:raise ValueError('gas-only atmosphere has equilibrium condensates')
            for field in ('atmosphere_input_sha256','diagnostic_files_sha256'):
                for name,value in row[field].items():inputs[str(Path(name).resolve())]=value
            chemistry.append(row)
    if len(chemistry)!=48 or {tuple(r['coordinates']) for r in chemistry}!={tuple(r['coordinates']) for r in plan['new_sources']}:
        raise ValueError('condensation checks do not cover all new structures')
    for path in manifest['capacity_replays']:load(path)
    for name,h in {**inputs,**reports}.items():
        if digest(name)!=h:raise ValueError('acceptance input changed: '+name)
    eos=Path(eos_acceptance['eos_family'])
    installed=Path('data/atmosphere/nongrey_gs98_z020_exhaustion_t6400_g600_v1.dat')
    with installed.open('xb') as stream:stream.write(candidate.read_bytes())
    if digest(installed)!=manifest['table_sha256']:raise ValueError('installed bytes differ')
    result={'scope':__doc__,'checked_utc':datetime.now(timezone.utc).isoformat(),
            'accepted_for_gas_trajectory':True,'table':str(installed),'table_sha256':digest(installed),
            'assembly_manifest':str(manifestpath),'assembly_manifest_sha256':digest(manifestpath),
            'assembly_specification':str(planpath),'assembly_specification_sha256':digest(planpath),
            'EOS':str(eos),'EOS_sha256':digest(eos),'previous_acceptance':plan['previous_acceptance'],
            'source_models':300,'retained_source_models':252,'retained_missing_masks':retention['retained_missing_states'],
            'source_values_and_masks_exact':True,'independent_matching_state_relative_criterion':.01,
            'maximum_independent_matching_state_relative_difference':maximum,
            'independent_comparisons':[{'coordinates':r['coordinates'],'relative_difference':r['relative_difference']} for r in runtime['heldout_comparisons']],
            'depth_comparisons':len(depths),'maximum_depth_matching_state_relative_difference':max(r['maximum_matching_state_relative_difference'] for r in depths),
            'condensation_profiles':len(chemistry),'profiles_with_condensates':0,'runtime_supported_source_nodes':300,
            'new_source_hydrogen':[.15,.2],'new_source_helium3':[0,.12],
            'new_source_temperature_range_K':[5400,6400],'new_source_logg':[5.7,6.0],
            'capacity_replays':manifest['capacity_replays'],'reports_sha256':reports,
            'input_sha256':inputs,'installation_script_sha256':digest(__file__)}
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'accepted':True,'table':str(installed),'source_models':300,'maximum_independent_relative_difference':maximum,'depth_checks':len(depths)}))


if __name__=='__main__':main()
