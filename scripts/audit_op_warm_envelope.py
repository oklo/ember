#!/usr/bin/env python3
"""Evaluate local OP spectra across the saved star's warm envelope.

Reuse the exact common envelope composition and previously calculated native
states. Preserve invalid states and coverage failures. This comparison does
not supply the four absent trace elements or select a stellar opacity.
"""
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import tarfile

import numpy as np

from audit_mesa_oplib_native import read_native, evaluate
from compare_mesa_oplib_means_v2 import stencil
from compare_op_native_mixtures import select_connected
from op_native_spectra_v2 import read_spectra, read_mesh, native_rosseland
from audit_op_mixture_corrections import native_task
from audit_op_plasma_transport import native_calculation, case_interpolations

ROOT = Path(__file__).resolve().parents[1]
WORK = Path('/tmp/ember-op-warm-envelope-v1')
OLD = Path('/tmp/ember-op-selected-spectra-v1/OP4STARS_1.3/mono')
ARCHIVE = Path('/tmp/ember-op-monochromatic-fetch-v2/OP4STARS_1.3.tar.xz')
MEAN_ARCHIVE = Path('/tmp/ember-mesa-native-input-fetch-v1/kap_input_data.tar.xz')
OUTPUT = ROOT/'docs/results/op_warm_envelope_v1.json'
DOMAIN = ROOT/'docs/results/op_native_stellar_domain_v2.json'
BASE = ROOT/'docs/results/op_native_mixture_comparison_v1.json'
CORRECTIONS = ROOT/'docs/results/op_mixture_corrections_v1.json'
PLASMA = ROOT/'docs/results/op_plasma_transport_v1.json'
MU = 1.66053906660e-24
A0 = 5.29177210903e-9


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def extract_requested(archive, expected_hash, requests):
    """Extract only explicitly named regular files, keeping existing files."""
    if digest(archive) != expected_hash:
        raise ValueError('archive identity mismatch')
    pending = {name:path for name,path in requests.items() if not path.exists()}
    extracted = []
    with tarfile.open(archive,'r:xz') as src:
        for member in src:
            if member.name not in pending:
                continue
            if not member.isfile():
                raise ValueError('requested member is not a regular file')
            path = pending.pop(member.name)
            body = src.extractfile(member).read()
            if len(body) != member.size:
                raise ValueError('truncated archive member')
            path.parent.mkdir(parents=True,exist_ok=True)
            with path.open('xb') as f:
                f.write(body)
            extracted.append(dict(archive_path=member.name,path=str(path),bytes=len(body),sha256=digest(path)))
    if pending:
        raise ValueError('missing archive members: '+str(list(pending)))
    return extracted


def calculate(task):
    key, native = task
    row, arrays = native_task(native)
    try:
        if not row['used_in_mean_interpolation']:
            raise ValueError('nonpositive corrected spectrum')
        plasma = native_calculation((row,arrays['total'],arrays['free_electron_scattering']))
    except ValueError as error:
        return key, row, None, str(error), arrays
    return key, row, plasma, None, arrays


def main():
    if OUTPUT.exists() or WORK.exists():
        raise FileExistsError('completed or partial warm-envelope work already exists')
    WORK.mkdir()
    inputs = {str(p):digest(p) for p in [Path(__file__),DOMAIN,BASE,CORRECTIONS,PLASMA,
        ROOT/'include/ember/gs98_mixture.hpp',ROOT/'docs/reports/2026-09-11/evolution_latest_profile.csv',
        ROOT/'scripts/audit_mesa_oplib_native.py',ROOT/'scripts/compare_mesa_oplib_means_v2.py',
        ROOT/'scripts/compare_op_native_mixtures.py',ROOT/'scripts/op_native_spectra_v2.py',
        ROOT/'scripts/op_mixture_corrections.py',ROOT/'scripts/op_plasma_transport.py',
        ROOT/'scripts/audit_op_mixture_corrections.py',ROOT/'scripts/audit_op_plasma_transport.py',
        ROOT/'scripts/audit_tops_electron_dispersion.py',OLD/'m01.mesh']}
    domain = json.loads(DOMAIN.read_text())
    stars = domain['stellar_records'][482:505]
    assert len(stars) == 23 and all(r['four_temperature']['supported'] for r in stars)
    composition = {(r['X'],r['Y3'],r['Y4'],r['Z']) for r in stars}
    assert len(composition) == 1
    x,y3,y4,z = next(iter(composition))
    base = json.loads(BASE.read_text())
    old_case = next(c for c in base['cases'] if c['label']=='stellar_zone493')
    assert old_case['temperature_K'] == stars[493-482]['temperature_K']
    weights = {int(k):v for k,v in old_case['included_atom_weights'].items()}
    assert weights[1] == x and weights[2] == y3/3+y4/4
    for p,h in base['input_sha256'].items():
        if digest(p) != h:
            raise ValueError('base mixture input changed: '+p)
        inputs[p] = h
    temperatures = sorted(set(it for r in stars for it in r['four_temperature']['temperature_indices']))
    paths = {}; requested = {}
    for it in temperatures:
        for element in weights:
            name = f'm{element:02d}.{it:03d}'
            path = OLD/name if (OLD/name).exists() else WORK/'spectra'/name
            paths[element,it] = path
            if not path.exists():
                requested['OP4STARS_1.3/mono/'+name] = path
    extraction = extract_requested(ARCHIVE,'aeae2b31e62c7cebc100be2813e9b976de0681a31faa5fa2716c766cc1a6e809',requested)
    inputs[str(ARCHIVE)] = 'aeae2b31e62c7cebc100be2813e9b976de0681a31faa5fa2716c766cc1a6e809'
    (WORK/'extraction.json').write_text(json.dumps(extraction,indent=2)+'\n')
    dv,u = read_mesh(OLD/'m01.mesh')
    planes = {}; reader_checks = []; masked = []
    for it in temperatures:
        elements = {}
        for element in weights:
            path = paths[element,it]; inputs[str(path)] = digest(path)
            data = read_spectra(path); elements[element] = data
            errors = []
            for state in data['states']:
                if state['positive_cross_section']:
                    mean = float(native_rosseland(state['cross_section_atomic'],dv))
                    error = mean/state['rosseland_atomic']-1
                    if abs(error) > data['packing_tolerance']+.00011:
                        raise ValueError('decoded native mean disagrees with source')
                    errors.append(error)
                else:
                    masked.append(dict(element=element,temperature_index=it,electron_index=state['electron_index']))
            reader_checks.append(dict(element=element,temperature_index=it,positive_states=len(errors),
                                      maximum_relative_mean_difference=max(map(abs,errors),default=None)))
        jns = np.array([s['electron_index'] for s in elements[1]['states']])
        sigma = np.zeros((len(jns),len(u))); charge = np.zeros(len(jns)); ions = np.zeros((len(jns),28))
        valid = np.ones(len(jns),dtype=bool)
        for element,w in weights.items():
            states = elements[element]['states']
            assert [s['electron_index'] for s in states] == jns.tolist()
            for j,state in enumerate(states):
                sigma[j] += w*state['cross_section_atomic']
                charge[j] += w*state['electrons_per_atom']
                valid[j] &= state['positive_cross_section']
                for ion,fraction in zip(state['ion_indices'],state['ion_fractions']):
                    q = element-1-int(ion)
                    if q > 0:
                        ions[j,q-1] += w*fraction
        assert np.all(charge > 0)
        density = MU*10**(.25*jns)/charge
        assert np.all(np.diff(density) > 0)
        planes[it] = dict(jns=jns,sigma=sigma,charge=charge,ions=ions,valid=valid,density=density)
    cases = []; needed = set(); coverage_failures = []
    for star in stars:
        case = dict(label='stellar_envelope',zone=star['zone'],temperature_K=star['temperature_K'],
                    density_g_cm3=star['density_baryonic_g_cm3'],planes=[])
        for it in star['four_temperature']['temperature_indices']:
            p = planes[it]; density_rows = []
            for number in (2,4):
                try:
                    selected = select_connected(np.log(p['density']),p['valid'],math.log(case['density_g_cm3']),number)
                except ValueError as error:
                    coverage_failures.append(dict(zone=star['zone'],temperature_index=it,density_points=number,reason=str(error)))
                    continue
                chosen = p['jns'][selected].tolist()
                needed.update((it,jn) for jn in chosen)
                density_rows.append(dict(density_points=number,electron_indices=chosen,native_densities=p['density'][selected].tolist()))
            case['planes'].append(dict(temperature_index=it,temperature_K=10**(.025*it),density_interpolations=density_rows))
        cases.append(case)
    old_corrected = json.loads(CORRECTIONS.read_text()); old_plasma = json.loads(PLASMA.read_text())
    for previous in (old_corrected,old_plasma):
        for p,h in previous['input_sha256'].items():
            if digest(p) != h:
                raise ValueError('reused correction input changed: '+p)
            inputs[p] = h
    old_rows = {(r['temperature_index'],r['electron_index']):r for r in old_corrected['native_states'] if r['label']=='stellar_zone493'}
    old_results = {(r['temperature_index'],r['electron_index']):r for r in old_plasma['native_states'] if r['label']=='stellar_zone493'}
    artifact = old_corrected['spectral_artifact']; assert digest(artifact['path']) == artifact['sha256']
    inputs[artifact['path']] = artifact['sha256']
    tasks = []; records = {}; result_rows = {}; arrays = {}; reused = []; failures = []
    for key in sorted(needed):
        it,jn = key
        if key in old_results:
            # Every included atom weight is identical to the saved calculation.
            # Reuse its correction and plasma integration as well as its mean.
            records[key] = {**old_rows[key],'label':'stellar_envelope'}
            result_rows[key] = {**old_results[key],'label':'stellar_envelope'}
            reused.append(list(key)); continue
        p = planes[it]; j = p['jns'].tolist().index(jn)
        tasks.append((key,('stellar_envelope',it,jn,p['sigma'][j],p['ions'][j],p['charge'][j],sum(weights.values()))))
    (WORK/'task_summary.json').write_text(json.dumps(dict(cases=len(cases),needed=len(needed),reused=reused,new=len(tasks)))+'\n')
    print(f'envelope: {len(cases)} layers, {len(tasks)} new native states, {len(reused)} reused',flush=True)
    with ProcessPoolExecutor(max_workers=4) as pool:
        for key,row,result,error,asset in pool.map(calculate,tasks):
            records[key] = row
            for name,value in asset.items():
                arrays[f't{key[0]}_n{key[1]}_{name}'] = value
            if error:
                failures.append(dict(temperature_index=key[0],electron_index=key[1],reason=error))
            else:
                result_rows[key] = result
    spectrum_file = WORK/'new_corrected_spectra.npz'
    np.savez_compressed(spectrum_file,**arrays)
    # MESA's original mean tables around the isotope-mapped X and Z.
    import re
    metals = [tuple(map(float,r.split(','))) for r in re.findall(r'^\s*\{([^{}]+)\}, //',(ROOT/'include/ember/gs98_mixture.hpp').read_text(),re.M)]
    metal_scale = sum(r[3]*r[2]/r[1] for r in metals)
    scale = 1.00794*x+4.002602*(y3/3+y4/4)+z*metal_scale
    xa = 1.00794*x/scale; za = z*metal_scale/scale
    xgrid = [.1,.15,.2,.25]; zgrid = [.018,.019,.02,.021]
    mean_requests = {}
    for xx in xgrid:
        for zz in zgrid:
            name = f'oplib_gs98_z{zz}_x{xx}.data'
            mean_requests['kap_input_data/oplib/OPLIB_GS98_LOG_1194_25E/'+name] = WORK/'means'/name
    mean_extraction = extract_requested(MEAN_ARCHIVE,'042eb369c8f15c808758bee7938ed138f686a7d31bb48475bbb719be24cec997',mean_requests)
    inputs[str(MEAN_ARCHIVE)] = '042eb369c8f15c808758bee7938ed138f686a7d31bb48475bbb719be24cec997'
    parsed_tables = [read_native(e) for e in mean_extraction]
    tables = {(t['X'],t['Z']):t for t in parsed_tables}
    for e in mean_extraction:
        inputs[e['path']] = e['sha256']
    samples = []; unavailable = []
    for case in cases:
        missing = [(p['temperature_index'],jn) for p in case['planes'] for d in p['density_interpolations'] for jn in d['electron_indices'] if (p['temperature_index'],jn) not in result_rows]
        if missing or any(len(p['density_interpolations']) != 2 for p in case['planes']):
            unavailable.append(dict(zone=case['zone'],missing_corrected_states=missing)); continue
        lookup = {('stellar_envelope',it,jn):r for (it,jn),r in result_rows.items()}
        row = case_interpolations(case,lookup); row['zone'] = case['zone']
        row['logR_atomic'] = math.log10(scale*case['density_g_cm3'])-3*math.log10(case['temperature_K'])+18
        comparisons = []
        for order in (2,4):
            try:
                ii,wx = stencil(xgrid,xa,order); jj,wz = stencil(zgrid,za,order)
                logk = math.fsum(a*b*evaluate(tables[xgrid[i],zgrid[j]],math.log10(case['temperature_K']),row['logR_atomic'],4)
                                for i,a in zip(ii,wx) for j,b in zip(jj,wz))
                mean = scale*10**logk
                op = next(v for v in row['interpolations'] if v['temperature_points']==4 and v['density_points']==4)
                comparisons.append(dict(composition_points=order,mesa_baryonic_cm2_g=mean,
                                        op_relative_to_mesa=op['uncut_finite_mean']/mean-1))
            except ValueError as error:
                row['mesa_unavailable_reason'] = str(error)
        row['mesa_comparisons'] = comparisons
        samples.append(row)
    reference = next(r for r in old_plasma['cases'] if r['label']=='stellar_zone493')
    replay = next((r for r in samples if r['zone']==493),None)
    assert replay is not None
    replay_error = max(abs(a['uncut_finite_mean']/b['uncut_finite_mean']-1) for a,b in zip(replay['interpolations'],reference['interpolations']))
    assert replay_error < 1e-14
    for p,h in inputs.items():
        assert digest(p) == h, p
    result = dict(outcome='completed_warm_envelope_comparison',accepted_for_stellar_opacity=False,
        input_sha256=inputs,composition=dict(X=x,Y3=y3,Y4=y4,Z=z,atomic_X=xa,atomic_Z=za,atomic_mass_scale=scale,
            included_atom_weights=weights,missing_mass_fraction=stars[0]['missing_metal_mass_fraction']),
        extraction=extraction,mean_extraction=mean_extraction,reader_checks=reader_checks,masked_source_states=masked,
        coverage_failures=coverage_failures,corrected_state_failures=failures,unavailable_layers=unavailable,
        reused_native_states=reused,new_native_states=len(tasks),native_corrections=list(records.values()),
        native_plasma_results=list(result_rows.values()),cases=cases,samples=samples,
        saved_layer493_replay_relative_error=replay_error,
        spectral_artifact=dict(path=str(spectrum_file),bytes=spectrum_file.stat().st_size,sha256=digest(spectrum_file)),
        limitations=['The four absent elements P, Cl, K and Ti retain mass but have zero charge and opacity in this conditional comparison.',
            'OP/MESA differences include different atomic calculations and metal distributions; they are not interpolation errors.',
            'MESA comparison maps isotope number densities and interpolates composition; it does not recover the exact 21-element TOPS mixture.',
            'Component intervals include the unspecified absorption/scattering split and signed interpolation weights only.',
            'Frequency tails, atomic data errors, bound-electron dispersion and state-interpolation errors are outside those intervals.',
            'No extrapolated mean, new opacity table, atmosphere or stellar trajectory is selected.'])
    with OUTPUT.open('x') as f:
        json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(report=str(OUTPUT),sha256=digest(OUTPUT),evaluated_layers=len(samples),
                         unavailable_layers=len(unavailable),new_states=len(tasks),reused_states=len(reused),
                         corrected_failures=len(failures),retained_spectra_bytes=spectrum_file.stat().st_size)),flush=True)


if __name__ == '__main__':
    main()
