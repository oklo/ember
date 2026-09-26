#!/usr/bin/env python3
"""Test warmer OP spectra at the current envelope composition and four corners.

Reuse the existing physical spectral evaluator without changing its source or
its configuration files. New source members come from the retained OP archive.
This pilot does not establish a complete composition family or select opacity.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.interpolate import PchipInterpolator
import generate_op_warm_family as physical
from audit_op_warm_envelope import extract_requested
from reduce_tops_group_factors import verify


def initialize(work):
    physical.WORK=Path(work)
    physical.configuration.cache_clear()
    physical.native_elements.cache_clear()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--work',type=Path,required=True)
    p.add_argument('--report',type=Path,required=True)
    a=p.parse_args()
    if a.work.exists() or a.report.exists():raise FileExistsError('preserve pilot outputs')
    a.work.mkdir();(a.work/'planes').mkdir()
    old=Path('/tmp/ember-op-warm-family-v1/configuration.json')
    config=json.loads(old.read_text())
    warm=Path('docs/results/op_warm_envelope_v1.json')
    composition=json.loads(warm.read_text())['composition']
    archive=Path('/tmp/ember-op-monochromatic-fetch-v2/OP4STARS_1.3.tar.xz')
    inputs={str(q.resolve()):physical.digest(q) for q in (Path(__file__),old,warm,archive,
        Path('scripts/audit_op_warm_envelope.py'),Path('scripts/generate_op_warm_family.py'),
        Path('scripts/op_native_spectra_v2.py'),Path('scripts/op_mixture_corrections.py'),
        Path('scripts/op_plasma_transport.py'),Path('scripts/audit_tops_electron_dispersion.py'),physical.MESH)}
    verify(config['input_sha256'])
    elements=sorted(next(iter(config['source_paths'].values())),key=int)
    temps=list(range(208,232,2))
    requests={f'OP4STARS_1.3/mono/m{int(e):02d}.{it}':a.work/'spectra'/f'm{int(e):02d}.{it}' for it in temps for e in elements}
    extracted=extract_requested(archive,'aeae2b31e62c7cebc100be2813e9b976de0681a31faa5fa2716c766cc1a6e809',requests)
    if sum(r['bytes'] for r in extracted)>250000000:raise ValueError('source storage cap exceeded')
    for r in extracted:inputs[r['path']]=r['sha256']
    config['source_paths']={str(it):{e:str(requests[f'OP4STARS_1.3/mono/m{int(e):02d}.{it}']) for e in elements} for it in temps}
    config['temperatures']=temps
    config['input_sha256']=inputs.copy()
    physical.atomic_json(a.work/'configuration.json',config)
    mixtures=[(composition['atomic_X'],composition['atomic_Z']), (0.,.01),(0.,.03),(.75,.01),(.75,.03)]
    tasks=[(it,x,z) for it in temps for x,z in mixtures]
    rows=[];start=time.monotonic()
    with ProcessPoolExecutor(max_workers=4,initializer=initialize,initargs=(str(a.work),)) as pool:
        futures={pool.submit(physical.plane_task,task):task for task in tasks}
        for future in as_completed(futures):
            task=futures[future];row=future.result();rows.append(row)
            name=physical.hashlib.sha256(json.dumps(task).encode()).hexdigest()+'.json'
            physical.atomic_json(a.work/'planes'/name,row)
            if len(rows)%10==0:print(json.dumps(dict(completed=len(rows),total=len(tasks),elapsed_seconds=time.monotonic()-start)),flush=True)
    rows.sort(key=lambda r:(r['temperature_index'],r['X'],r['Z']))
    unsupported=[];summaries=[]
    for r in rows:
        try:
            states=physical.connected_states(r)
            summaries.append(dict(temperature_index=r['temperature_index'],X=r['X'],Z=r['Z'],
                positive_native_states=len(states),logR_min=states[0]['logR'],logR_max=states[-1]['logR'],
                source_masks=r['source_masks'],failures=r['failures']))
        except ValueError as error:unsupported.append(dict(task=[r['temperature_index'],r['X'],r['Z']],reason=str(error)))
    connection_path=Path('docs/results/opacity_connections_v1.json')
    inputs[str(connection_path.resolve())]=physical.digest(connection_path)
    connection=json.loads(connection_path.read_text())
    profile=sorted(connection['rows'],key=lambda r:r['query'][4])
    logt=np.log10([r['query'][4] for r in profile]);logrho=np.log10([r['query'][5] for r in profile])
    locus=PchipInterpolator(logt,logrho,extrapolate=False)
    comparisons=[];queries=[]
    for row in rows:
        if (row['X'],row['Z'])!=mixtures[0]:continue
        t=10**(.025*row['temperature_index']);rho=10**float(locus(math.log10(t)))
        logr=math.log10(rho*composition['atomic_mass_scale'])-3*math.log10(t)+18
        try:states=physical.connected_states(row)
        except ValueError:continue
        values={field:float(10**PchipInterpolator([r['logR'] for r in states],np.log10([r[field] for r in states]),extrapolate=False)(logr))*composition['atomic_mass_scale']
                for field in ('ordinary','absorption','scattering')}
        if not all(math.isfinite(v) and v>0 for v in values.values()):continue
        comparisons.append(dict(T=t,rho=rho,logR=logr,op=values))
        queries.append([1,composition['X'],composition['Y3'],composition['Z'],t,rho])
    probe=Path('/tmp/ember-opacity-family-probe-build-v1/opacity_family_probe')
    known=json.loads(Path('docs/results/op_refined_compositions_v3.json').read_text())
    if known['input_sha256'][str(probe)]!=physical.digest(probe):raise ValueError('changed source probe')
    inputs[str(probe)]=physical.digest(probe)
    for label,manifest in [('tops_bridge',Path('/tmp/ember-tops-bridge-v2/tables/tops_bridge.dat')),
                           ('tops_hot',Path('/tmp/ember-refractive-hot-family-refined-v1/tops_gs98_mixture_high.dat'))]:
        for q in [manifest]+[manifest.parent/line.split('"')[1] for line in manifest.read_text().splitlines()[1:]]:
            inputs[str(q.resolve())]=physical.digest(q)
        text=''.join(' '.join(format(v,'.17g') for v in q)+'\n' for q in queries)
        reply=subprocess.run([str(probe),str(manifest)],input=text,text=True,capture_output=True,timeout=30)
        physical.atomic_json(a.work/(label+'_runtime.json'),dict(request=text,stdout=reply.stdout,stderr=reply.stderr,returncode=reply.returncode))
        reply.check_returncode();answers=[json.loads(line) for line in reply.stdout.splitlines()]
        if len(answers)!=len(queries):raise ValueError('missing runtime replies')
        for row,answer in zip(comparisons,answers):
            row[label]=answer
            if answer['covered']:row['op_ordinary_relative_to_'+label]=row['op']['ordinary']/answer['kappa']-1
    verify(inputs)
    result=dict(scope=__doc__,outcome='completed_conditional_pilot',accepted_for_stellar_opacity=False,
        source_members_extracted=len(extracted),source_bytes=sum(r['bytes'] for r in extracted),
        new_atomic_calculations=0,new_mixture_isotherms=len(rows),native_states=sum(len(r['native_states']) for r in rows),
        failed_native_states=sum(len(r['failures']) for r in rows),unsupported_planes=unsupported,
        planes=summaries,source_records=[dict(path=str(q),sha256=physical.digest(q)) for q in sorted((a.work/'planes').glob('*.json'))],
        current_composition_comparisons=comparisons,elapsed_mixture_seconds=time.monotonic()-start,
        limitations=['Five compositions are a pilot, not the full hydrogen/metallicity family.',
            'The density comparison follows an interpolated saved profile; no stellar evolution was performed.',
            'OP missing trace elements and plasma-transport alternatives retain their prior conditional status.',
            'No independent temperature/composition interpolation acceptance follows from native isotherms.'],
        input_sha256=inputs,output_sha256={str(q.resolve()):physical.digest(q) for q in a.work.rglob('*') if q.is_file()})
    physical.atomic_json(a.report,result)
    print(json.dumps({k:result[k] for k in ('outcome','source_members_extracted','source_bytes','new_mixture_isotherms','native_states','failed_native_states','unsupported_planes')}),flush=True)


if __name__=='__main__':main()
