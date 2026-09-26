#!/usr/bin/env python3
"""Compare independently mixed OP compositions with a generated C++ family.

Mixtures used as family nodes are never counted as independent comparisons.
Retain each additional calculation separately so a paused audit can resume.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

from generate_op_warm_family import plane_task, atomic_json

ROOT=Path(__file__).resolve().parents[1]
PROBE=Path('/tmp/ember-opacity-family-probe-build-v1/opacity_family_probe')
SNAPSHOT=Path('/tmp/ember-op-family-source-v1/manifest.json')
ENVELOPE=ROOT/'docs/results/op_warm_envelope_v1.json'
OLD=ROOT/'docs/results/op_warm_runtime_v1.json'


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def key(row):return row['temperature_index'],row['X'],row['Z']


def plane_path(work,task):
    return work/'planes'/('mixture_'+hashlib.sha256(json.dumps(task,separators=(',',':')).encode()).hexdigest()+'.json')


def query(work,label,manifest,queries,inputs):
    source=''.join(' '.join(str(v) for v in q)+'\n' for q in queries)
    (work/(label+'_queries.txt')).write_text(source)
    reply=subprocess.run([str(PROBE),str(manifest)],input=source,text=True,capture_output=True,timeout=45)
    (work/(label+'_stdout.jsonl')).write_text(reply.stdout)
    (work/(label+'_stderr.txt')).write_text(reply.stderr)
    reply.check_returncode()
    rows=[json.loads(line) for line in reply.stdout.splitlines()]
    assert len(rows)==len(queries) and all(r['query']==list(q) for r,q in zip(rows,queries))
    inputs[str(manifest)]=digest(manifest)
    for line in manifest.read_text().splitlines()[1:]:
        p=manifest.parent/line.split('"')[1];inputs[str(p)]=digest(p)
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--family',type=Path,required=True)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();work=args.work
    if args.output.exists():raise FileExistsError(args.output)
    work.mkdir(exist_ok=True);(work/'planes').mkdir(exist_ok=True)
    inputs={str(p):digest(p) for p in [Path(__file__),args.family,PROBE,SNAPSHOT,ENVELOPE,OLD]}
    snap=json.loads(SNAPSHOT.read_text());assert digest(PROBE)==snap['probe_sha256']
    for row in snap['sources']:
        assert digest(row['source'])==row['sha256'];inputs[row['source']]=row['sha256']
    family=json.loads(args.family.read_text())
    assert family['outcome']=='assembled_refined_conditional_family'
    for p,h in family['input_sha256'].items():assert digest(p)==h;inputs[p]=h
    for f in family['table_files']:assert digest(f['path'])==f['sha256']
    config_path=Path(family['configuration'])
    assert digest(config_path)==family['configuration_sha256']
    inputs[str(config_path)]=digest(config_path)
    config=json.loads(config_path.read_text());xs,zs=config['X'],config['Z']
    tables=config_path.parent/'density0025'
    env=json.loads(ENVELOPE.read_text());c=env['composition']
    actual_x=c['atomic_X'];actual_z=c['atomic_Z']
    midpoints=[.5*(a+b) for a,b in zip(xs,xs[1:])]
    its=[176,184,192,200]
    tasks=[(it,x,z) for it in its for z in zs+[actual_z] for x in midpoints]
    tasks += [(it,actual_x,actual_z) for it in its]
    assert len(set(tasks))==len(tasks)
    training={tuple(r['task']) for r in family['source_records']}
    assert all(t not in training for t in tasks),'training point included in independent comparisons'
    configuration=dict(family=str(args.family),family_sha256=digest(args.family),temperatures=its,
        tasks=tasks,target_relative_error=.005,input_sha256=inputs.copy())
    configuration_path=work/'configuration.json'
    # JSON normalization preserves exact float coordinates while making
    # resumption comparisons independent of tuple/list serialization.
    configuration=json.loads(json.dumps(configuration))
    if configuration_path.exists():assert json.loads(configuration_path.read_text())==configuration
    else:atomic_json(configuration_path,configuration)
    prior=json.loads(OLD.read_text())
    direct_file=Path('/tmp/ember-op-warm-runtime-v1/withheld_composition_means.json')
    assert digest(direct_file)==prior['input_sha256'][str(direct_file)]
    inputs[str(direct_file)]=digest(direct_file)
    retained={key(r):r for r in json.loads(direct_file.read_text())}
    direct=[];pending=[];records=[]
    for task in tasks:
        path=plane_path(work,task)
        if task in retained:
            direct.append(retained[task]);records.append(dict(task=task,kind='retained_independent_mixture',path=str(direct_file)))
        elif path.exists():
            row=json.loads(path.read_text());assert key(row)==task;direct.append(row)
            records.append(dict(task=task,kind='resumed_independent_mixture',path=str(path),sha256=digest(path)))
        else:pending.append(task)
    reused=len(direct);start=time.monotonic()
    print(json.dumps(dict(total=len(tasks),reused=reused,new=len(pending))),flush=True)
    with ProcessPoolExecutor(max_workers=8) as pool:
        futures={pool.submit(plane_task,t):t for t in pending}
        for future in as_completed(futures):
            task=futures[future];row=future.result();assert key(row)==task
            path=plane_path(work,task);atomic_json(path,row);direct.append(row)
            records.append(dict(task=task,kind='new_independent_mixture',path=str(path),sha256=digest(path)))
            if (len(direct)-reused)%36==0 or len(direct)==len(tasks):
                progress=dict(completed=len(direct),total=len(tasks),elapsed_seconds=time.monotonic()-start)
                atomic_json(work/'progress.json',progress);print(json.dumps(progress),flush=True)
    direct.sort(key=lambda r:(r['temperature_index'],r['Z'],r['X']))
    atomic_json(work/'withheld_composition_means.json',direct)
    inputs[str(work/'withheld_composition_means.json')]=digest(work/'withheld_composition_means.json')
    queries=[];source_values=[];invalid=[]
    for row in direct:
        if row['failures'] or row['status']!='evaluated':invalid.append(dict(task=key(row),status=row['status'],failures=row['failures']))
        for s in row['native_states']:
            if config['logR_min']<=s['logR']<=config['logR_max']:
                queries.append((0,row['X'],0.,row['Z'],10**row['logT'],s['rho']))
                source_values.append(dict(X=row['X'],Z=row['Z'],logT=row['logT'],logR=s['logR'],source=s))
    response=query(work,'withheld',tables/'op_gs98_ordinary.dat',queries,inputs)
    held=[];unsupported=[]
    for a,b in zip(source_values,response):
        if not b['covered']:
            unsupported.append(dict(**{k:a[k] for k in ['X','Z','logT','logR']},runtime=b));continue
        held.append(dict(**a,runtime=b,relative_opacity_difference=b['kappa']/a['source']['ordinary']-1))
    intervals=[]
    for a,b in zip(xs,xs[1:]):
        for z in zs+[actual_z]:
            rows=[r for r in held if r['X']==.5*(a+b) and r['Z']==z]
            if not rows:raise ValueError('no supported comparisons in interval')
            worst=max(rows,key=lambda r:abs(r['relative_opacity_difference']))
            intervals.append(dict(X_lower=a,X_upper=b,Z=z,samples=len(rows),
                maximum_relative_difference=abs(worst['relative_opacity_difference']),
                worst={k:worst[k] for k in ['X','Z','logT','logR','relative_opacity_difference']}))
    profile=[(1,c['X'],c['Y3'],c['Z'],r['temperature_K'],r['density_g_cm3']) for r in env['samples']]
    h=1e-5;hx=1e-6;shifted=[]
    for q in profile:
        shifted.append(q)
        for index in [4,5,1]:
            for sign in [-1,1]:
                trial=list(q)
                if index==1:trial[index]+=sign*hx
                else:trial[index]*=math.exp(sign*h)
                shifted.append(tuple(trial))
    current=query(work,'profile',tables/'op_gs98_ordinary.dat',shifted,inputs)
    absorption=query(work,'absorption',tables/'op_gs98_absorption.dat',profile,inputs)
    scattering=query(work,'scattering',tables/'op_gs98_scattering.dat',profile,inputs)
    original=query(work,'original',Path('/tmp/ember-op-warm-family-v1/density0025/op_gs98_ordinary.dat'),profile,inputs)
    samples=[];derivatives=[];uncovered=[]
    for i,q in enumerate(profile):
        group=current[7*i:7*i+7];zone=env['samples'][i]['zone'];a=group[0]
        if not all(r['covered'] for r in group):uncovered.append(zone);continue
        for index,name,delta in [(1,'dlnk_dlnT',h),(3,'dlnk_dlnrho',h),(5,'dlnk_dX',hx)]:
            fd=(math.log(group[index+1]['kappa'])-math.log(group[index]['kappa']))/(2*delta)
            derivatives.append(dict(zone=zone,derivative=name,analytic=a[name],numerical=fd,
                normalized_difference=abs(fd-a[name])/max(1.,abs(fd),abs(a[name]))))
        samples.append(dict(zone=zone,query=q,runtime=a,original=original[i],absorption=absorption[i],scattering=scattering[i],
            relative_to_original=a['kappa']/original[i]['kappa']-1,
            endpoint_order_preserved=all(r['covered'] for r in [absorption[i],scattering[i]]) and
                a['kappa']<=scattering[i].get('kappa',-1)<=absorption[i].get('kappa',-1)))
    max_h=max(abs(r['relative_opacity_difference']) for r in held if r['Z'] in zs)
    max_actual_z=max(abs(r['relative_opacity_difference']) for r in held if r['Z']==actual_z)
    max_derivative=max(r['normalized_difference'] for r in derivatives)
    result=dict(outcome='completed_independent_composition_comparison',accepted_for_stellar_opacity=False,
        input_sha256=inputs,configuration=str(configuration_path),configuration_sha256=digest(configuration_path),
        independent_plane_records=records,independent_planes=len(direct),reused_planes=reused,new_planes=len(pending),
        native_queries=len(queries),heldout_supported=held,heldout_unsupported=unsupported,invalid_planes=invalid,
        interval_errors=intervals,maximum_hydrogen_relative_difference=max_h,
        maximum_actual_metallicity_relative_difference=max_actual_z,hydrogen_interpolation_target=.005,
        hydrogen_interpolation_passed=max_h<.005 and not invalid,
        actual_metallicity_interpolation_passed=max_actual_z<.005 and not invalid,
        profile_uncovered=uncovered,profile_samples=samples,derivative_checks=derivatives,
        maximum_derivative_normalized_difference=max_derivative,derivative_check_passed=max_derivative<1e-5 and not uncovered,
        limitations=['Independent midpoints at four native temperatures test composition interpolation; they do not bound all thermodynamic locations or atomic-data accuracy.',
            'Queries outside the complete temperature/density/composition stencil remain excluded rather than extrapolated.',
            'Density discretization, source transitions, missing elements and partially ionized transport retain their separate requirements.',
            'No selected stellar input, atmosphere or evolutionary state is changed.'])
    for p,hsh in inputs.items():assert digest(p)==hsh,p
    atomic_json(args.output,result)
    print(json.dumps(dict(report=str(args.output),sha256=digest(args.output),supported=len(held),unsupported=len(unsupported),
        maximum_hydrogen_error=max_h,maximum_actual_metallicity_error=max_actual_z,maximum_derivative_error=max_derivative)),flush=True)


if __name__=='__main__':main()
