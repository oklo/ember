#!/usr/bin/env python3
"""Check the C++ OP family and its source/composition interpolation.

Use independent composition mixtures at native OP temperatures and densities,
then check runtime derivatives and neighboring opacity sources on the star.
No evolving stellar model or selected atmosphere is changed.
"""
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import subprocess

from generate_op_warm_family import plane_task, XGRID

ROOT=Path(__file__).resolve().parents[1]
WORK=Path('/tmp/ember-op-warm-runtime-v1')
FAMILY=ROOT/'docs/results/op_warm_family_v1.json'
ENVELOPE=ROOT/'docs/results/op_warm_envelope_v1.json'
OUTPUT=ROOT/'docs/results/op_warm_runtime_v1.json'
PROBE=Path('/tmp/ember-opacity-family-probe-build-v1/opacity_family_probe')
TABLES=Path('/tmp/ember-op-warm-family-v1')
SOURCE=Path('/tmp/ember-refractive-hot-family-refined-v1')


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def query(label,manifest,queries,inputs):
    text=''.join(' '.join(str(v) for v in q)+'\n' for q in queries)
    qfile=WORK/(label+'_queries.txt');qfile.write_text(text)
    reply=subprocess.run([str(PROBE),str(manifest)],input=text,text=True,capture_output=True,timeout=30)
    (WORK/(label+'_stdout.jsonl')).write_text(reply.stdout)
    (WORK/(label+'_stderr.txt')).write_text(reply.stderr)
    reply.check_returncode()
    rows=[json.loads(line) for line in reply.stdout.splitlines()]
    assert len(rows)==len(queries) and all(r['query']==list(q) for r,q in zip(rows,queries))
    inputs[str(manifest)]=digest(manifest)
    for line in manifest.read_text().splitlines()[1:]:
        p=manifest.parent/line.split('"')[1];inputs[str(p)]=digest(p)
    return rows


def main():
    if OUTPUT.exists() or WORK.exists():raise FileExistsError('runtime comparison exists')
    WORK.mkdir()
    inputs={str(p):digest(p) for p in [Path(__file__),PROBE,FAMILY,ENVELOPE,
        ROOT/'scripts/generate_op_warm_family.py',ROOT/'scripts/opacity_family_probe.cpp',
        ROOT/'src/opacity_mixture.cpp',ROOT/'src/opacity_table.cpp',ROOT/'src/interp.cpp',ROOT/'src/composition.cpp']}
    family=json.loads(FAMILY.read_text());assert family['outcome']=='assembled_conditional_warm_family'
    for p,h in family['input_sha256'].items():assert digest(p)==h
    for f in family['table_files']:assert digest(f['path'])==f['sha256']
    old=json.loads(ENVELOPE.read_text());c=old['composition']
    # Every H interval is tested, with extra metallicity comparisons at two
    # H fractions. These mixtures were not used to construct the tables.
    tasks=[(it,.5*(a+b),.02) for it in [176,184,192,200] for a,b in zip(XGRID,XGRID[1:])]
    tasks += [(it,x,z) for it in [176,192] for x in [.175,.5] for z in [.015,.025]]
    with ProcessPoolExecutor(max_workers=8) as pool:direct=list(pool.map(plane_task,tasks))
    direct_file=WORK/'withheld_composition_means.json';direct_file.write_text(json.dumps(direct,indent=2,allow_nan=False)+'\n')
    inputs[str(direct_file)]=digest(direct_file)
    held_queries=[];held_values=[]
    for row in direct:
        for s in row['native_states']:
            if -3<=s['logR']<=2.25:
                held_queries.append((0,row['X'],0.,row['Z'],10**row['logT'],s['rho']))
                held_values.append(dict(X=row['X'],Z=row['Z'],logT=row['logT'],logR=s['logR'],source=s))
    response=query('withheld',TABLES/'density0025/op_gs98_ordinary.dat',held_queries,inputs)
    held=[];unsupported=[]
    for a,b in zip(held_values,response):
        if not b['covered']:unsupported.append({k:a[k] for k in ['X','Z','logT','logR']});continue
        held.append({**a,'runtime':b,'relative_opacity_difference':b['kappa']/a['source']['ordinary']-1})
    profile_queries=[(1,c['X'],c['Y3'],c['Z'],r['temperature_K'],r['density_g_cm3']) for r in old['samples']]
    h=1e-5;queries=[]
    for q in profile_queries:
        queries.append(q)
        for index in [4,5]:
            for sign in [-1,1]:
                shifted=list(q);shifted[index]*=math.exp(sign*h);queries.append(tuple(shifted))
    fine=query('fine',TABLES/'density0025/op_gs98_ordinary.dat',queries,inputs)
    coarse=query('coarse',TABLES/'density005/op_gs98_ordinary.dat',profile_queries,inputs)
    absorption=query('absorption',TABLES/'density0025/op_gs98_absorption.dat',profile_queries,inputs)
    scattering=query('scattering',TABLES/'density0025/op_gs98_scattering.dat',profile_queries,inputs)
    molecular=query('molecular',SOURCE/'aesopus21_gs98_mixture.dat',profile_queries,inputs)
    tops_low=query('tops_low',SOURCE/'tops_gs98_mixture_low.dat',profile_queries,inputs)
    sampled=[];derivatives=[];uncovered=[]
    for i,q in enumerate(profile_queries):
        group=fine[5*i:5*i+5];a=group[0]
        zone=old['samples'][i]['zone']
        if not all(s['covered'] for s in group):
            uncovered.append(zone);continue
        for index,key in [(1,'dlnk_dlnT'),(3,'dlnk_dlnrho')]:
            fd=(math.log(group[index+1]['kappa'])-math.log(group[index]['kappa']))/(2*h)
            derivatives.append(dict(zone=zone,derivative=key,analytic=a[key],numerical=fd,
                                    normalized_difference=abs(fd-a[key])/max(1.,abs(fd),abs(a[key]))))
        out=dict(zone=zone,T=q[4],rho=q[5],ordinary=a['kappa'],dlnk_dlnT=a['dlnk_dlnT'],dlnk_dlnrho=a['dlnk_dlnrho'],
                 coarse=coarse[i],absorption=absorption[i],scattering=scattering[i],molecular=molecular[i],tops_low=tops_low[i])
        if coarse[i]['covered']:out['coarse_relative_to_fine']=coarse[i]['kappa']/a['kappa']-1
        if molecular[i]['covered']:out['op_relative_to_molecular']=a['kappa']/molecular[i]['kappa']-1
        if tops_low[i]['covered']:out['op_relative_to_tops_low']=a['kappa']/tops_low[i]['kappa']-1
        if absorption[i]['covered'] and scattering[i]['covered']:
            out['endpoint_order_preserved']=a['kappa']<=scattering[i]['kappa']<=absorption[i]['kappa']
        sampled.append(out)
    max_derivative=max(r['normalized_difference'] for r in derivatives)
    result=dict(outcome='completed_runtime_comparison',accepted_for_stellar_opacity=False,input_sha256=inputs,
        heldout_composition_planes=len(direct),heldout_native_queries=len(held_queries),heldout_supported=held,
        heldout_unsupported=unsupported,profile_uncovered=uncovered,profile_samples=sampled,derivative_checks=derivatives,
        maximum_derivative_normalized_difference=max_derivative,
        maximum_heldout_relative_difference=max(abs(r['relative_opacity_difference']) for r in held),
        maximum_density_regrid_difference=max(abs(r['coarse_relative_to_fine']) for r in sampled if 'coarse_relative_to_fine' in r),
        derivative_check_passed=max_derivative<1e-5,
        limitations=['Withheld comparisons measure interpolation within the adopted OP physics; they do not measure absolute atomic-data accuracy.',
            'H fractions only cover the generated 0–0.75 range; later hydrogen-rich diffusion atmospheres require additional composition coverage.',
            'Molecular and TOPS connections are comparisons at saved states; neither smooth blended transport nor complete stellar evolution is tested here.',
            'Source masks and rejected queries remain explicit. No selected stellar opacity is changed.'])
    for p,hsh in inputs.items():assert digest(p)==hsh,p
    with OUTPUT.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(report=str(OUTPUT),sha256=digest(OUTPUT),heldout=len(held),profile_covered=len(sampled),
        derivative=max_derivative,density_regrid=result['maximum_density_regrid_difference'],
        composition_error=result['maximum_heldout_relative_difference'])),flush=True)


if __name__=='__main__':main()
