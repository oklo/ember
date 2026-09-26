#!/usr/bin/env python3
"""Test a retained-data TOPS connection through C++, without selecting it."""
import argparse
from collections import defaultdict
import gzip
import json
import math
from pathlib import Path
import subprocess

from audit_tops_electron_dispersion import KEV, KB
from fetch_tops_composition import digest
from reduce_tops_group_factors import add_inputs, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('assembly', 'work', 'report'):
        parser.add_argument('--'+name,type=Path,required=True)
    args = parser.parse_args()
    if args.work.exists() or args.report.exists():
        raise FileExistsError('preserve previous audit outputs')
    args.work.mkdir()
    assembly = json.loads(args.assembly.read_text())
    inputs = dict(assembly['input_sha256'])
    add_inputs(inputs,assembly['output_sha256'])
    probe = Path('/tmp/ember-opacity-family-probe-build-v1/opacity_family_probe')
    runtime = Path('docs/results/op_refined_compositions_v3.json')
    pinned_runtime = json.loads(runtime.read_text())
    if pinned_runtime['input_sha256'].get(str(probe)) != digest(probe):
        raise ValueError('family probe does not match retained audit')
    for p in (args.assembly,Path(__file__),probe,runtime):
        add_inputs(inputs,{str(p.resolve()):digest(p)})
    manifest = next(Path(p) for p in assembly['output_sha256'] if Path(p).name == 'tops_bridge.dat')

    def run(label, family, points):
        for p in [family]+[family.parent/line.split('"')[1] for line in family.read_text().splitlines()[1:]]:
            add_inputs(inputs,{str(p.resolve()):digest(p)})
        request = ''.join(' '.join(format(v,'.17g') for v in q)+'\n' for q in points)
        command = [str(probe),str(family)]
        reply = subprocess.run(command,input=request,text=True,capture_output=True,timeout=60)
        archive = dict(command=command,request=request,returncode=reply.returncode,stdout=reply.stdout,stderr=reply.stderr)
        (args.work/(label+'.json.gz')).write_bytes(gzip.compress(json.dumps(archive).encode(),mtime=0))
        reply.check_returncode()
        rows = [json.loads(line) for line in reply.stdout.splitlines()]
        if len(rows)!=len(points) or any(r['query']!=list(q) for r,q in zip(rows,points,strict=True)):
            raise ValueError('runtime changed queries')
        if any(r['covered'] and (not math.isfinite(r['kappa']) or r['kappa']<=0) for r in rows):
            raise ValueError('invalid supported opacity')
        return rows

    tt,rr = assembly['temperatures_keV'],assembly['densities_atomic_g_cm3']
    counts = [min(p['density_prefix_sizes'][i] for p in assembly['planes']) for i in range(len(tt))]
    safe = {t:rr[min(counts[max(0,i-2):min(len(tt),i+3)])-1] for i,t in enumerate(tt)}
    source_states = next(Path(p) for p in assembly['output_sha256'] if Path(p).name=='assembled_states.jsonl')
    nodes = [json.loads(line) for line in source_states.read_text().splitlines()]
    checked_nodes = [r for r in nodes if r['density_atomic_g_cm3']<=safe[r['temperature_keV']]]
    q = [(0,r['X'],0,r['Z'],r['temperature_keV']*KEV/KB,r['density_atomic_g_cm3']) for r in checked_nodes]
    replies = run('native_nodes',manifest,q)
    node_unsupported = sum(not r['covered'] for r in replies)
    node_max = max(abs(v['kappa']/r['kappa']-1) for r,v in zip(checked_nodes,replies) if v['covered'])
    references = []
    paths = [Path(f'docs/results/tops_cool_ratio_x{x}-z020-control_reduction_v1.json') for x in ('010','030','070')]
    paths += [Path('docs/results/tops_cool_density_pilot_x010-z020_reduction_v1.json')]
    paths += sorted(Path('docs/results').glob('tops_cool_composition_pilot_*_reduction_v2.json'))
    reference_inputs = {}
    for path in paths:
        report = json.loads(path.read_text())
        add_inputs(inputs,{str(path.resolve()):digest(path)})
        add_inputs(reference_inputs,report['input_sha256'])
        for r in report['records']:
            if tt[0]<=r['temperature_keV']<=tt[-1]:
                references.append(dict(report=str(path),X=report['X'],Z=report['Z'],**r))
    # Historical source references retain their old source hashes. Relocate only
    # byte-identical copies explicitly identified by the assembly.
    for relocation in assembly['source_relocations']:
        p=relocation['original']
        if p in reference_inputs:
            if reference_inputs[p]!=relocation['sha256']:
                raise ValueError('reference source version differs')
            del reference_inputs[p]
            reference_inputs[relocation['retained']]=relocation['sha256']
    add_inputs(inputs,reference_inputs)
    verify(inputs)
    q = [(0,r['X'],0,r['Z'],r['temperature_keV']*KEV/KB,r['density_atomic_g_cm3']) for r in references]
    values = run('independent_groups',manifest,q)
    comparisons = []
    for r,v in zip(references,values,strict=True):
        row=dict(X=r['X'],Z=r['Z'],temperature_keV=r['temperature_keV'],
            density_atomic_g_cm3=r['density_atomic_g_cm3'],source=r['source'],report=r['report'],runtime=v,
            native_uncut=r['uncut_source_rosseland'],direct_ratio=r['refractive_ratio'],
            spectral_mean=r['rosseland_atomic_cm2_g'],normalized_mean=r['native_uncut_times_ratio_atomic_cm2_g'],
            source_uncut_recovery_error=r['uncut_recovery_relative_error'])
        if v['covered']:
            row['normalized_relative_error']=v['kappa']/row['normalized_mean']-1
            row['spectral_relative_error']=v['kappa']/row['spectral_mean']-1
        comparisons.append(row)
    connection_path=Path('docs/results/opacity_connections_v1.json')
    add_inputs(inputs,{str(connection_path.resolve()):digest(connection_path)})
    connections=json.loads(connection_path.read_text())
    profile_q=[r['query'] for r in connections['rows']]
    bridge=run('profile',manifest,profile_q)
    profile=[]
    for r,v in zip(connections['rows'],bridge,strict=True):
        row=dict(zone=r['zone'],query=r['query'],bridge=v,sources=r['sources'])
        if v['covered']:
            row['relative_to_source']={name:v['kappa']/s['kappa']-1 for name,s in r['sources'].items() if s['covered']}
            means=r['mesa_original_means']
            if 'kappa_baryonic' in means[0]: row['relative_to_mesa']=v['kappa']/means[0]['kappa_baryonic']-1
        profile.append(row)
    centers=[r['query'] for r in profile if 467<=r['zone']<=483 and r['bridge']['covered']]
    # Independent composition interiors complement the actual isotope-mapped
    # stellar composition. They test derivative arithmetic, not atomic physics.
    centers += [[0,x,0,z,math.sqrt(tt[i]*tt[i+1])*KEV/KB,rho]
                for x in (.1125,.1625,.225,.55,.725) for z in (.015,.025)
                for i,rho in ((2,.03),(4,.1),(6,.3))]
    derivative_q=[]
    h=2e-5
    for q in centers:
        derivative_q.append(q)
        for axis in (4,5,1):
            for sign in (-1,1):
                altered=list(q)
                altered[axis]=q[axis]*math.exp(sign*h) if axis!=1 else q[axis]+sign*h
                derivative_q.append(altered)
    replies=run('derivatives',manifest,derivative_q)
    derivatives=[]
    for i,q in enumerate(centers):
        v=replies[7*i:7*i+7]
        if any(not r['covered'] for r in v):raise ValueError('derivative stencil unsupported')
        for j,key in enumerate(('dlnk_dlnT','dlnk_dlnrho','dlnk_dX')):
            fd=(math.log(v[2*j+2]['kappa'])-math.log(v[2*j+1]['kappa']))/(2*h)
            error=abs(fd-v[0][key])/max(1,abs(fd),abs(v[0][key]))
            derivatives.append(dict(query=q,axis=key,runtime=v[0][key],finite_difference=fd,scaled_error=error))
    guards=[[0,.173,0,.02,tt[0]*KEV/KB*.99,.1],
            [0,.173,0,.02,tt[-1]*KEV/KB*1.01,.1],
            [0,.173,0,.02,.01*KEV/KB,1e5],
            [0,.173,0,.02,.01*KEV/KB,1e-11],
            [0,.8,0,.02,.02*KEV/KB,.1]]
    guard_values=run('domain_guards',manifest,guards)
    summary=[]
    for name,predicate in [('all_retained',lambda r:True),
                          ('bridge_controls_rho_001_1_T_001_003',lambda r:.01<=r['temperature_keV']<=.03 and .01<=r['density_atomic_g_cm3']<=1)]:
        selected=[r for r in comparisons if predicate(r)]
        covered=[r for r in selected if r['runtime']['covered']]
        worst=max(covered,key=lambda r:abs(r['normalized_relative_error']))
        summary.append(dict(domain=name,queries=len(selected),covered=len(covered),
            maximum_normalized_relative_error=abs(worst['normalized_relative_error']),worst=worst,
            maximum_spectral_relative_error=max(abs(r['spectral_relative_error']) for r in covered),
            comparisons_passing_005=sum(abs(r['normalized_relative_error'])<=.005 for r in covered)))
    numerical_pass=node_unsupported==0 and node_max<1e-10 and max(r['scaled_error'] for r in derivatives)<1e-5 and not any(r['covered'] for r in guard_values)
    verify(inputs)
    result=dict(scope=__doc__,outcome='numerical_pass_physical_comparison_pending' if numerical_pass else 'numerical_failure',
        accepted_for_stellar_opacity=False,native_nodes=len(checked_nodes),native_nodes_excluded_for_stencil=len(nodes)-len(checked_nodes),
        native_nodes_unsupported=node_unsupported,maximum_node_relative_error=node_max,numerical_checks_passed=numerical_pass,
        derivative_checks=derivatives,maximum_scaled_derivative_error=max(r['scaled_error'] for r in derivatives),
        domain_guards=guard_values,independent_comparisons=comparisons,independent_summaries=summary,
        profile=profile,gap_zones_covered=[r['zone'] for r in profile if 470<=r['zone']<=480 and r['bridge']['covered']],
        limitations=['Independent references sample native temperatures, not temperature interpolation between native source rows.',
            'Reference hydrogen fractions lie on the native grey composition grid; ratio compositions and many densities are independent.',
            'The restricted control region does not establish coverage or accuracy of an evolving whole track.',
            'Source overlap and smooth blending need assessment before selection.'],
        input_sha256=inputs,output_sha256={str(p.resolve()):digest(p) for p in args.work.iterdir() if p.is_file()})
    args.report.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('outcome','native_nodes','native_nodes_unsupported','maximum_node_relative_error',
        'maximum_scaled_derivative_error','gap_zones_covered','independent_summaries')}),flush=True)
    if not numerical_pass:raise ValueError('runtime numerical check failed')


if __name__=='__main__':main()
