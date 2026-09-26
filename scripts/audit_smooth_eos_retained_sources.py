#!/usr/bin/env python3
"""Reuse the accepted base, dense and stellar FreeEOS comparison targets.

No direct source calculation is repeated. The new interpolation is compared
with every retained target using the original pressure/energy/thermal limits.
Chemical curvature is recorded where abundances are positive; comparison of
that curvature with independent source derivatives remains separate work.
"""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import numpy as np
from audit_smooth_eos_runtime import query,sha
from composition_potential_v2 import RGAS


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['probe','family','work','output']:p.add_argument(name,type=Path)
    a=p.parse_args();a.work.mkdir()
    accepted=[Path('docs/results/numerical_electron_base_acceptance_v1.json'),
              Path('docs/results/numerical_electron_hot_dense_exhaustion_acceptance_v1.json')]
    trusted={}
    for path in accepted:trusted.update(json.loads(path.read_text())['validation_reports_sha256'])
    sources=[Path('docs/results/numerical_electron_base_general_v4.json'),
             Path('docs/results/numerical_electron_hot_dense_exhaustion_general_v1.json'),
             Path('docs/results/numerical_electron_hot_dense_profile_3890gyr_v1.json')]
    inputs={str(p.resolve()):sha(p) for p in [Path(__file__),a.probe,a.family,*accepted]}
    manifest_path=a.family.parent/'sources/freeeos300_gs98_manifest.json'
    manifest=json.loads(manifest_path.read_text())
    if sha(a.family)!=manifest['family_sha256']:raise ValueError('assembled family changed')
    inputs[str(manifest_path.resolve())]=sha(manifest_path)
    for row in manifest['planes']:
        path=(a.family.parent/row['potential']).resolve()
        if sha(path)!=row['potential_sha256']:raise ValueError('potential changed')
        inputs[str(path)]=row['potential_sha256']
    records=[];points=[];indices={};groups=[]
    for path in sources:
        if sha(path)!=trusted[str(path)]:raise ValueError('retained source comparison changed')
        inputs[str(path.resolve())]=sha(path);report=json.loads(path.read_text())
        for i,r in enumerate(report['records']):
            q=tuple(r['query']);target=r.get('source',r.get('direct'))
            if q not in indices:
                indices[q]=len(points);points.append([*q,int(q[0]>0 and q[1]>0)])
            records.append(dict(source_report=str(path),source_index=i,query=list(q),
                                target=target,query_index=indices[q]))
        groups.append(dict(source_report=str(path),records=len(report['records'])))
    answers=query(a.probe,a.family,points,a.work,'retained')
    columns=dict(P=0,E=1,cv=3,cp=4,grad_ad=8)
    limits=dict(P=.001,E=.002,cv=.003,cp=.003,grad_ad=.003)
    for group in groups:
        maximum={k:0. for k in limits};failures=[]
        for r in records:
            if r['source_report']!=group['source_report']:continue
            answer=answers[r['query_index']]
            if not answer['ok']:
                failures.append(dict(query=r['query'],reason=answer['error']));continue
            values=answer['values'];error={k:values[j]/r['target'][k]-1 for k,j in columns.items()}
            r['relative_difference']=error
            for k,v in error.items():maximum[k]=max(maximum[k],abs(v))
            if any(abs(error[k])>=limits[k] for k in limits):failures.append(dict(query=r['query'],relative_difference=error))
        group.update(maximum_relative_source_differences=maximum,failed_comparisons=failures,passed=not failures)
    negative=[];minimum=None
    for q,answer in zip(points,answers):
        if not q[4] or not answer['ok']:continue
        eigen=np.linalg.eigvalsh(np.array(answer['values'][24:28]).reshape(2,2)/RGAS)
        minimum=float(eigen[0]) if minimum is None else min(minimum,float(eigen[0]))
        if eigen[0]<=0:negative.append(dict(query=q[:4],Hessian_eigenvalues_over_Rgas=eigen.tolist()))
    passed=all(g['passed'] for g in groups)
    # Recheck inputs after evaluation. Probe replies and all retained target
    # references allow a failed comparison to be diagnosed without rerunning.
    if any(sha(p)!=h for p,h in inputs.items()):raise ValueError('input changed during comparison')
    report=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),
                outcome='source_value_comparisons_passed' if passed else 'failed',
                accepted_for_evolution=False,new_FreeEOS_queries=0,retained_source_records=len(records),
                unique_runtime_queries=len(points),source_comparison_limits=limits,groups=groups,
                minimum_chemical_Hessian_eigenvalue_over_Rgas=minimum,nonpositive_chemical_curvature=negative,
                records=records,inputs_sha256=inputs,
                artifacts_sha256={str(p):sha(p) for p in a.work.iterdir() if p.is_file()},
                remaining_work=['Independent chemical-curvature checks in ionization regions',
                                'Composition boundaries and source masks','Atmosphere EOS support',
                                'Consistent diffusion forces and stellar integration'])
    with a.output.open('x') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:report[k] for k in ['outcome','retained_source_records','unique_runtime_queries',
        'minimum_chemical_Hessian_eigenvalue_over_Rgas']}))
    print(json.dumps({'groups':[{k:v for k,v in g.items() if k!='failed_comparisons'} for g in groups],
                      'nonpositive_chemical_curvature':len(negative)}))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
