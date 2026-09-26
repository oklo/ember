#!/usr/bin/env python3
"""Compare a changed native EOS against all retained chemical references.

The 102 reference gradients and Hessians were independently checked with two
source quadratures and two composition spacings. Their saved coordinates,
criteria and source data are reused without any new FreeEOS calculation.
"""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import numpy as np
from audit_smooth_eos_runtime import query,sha
from audit_diffusion_eos_two_compositions import normalized_error
from composition_potential_v2 import RGAS


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['probe','family','work','output']:p.add_argument(name,type=Path)
    p.add_argument('--cached-native',type=Path,
                   help='Replay original saved replies solely to verify this comparison helper')
    a=p.parse_args()
    reference_path=Path('docs/results/smooth_eos_chemical_domain_v1.json')
    reference=json.loads(reference_path.read_text())
    specification_path=Path('/tmp/ember-smooth-eos-chemical-domain-run-v1/calculations/specification.json')
    specification=json.loads(specification_path.read_text())
    inputs={str(path.resolve()):sha(path) for path in [Path(__file__),a.probe,a.family,
            reference_path,Path('scripts/composition_potential_v2.py'),
            Path('scripts/audit_diffusion_eos_two_compositions.py')]}
    for path,digest in reference['artifacts_sha256'].items():
        if sha(path)!=digest:raise ValueError('changed retained chemical reference: '+path)
        inputs[str(Path(path).resolve())]=digest
    if reference['source_failures'] or not all(r['source_reference_passes'] for r in reference['comparisons']):
        raise ValueError('independent source references did not pass')
    if (specification['controls']!=[{k:r[k] for k in ['material','X','Y3','hx','hy']}
                                   for r in reference['comparisons']] or
            specification['materials']!=reference['materials'] or
            specification['criteria']!=reference['criteria']):
        raise ValueError('reference coordinates or criteria changed')
    family_manifest=a.family.parent/'sources/freeeos300_gs98_manifest.json'
    manifest=json.loads(family_manifest.read_text())
    if sha(a.family)!=manifest['family_sha256']:raise ValueError('assembled family changed')
    inputs[str(family_manifest.resolve())]=sha(family_manifest)
    for row in manifest['planes']:
        path=(a.family.parent/row['potential']).resolve()
        if sha(path)!=row['potential_sha256']:raise ValueError('changed potential')
        inputs[str(path)]=row['potential_sha256']
    points=[]
    for c in specification['controls']:
        m=specification['materials'][c['material']]
        points.append([c['X'],c['Y3'],m['T'],m['rho'],1])
    a.work.mkdir()
    if a.cached_native:
        if (sha(a.probe)!=reference['inputs_sha256'][str(a.probe.resolve())] or
                sha(a.family)!=reference['inputs_sha256'][str(a.family.resolve())] or
                sha(a.cached_native)!=reference['artifacts_sha256'][str(a.cached_native)]):
            raise ValueError('cached replies belong to a different executable or family')
        answers=[json.loads(line) for line in a.cached_native.read_text().splitlines()]
    else:
        answers=query(a.probe,a.family,points,a.work,'native')
    if len(answers)!=len(points):raise ValueError('incomplete native replies')
    comparisons=[];criteria=reference['criteria']
    for c,native in zip(reference['comparisons'],answers,strict=True):
        result={k:c[k] for k in ['material','material_name','X','Y3','hx','hy']}
        result.update(source_gradient_over_Rgas=c['source_gradient_over_Rgas'],
                      source_Hessian_over_Rgas=c['source_Hessian_over_Rgas'],
                      source_reference_passes=c['source_reference_passes'])
        if not native['ok']:
            result.update(passed=False,error=native['error'])
        else:
            v=np.array(native['values']);g=np.array(c['source_gradient_over_Rgas'])*RGAS
            h=np.array(c['source_Hessian_over_Rgas'])*RGAS
            gradient=float(np.max(np.abs(v[22:24]-g))/RGAS)
            herror=normalized_error(v[24:28].reshape(2,2)-h,h)
            result.update(native_gradient_error_over_Rgas=gradient,native_Hessian_error=herror,
                          passed=bool(gradient<criteria['gradient_over_Rgas'] and
                                      herror<criteria['normalized_Hessian']))
        comparisons.append(result)
    if any(sha(path)!=digest for path,digest in inputs.items()):
        raise ValueError('input changed during comparison')
    passed=all(r['passed'] for r in comparisons)
    report=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),
                outcome='passed' if passed else 'failed',accepted_for_evolution=False,
                new_FreeEOS_queries=0,new_runtime_queries=0 if a.cached_native else len(points),
                cached_native=str(a.cached_native) if a.cached_native else None,
                criteria=criteria,comparisons=comparisons,inputs_sha256=inputs,
                artifacts_sha256={str(path):sha(path) for path in a.work.iterdir() if path.is_file()},
                remaining_work=['Thermodynamic, atmosphere and source-mask checks',
                                'Diffusion forces and conservative stellar coupling'])
    with a.output.open('x') as stream:json.dump(report,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(dict(outcome=report['outcome'],comparisons=len(comparisons),
                         passed=sum(r['passed'] for r in comparisons),new_FreeEOS_queries=0)))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
