#!/usr/bin/env python3
"""Compare variable-metal material values along a saved physical composition profile."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time
import numpy as np
from metal_eos_composition import mixture
from audit_variable_metal_eos import isotope_mixing, ARAD
from write_scientific_result import write_result


def direct(item):
    probe, q = item
    x, y, z, T, rho = q[:5]
    m = mixture(x, y, metallicity=z)
    scale = m['source_mass_scale']
    request = ' '.join(format(v, '.17g') for v in m['eps'])+'\n3 223 -2\n'+f'{math.log(rho*scale):.17g} {math.log(T):.17g}\n'
    run = subprocess.run([probe], input=request, text=True, capture_output=True, check=True, timeout=60)
    row = np.array([float(v) for v in run.stdout.split()])
    if len(row) != 22 or row[0] != 0 or not np.isfinite(row).all():
        raise ValueError(f'failed direct source: {q}')
    P, E, S = row[4], row[5]*scale, row[6]*scale-isotope_mixing(x,y,z)
    pr = ARAD*T**4/3
    pt = P+pr
    cr = P*row[7]/pt
    ct = (P*row[8]+4*pr)/pt
    delta = ct/cr
    cv = row[10]*scale/T+12*pr/(rho*T)
    cp = cv+pt/(rho*T)*ct*delta
    return dict(query=q, request=request, raw=row.tolist(),
                state=[pt,E+3*pr/rho,S+4*pr/(rho*T),cv,cp,ct,cr,delta,pt/(rho*T)*delta/cp,cr*cp/cv], phi=E/T-S)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('queries', type=Path);p.add_argument('native', type=Path)
    p.add_argument('source', type=Path);p.add_argument('output', type=Path)
    p.add_argument('--scratch', type=Path, required=True)
    p.add_argument('--jobs', type=int, default=8)
    a=p.parse_args()
    if a.output.exists() or a.scratch.exists():raise FileExistsError('preserve previous controls')
    a.scratch.mkdir(parents=True)
    sha=lambda v:hashlib.sha256(Path(v).read_bytes()).hexdigest()
    paths=[a.queries,a.native,a.source,Path(__file__),Path(__file__).with_name('metal_eos_composition.py'),Path(__file__).with_name('audit_variable_metal_eos.py')]
    hashes={str(v.resolve()):sha(v) for v in paths}
    queries=json.loads(a.queries.read_text())
    native=[json.loads(v) for v in gzip.decompress(a.native.read_bytes()).decode().splitlines()]
    if len(native)!=len(queries):raise ValueError('profile lengths disagree')
    start=time.monotonic()
    with ProcessPoolExecutor(a.jobs) as pool:
        sources=list(pool.map(direct,[(str(a.source),q) for q in queries]))
    errors=[];failures=[]
    for i,(q,n,s) in enumerate(zip(queries,native,sources)):
        if 'error' in n:failures.append(dict(zone=i,error=n['error']));continue
        values=np.array(s['state'])
        e=np.abs(np.array(n['state'])-values)/np.maximum(np.abs(values),1)
        errors.append(dict(zone=i,query=q[:5],state_errors=e.tolist(),maximum=float(max(e)),potential_error=abs(n['phi']-s['phi'])/max(abs(s['phi']),1)))
        if max(e)>1e-3:failures.append(dict(zone=i,maximum=float(max(e)),quantity=int(np.argmax(e))))
    (a.scratch/'direct_sources.json.gz').write_bytes(gzip.compress(json.dumps(sources,allow_nan=False).encode(),mtime=0))
    for v,h in hashes.items():
        if sha(v)!=h:raise ValueError(f'changed input: {v}')
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),
        outcome='passed_current_profile_source_controls' if not failures else 'failed_current_profile_source_controls',
        selected_for_evolution=False,source_queries=len(sources),elapsed_seconds=time.monotonic()-start,
        maximum_state_error=max(v['maximum'] for v in errors),
        maximum_potential_error=max(v['potential_error'] for v in errors),
        errors=errors,failures=failures,input_sha256=hashes,
        artifacts_sha256={str(v):sha(v) for v in a.scratch.iterdir()},
        limitations=['Current physical composition profile only; later structures need their own coverage checks.',
                    'State values compared with the original FreeEOS source at exact profile states; this does not independently test composition derivatives.',
                    'Fixed GS98 relative metal pattern remains an explicit material approximation.'])
    write_result(a.output,report)
    print(json.dumps({k:v for k,v in report.items() if k not in ['errors','input_sha256','artifacts_sha256']},indent=2))
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
