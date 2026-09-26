#!/usr/bin/env python3
"""Check native CN quadrature against independent Python integration and measure its stellar source."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('output', type=Path)
    ap.add_argument('--scratch', type=Path, required=True)
    a = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    a.scratch.mkdir()
    assert not a.output.exists()
    native = Path('/tmp/ember-cn-native-build-v1/cn_probe')
    worker = root / 'docs/research/fable/results/fable-cno-rate-v4/cno_rate_v2.py'
    checkpoint = Path('/tmp/ember-adaptive-diffusion-v7/checkpoint.json')
    identities = {}

    def pin(p):
        p = Path(p)
        identities[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()

    for p in (native, worker, checkpoint, __file__, root/'src/nuclear_pp.cpp', root/'include/ember/nuclear.hpp'):
        pin(p)
    spec = importlib.util.spec_from_file_location('independent_cn_reference', worker)
    ref = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ref)
    model = np.asarray(json.loads(checkpoint.read_text())['model_record']['model'])
    queries = []
    # Atomic-basis checks match the worker's documented abundance convention.
    for T, rho, X, Y3 in ((5e6,300.,.5,.003), (8e6,1000.,.1,.01), (1.32e7,1e4,.0008,1e-9), (1.9e7,300.,.01,.002)):
        queries.append([0,1,T,rho,X,Y3,1-X-Y3-.02])
    for carbon in (0.,1.):
        queries.extend([1,carbon,T,rho,X,Y3,1-X-Y3-.02] for m,r,rho,T,L,X,Y3 in model)
    raw = ''.join(' '.join(format(v,'.17g') for v in row)+'\n' for row in queries)
    (a.scratch/'queries.txt').write_text(raw)
    result = subprocess.run([str(native)], input=raw, text=True, capture_output=True, timeout=60, check=True)
    (a.scratch/'responses.jsonl').write_text(result.stdout)
    (a.scratch/'stderr.txt').write_text(result.stderr)
    records = [json.loads(line) for line in result.stdout.splitlines()]
    assert len(records) == len(queries) and not any('error' in r for r in records)
    comparisons = []
    for q, r in zip(queries[:4], records[:4]):
        basis, carbon, T, rho, X, Y3, Y4 = q
        bare, slope = ref.bare_rate(T,ref.SF3['N14(p,g)'],npts=24001)
        screen = ref.screening(T,rho,X,Y3,ref.SF3['N14(p,g)'])
        catalyst = .02*(.171836/12+.050335/14)/ref.GS98_SCALE
        molar_cycles = rho*bare*math.exp(screen['logf'])*X/ref.A_H*catalyst
        q_molar = (4*ref.A_H-ref.A_HE4)*ref.C_L**2
        expected = dict(bare_rate=bare,bare_slope=slope,screen=screen['logf'],zeta=screen['zeta'],
                        eps_cn=molar_cycles*(q_molar-ref.Q_NU*ref.MEV*ref.NA),
                        neutrino_cn=molar_cycles*ref.Q_NU*ref.MEV*ref.NA)
        differences = {k:abs(r[k]/v-1) for k,v in expected.items()}
        comparisons.append(dict(query=q,expected=expected,measured=r,relative_differences=differences))
    maximum = max(v for row in comparisons for v in row['relative_differences'].values())
    assert maximum < 2e-8
    m=model[:,0];weights=np.empty(len(m));weights[0]=.5*(m[0]+m[1]);weights[-1]=.5*(m[-1]-m[-2]);weights[1:-1]=.5*(m[2:]-m[:-2])
    stars = {}
    for j,name in enumerate(('nitrogen_only','CN')):
        rs=records[4+j*len(m):4+(j+1)*len(m)]
        sums={k:float(np.dot(weights,[r[k] for r in rs])) for k in ('eps_cn','eps_pp','neutrino_cn','dXH_cn','dXHe4_cn')}
        stars[name]=dict(integrated_sources=sums,ratio_CN_to_actual_pp=sums['eps_cn']/sums['eps_pp'],
                         max_zeta=max(r['zeta'] for r in rs),central_sources=rs[0])
    for p in a.scratch.iterdir():
        if p.is_file(): pin(p)
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed',native_queries=len(queries),
                maximum_independent_relative_difference=maximum,independent_comparisons=comparisons,
                stellar_sources=stars,input_sha256=identities,
                limitations=['Fixed GS98 catalyst approximation; prior carbon-conversion fuel is not restored.',
                             'Source luminosity on an unchanged structure; a coupled stellar response is required.',
                             'No oxygen conversion, catalyst diffusion or CN leakage is included.'])
    a.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('outcome','native_queries','maximum_independent_relative_difference','stellar_sources')},indent=2))


if __name__ == '__main__':
    main()
