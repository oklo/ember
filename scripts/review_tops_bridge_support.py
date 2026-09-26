#!/usr/bin/env python3
"""Resolve the bridge audit's outer-temperature stencil accounting from saved replies.

The first audit counted three rows at the high-temperature endpoint, while
the C++ interpolation uses four. No query is rerun and no opacity is changed.
"""
import argparse
import gzip
import json
from pathlib import Path

from audit_tops_electron_dispersion import KEV, KB
from fetch_tops_composition import digest
from reduce_tops_group_factors import verify


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    assembly_path=Path('docs/results/tops_bridge_v2.json')
    previous_path=Path('docs/results/tops_bridge_runtime_v1.json')
    assembly=json.loads(assembly_path.read_text());previous=json.loads(previous_path.read_text())
    verify(previous['input_sha256']);verify(previous['output_sha256'])
    raw_path=Path('/tmp/ember-tops-bridge-audit-v1/native_nodes.json.gz')
    raw=json.loads(gzip.decompress(raw_path.read_bytes()))
    replies=[json.loads(line) for line in raw['stdout'].splitlines()]
    temps=assembly['temperatures_keV'];densities=assembly['densities_atomic_g_cm3']
    counts=[min(r['density_prefix_sizes'][j] for r in assembly['planes']) for j in range(len(temps))]
    limits={}
    for i,t in enumerate(temps):
        # At a knot either adjacent interval can be chosen after roundoff.
        rows=set()
        for interval in (i-1,i):
            start=max(0,min(interval-1,len(temps)-4))
            rows.update(range(start,start+4))
        limits[t*KEV/KB]=densities[min(counts[j] for j in rows)-1]
    classifications=[]
    for row in replies:
        q=row['query'];expected=q[5]<=limits[q[4]]
        classifications.append((expected,row['covered']))
    disagreements=sum(expected!=covered for expected,covered in classifications)
    supported=sum(expected for expected,covered in classifications)
    proper_rejections=sum(not expected and not covered for expected,covered in classifications)
    numerical=(disagreements==0 and previous['maximum_node_relative_error']<1e-10
        and previous['maximum_scaled_derivative_error']<1e-5
        and not any(r['covered'] for r in previous['domain_guards']))
    result=dict(scope=__doc__,outcome='supported_runtime_checks_pass' if numerical else 'unresolved_runtime_failure',
        accepted_for_stellar_opacity=False,numerical_checks_passed=numerical,reused_native_queries=len(replies),
        supported_native_queries=supported,properly_rejected_native_queries=proper_rejections,
        support_classification_disagreements=disagreements,new_runtime_queries=0,
        maximum_node_relative_error=previous['maximum_node_relative_error'],
        derivative_comparisons=len(previous['derivative_checks']),
        maximum_scaled_derivative_error=previous['maximum_scaled_derivative_error'],
        independent_summaries=previous['independent_summaries'],gap_zones_covered=previous['gap_zones_covered'],
        resolution='The earlier failure was the audit support prediction at the hottest endpoint; the saved C++ rejections agree with the full four-temperature stencil.',
        limitations=['The original audit and its failed receipt remain preserved.',
            'Density interpolation errors remain above the 0.5% target; this does not select a stellar opacity.',
            'For future native-node audits use the union of both clamped four-temperature stencils, including at both endpoints.'],
        input_sha256={str(p.resolve()):digest(p) for p in (Path(__file__),assembly_path,previous_path,raw_path)})
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('scope','independent_summaries','input_sha256','limitations')}),flush=True)
    if not numerical:raise ValueError('support classification remains unresolved')


if __name__=='__main__':main()
