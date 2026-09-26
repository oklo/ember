#!/usr/bin/env python3
"""Review retained Claude/Fable response artifacts without rerunning simulations."""
import argparse
import hashlib
import json
from pathlib import Path


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    root=Path('docs/research/fable/results')
    names=['fable-thermal-phase-correction-v1','fable-envelope-transport-v1',
           'fable-envelope-transport-v2','fable-mode-damping-v1','fable-encounter-rate-v1']
    batches=[];inputs={str(Path(__file__).resolve()):digest(__file__)}
    for name in names:
        directory=root/name;ready=directory/'READY.json';data=json.loads(ready.read_text())
        inputs[str(ready.resolve())]=digest(ready)
        checks=[]
        for relative,expected in data['files_sha256'].items():
            p=directory/relative
            actual=digest(p)
            checks.append(dict(path=str(p),expected=expected,actual=actual,matched=actual==expected))
            inputs[str(p.resolve())]=actual
        batches.append(dict(batch=name,files=len(checks),hash_checks=checks,all_hashes_match=all(r['matched'] for r in checks)))
    thermal=json.loads((root/names[0]/'READY.json').read_text())
    rates=json.loads((root/names[-1]/'encounter_outcomes_v1.json').read_text())
    partitions=[]
    for r in rates['scenarios']:
        total=r['rate_total_closed_per_yr'];parts=sum(r['rate_'+name+'_per_yr'] for name in ('contact','capture','survivor'))
        residual=parts/total-1
        partitions.append(dict(scenario=r['scenario'],unequal=r['unequal'],total=total,partition_sum=parts,
            relative_residual=residual,requested_relative_tolerance=1e-6,passes_requested_tolerance=abs(residual)<=1e-6))
    report=dict(scope=__doc__,accepted_new_stellar_physics=False,batches=batches,
        all_file_hashes_match=all(r['all_hashes_match'] for r in batches),
        files_checked=sum(r['files'] for r in batches),
        thermal_inverse=dict(source_review='Explicit liquid fraction spans the latent-energy interval and reconstructs energy at fixed melting temperature.',
            retained_worker_controls=len(thermal['tests']),retained_controls_all_pass=thermal['all_tests_passed'],
            primary_reran_controls=False,physical_liquid_heat_capacity_accepted=False),
        encounter_partitions=partitions,maximum_absolute_partition_residual=max(abs(r['relative_residual']) for r in partitions),
        envelope_review=dict(latest_reviewed_batch='fable-envelope-transport-v2',accepted=False,
            outstanding=['Use retained faces397–420 rather than asserting their kinetic data are absent.',
                'Correct rho*D*w scaling when w is already a velocity.',
                'Neutral cross section, drift scale and transported enthalpy need sensitivity tests before being called bounds.',
                'Compare conductive prescriptions and conserve both sides of every internal face.']),
        next_directives='docs/research/fable/PRIMARY_NOTES.md P088–P090',
        limitations=['File identity and scalar outcome closure do not independently validate SPH trajectories.',
            'No new particle simulation or source microphysics calculation was run.'],input_sha256=inputs)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(files_checked=report['files_checked'],all_hashes_match=report['all_file_hashes_match'],
        maximum_partition_residual=report['maximum_absolute_partition_residual'],
        partition_tolerance_failures=sum(not r['passes_requested_tolerance'] for r in partitions))),flush=True)
    if not report['all_file_hashes_match']:raise ValueError('sealed artifact changed')


if __name__=='__main__':main()
