#!/usr/bin/env python3
"""Compare existing full EOS replies byte for byte after active-API changes.

The old probe source is unchanged. Reuse the old synthetic fixtures, complete
physical family and recorded input/output streams. No physical source is
regenerated. Replaying several physical streams in one process avoids loading
the same family repeatedly.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import subprocess


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('probe','work','report'):p.add_argument(n,type=Path)
    a=p.parse_args();a.work.mkdir();groups=[]
    base=Path('/tmp/ember-smooth-eos-four-plane-build-run-v1')
    for count in ('three','four'):
        w=base/(count+'_plane')
        groups.append((count+'_complete',w/'complete/family.dat',[w/'analytic',w/'guards']))
        groups.append((count+'_masked',w/'masked/family.dat',[w/'masks']))
        if count=='four':groups.append((count+'_helium_masked',w/'helium_masked/family.dat',[w/'helium_masks']))
    w=Path('/tmp/ember-smooth-eos-refined-physical-run-v1')
    groups.append(('physical',Path('/tmp/ember-smooth-eos-refined-family-v1/freeeos300_gs98_z020.dat'),
                   [w/'chemical/native',w/'retained/retained',w/'profile/physical']))
    comparisons=[];inputs={str(a.probe):sha(a.probe),str(Path(__file__)):sha(__file__)}
    for name,family,prefixes in groups:
        inp=a.work/(name+'.input');out=a.work/(name+'.output');err=a.work/(name+'.stderr')
        inp.write_bytes(b''.join(p.with_suffix('.input').read_bytes() for p in prefixes))
        target=b''.join(p.with_suffix('.output').read_bytes() for p in prefixes)
        for prefix in prefixes:
            for suffix in ('.input','.output'):inputs[str(prefix.with_suffix(suffix))]=sha(prefix.with_suffix(suffix))
        inputs[str(family)]=sha(family)
        for line in family.read_text().splitlines():
            if line.startswith('"'):
                table=family.parent/json.loads(line);inputs[str(table)]=sha(table)
        with inp.open('rb') as i,out.open('xb') as o,err.open('xb') as e:
            proc=subprocess.run([str(a.probe),str(family),'smooth'],stdin=i,stdout=o,stderr=e,timeout=180)
        actual=out.read_bytes();same=proc.returncode==0 and actual==target
        comparison=dict(name=name,queries=len(inp.read_bytes().splitlines()),exit_code=proc.returncode,
                        replies_identical=same,expected_sha256=hashlib.sha256(target).hexdigest(),actual_sha256=sha(out))
        if not same:
            pairs=zip(actual.splitlines(),target.splitlines())
            comparison['first_differing_replies']=[dict(row=i,actual=x.decode(),expected=y.decode()) for i,(x,y) in enumerate(pairs) if x!=y][:3]
        comparisons.append(comparison)
    ok=all(c['replies_identical'] for c in comparisons)
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),scope=__doc__,outcome='passed' if ok else 'failed',
                new_FreeEOS_queries=0,accepted_for_evolution=False,comparisons=comparisons,
                inputs_sha256=inputs,artifacts_sha256={str(q):sha(q) for q in a.work.iterdir() if q.is_file()})
    a.report.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(outcome=report['outcome'],comparisons=comparisons)))
    if not ok:raise SystemExit(1)


if __name__=='__main__':main()
