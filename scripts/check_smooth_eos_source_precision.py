#!/usr/bin/env python3
"""Re-evaluate only first-law-flagged isotherms at tighter source precision.

The original input and returned data are preserved. A correction is only a
candidate: all source flags, first-law identities and precision changes must
pass, and any later import must record both original and actual identities.
This does not modify the running collector or select an EOS for evolution.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import numpy as np


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('specification',type=Path);p.add_argument('work',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();spec=json.loads(a.specification.read_text());a.work.mkdir()
    source=spec['source'];receipt_path=Path(source['build_receipt'])
    receipt=json.loads(receipt_path.read_text())
    if (sha(receipt_path)!=source['build_receipt_sha256'] or
            sha(source['probe'])!=receipt['probe_sha256'] or
            sha(receipt['library'])!=receipt['library_sha256']):
        raise ValueError('tighter source build changed')
    def run(task):
        original_path=Path(task['original_source_file'])
        if sha(original_path)!=task['original_source_sha256']:raise ValueError('original source changed')
        old=json.loads(gzip.decompress(original_path.read_bytes()))
        if old['request']['request_index']!=task['task_id']:raise ValueError('source task mismatch')
        stem=a.work/f"request-{task['task_id']:05d}"
        inp=stem.with_suffix('.input');out=stem.with_suffix('.stdout');err=stem.with_suffix('.stderr')
        inp.write_text(old['input'])
        with inp.open('rb') as i,out.open('xb') as o,err.open('xb') as e:
            process=subprocess.run([source['probe']],stdin=i,stdout=o,stderr=e,timeout=120)
        result=dict(task_id=task['task_id'],original_source_file=str(original_path),
                    original_source_sha256=task['original_source_sha256'],returncode=process.returncode)
        data=[list(map(float,line.split())) for line in out.read_text().splitlines()]
        if process.returncode or len(data)!=len(old['data']) or not all(len(r)==22 and all(map(math.isfinite,r)) for r in data):
            result.update(passed=False,reason='source process or output failed');return result
        scale=old['mixture']['source_mass_scale'];flags=[];defects=[];max_defect=0.
        for index,row in enumerate(data):
            row[2]/=scale
            for col in [5,6,9,10,11,12,13,16,17]:row[col]*=scale
            if row[0]!=0:flags.append(dict(density_index=index,source_info=row[0]));continue
            if (abs(row[2]/old['data'][index][2]-1)>1e-9 or
                    abs(row[3]/old['data'][index][3]-1)>1e-10):
                raise ValueError('tighter source coordinates changed')
            rho,T,P=row[2:5];ct,er,et,sr,st=row[8:13]
            if P<=0 or et==0:raise ValueError('nonpositive source pressure or zero heat response')
            defect=max(abs(rho*er/P+ct-1),abs(T*st/et-1),abs(rho*T*sr/P+ct))
            max_defect=max(max_defect,defect)
            if defect>=spec['criteria']['first_law']:defects.append(dict(density_index=index,maximum_absolute_defect=defect))
        before=np.array(old['data']);after=np.array(data)
        columns=dict(P=4,E=5,S=6,chiRho=7,chiT=8,E_lnrho=9,E_lnT=10,S_lnrho=11,S_lnT=12)
        differences={name:float(np.max(np.abs(after[:,col]-before[:,col])/np.maximum(np.abs(before[:,col]),1.))) for name,col in columns.items()}
        record={**old,'actual_probe_sha256':receipt['probe_sha256'],
                'relative_electron_integral_target':source['relative_electron_integral_target'],
                'returncode':process.returncode,'stderr':err.read_text(),'data':data,
                'source_accuracy_correction':dict(original_source_file=str(original_path),
                    original_source_sha256=task['original_source_sha256'],source_build=source)}
        path=stem.with_suffix('.json.gz')
        with path.open('xb') as f:f.write(gzip.compress((json.dumps(record,separators=(',',':'),allow_nan=False)+'\n').encode(),mtime=0))
        passed=not flags and not defects and all(v<spec['criteria']['maximum_relative_source_change'] for v in differences.values())
        result.update(passed=passed,states=len(data),source_file=str(path),source_sha256=sha(path),
                      new_source_info_flags=flags,new_first_law_failures=defects,
                      maximum_first_law_defect=max_defect,maximum_relative_changes=differences)
        return result
    with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,spec['tasks']))
    accepted=all(v['passed'] for v in results)
    report=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),
                outcome='precision_correction_candidates_passed' if accepted else 'failed',
                accepted_for_evolution=False,specification=str(a.specification),
                specification_sha256=sha(a.specification),source=source,corrections=results,
                criteria=spec['criteria'],inputs_sha256={str(p):sha(p) for p in [Path(__file__),a.specification,receipt_path]},
                artifacts_sha256={str(p):sha(p) for p in a.work.iterdir() if p.is_file()},
                remaining_work=['Explicit source-identity-preserving import of any accepted corrections',
                                'Physical and thermodynamic comparisons of the complete EOS family'])
    with a.output.open('x') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(outcome=report['outcome'],corrections=results)))
    if not accepted:raise SystemExit(1)


if __name__=='__main__':main()
