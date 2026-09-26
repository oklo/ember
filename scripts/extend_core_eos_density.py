#!/usr/bin/env python3
"""Append actual hot, dense FreeEOS states while preserving all old source rows.

This prepares a candidate family. Independent source, interpolation and finite
stellar checks are required before using it for evolution.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import shlex
import subprocess
import sys

from eos_source_coverage import absent_source_rows, inconsistent_source_rows


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['manifest', 'probe', 'work', 'report']:
        parser.add_argument(name, type=Path)
    parser.add_argument('--maximum-log-q', type=float, default=4.3)
    parser.add_argument('--jobs', type=int, default=4)
    args = parser.parse_args()
    if args.work.exists() or not 1 <= args.jobs <= 8:
        raise ValueError('Use a fresh directory and one to eight workers')
    args.work.mkdir()
    lines = args.manifest.read_text().splitlines()
    if lines[0] != 'EMBER_VARIABLE_METAL_HELMHOLTZ 1':
        raise ValueError('Expected variable-metal source family')
    potentials = [(args.manifest.parent/shlex.split(line)[0]).resolve() for line in lines[4:]]
    unique = list(dict.fromkeys(potentials))
    source_records = []
    probe_sha = digest(args.probe)
    for potential in unique:
        name = potential.name
        if name.startswith('potential-'):
            source = potential.parent/('plane-'+name[10:13])/'source.json.gz'
        elif name == 'potential.dat':
            source = potential.parent/'plane-000/source.json.gz'
        else:
            source = potential.with_name(potential.stem+'-explicit-masks.json.gz')
        with potential.open() as stream:
            next(stream)
            header = next(stream)
        expected = re.search(r'direct-source SHA256 ([a-f0-9]{64})', header)
        if expected is None or digest(source) != expected.group(1):
            raise ValueError(f'Raw source identity differs: {potential}')
        source_records.append(dict(potential=str(potential), source=str(source),
                                   potential_sha256=digest(potential), source_sha256=digest(source)))
    plan = dict(created_utc=datetime.now(timezone.utc).isoformat(), jobs=args.jobs,
                maximum_log_q=args.maximum_log_q, unique_planes=len(unique),
                selected_for_evolution=False, sources=source_records,
                input_sha256={str(p):digest(p) for p in [args.manifest,args.probe,Path(__file__),
                    Path(__file__).with_name('import_freeeos_potential.py'),
                    Path(__file__).with_name('eos_source_coverage.py')]})
    save(args.work/'plan.json', plan)

    def extend(item):
        index, record = item
        work = args.work/f'plane-{index:03d}'
        work.mkdir()
        old = json.loads(gzip.decompress(Path(record['source']).read_bytes()))
        if old['probe_sha256'] != probe_sha or old['options'] != [3,223,-2]:
            raise ValueError('Different source executable or physical options')
        absent_source_rows(old)
        inconsistent_source_rows(old)
        coverage = old['source_coverage']
        ts, qs = old['logT'], old['logQ']
        step = (qs[-1]-qs[0])/(len(qs)-1)
        count = round((args.maximum_log_q-qs[-1])/step)
        if count < 4 or abs(qs[-1]+count*step-args.maximum_log_q) > 1e-10:
            raise ValueError('New density range must align with the old grid')
        added_q = [qs[-1]+i*step for i in range(count+1)]
        new_q = qs+added_q[1:]
        scale = old['source_mass_scale']
        rows, overlaps, requests = [], [], []
        for it, t in enumerate(ts):
            retained = old['data'][it*len(qs):(it+1)*len(qs)]
            if (t < coverage['minimum_added_logT'] or
                    old['hydrogen'] > coverage.get('maximum_added_hydrogen',1.)):
                rows.extend(retained+[None]*count)
                continue
            request = (' '.join(map(str,old['eps']))+'\n3 223 -2\n'+''.join(
                f'{math.log(scale)+math.log(10)*(q+1.5*(t-6)):.17g} {math.log(10)*t:.17g}\n'
                for q in added_q))
            result = subprocess.run([str(args.probe)],input=request,text=True,
                                    capture_output=True,timeout=120)
            if result.returncode:
                (work/f'temperature-{it}.failure.log').write_text(result.stdout+'\n'+result.stderr)
                raise ValueError(f'FreeEOS process failed at logT={t}')
            added = [list(map(float,line.split())) for line in result.stdout.splitlines()]
            if len(added) != len(added_q) or any(len(r)!=22 or not all(map(math.isfinite,r)) for r in added):
                raise ValueError('Incomplete or nonfinite new source response')
            for r in added:
                r[2] /= scale
                for k in [5,6,9,10,11,12,13,16,17]:
                    r[k] *= scale
            left,right = retained[-1],added[0]
            difference = None
            if left[0] == right[0] == 0:
                difference = max(abs(a-b)/max(1.,abs(a),abs(b)) for a,b in zip(left[2:],right[2:]))
                if difference > 1e-8:
                    raise ValueError(f'Independent source overlap differs at logT={t}: {difference}')
            overlaps.append(dict(logT=t,flags=[left[0],right[0]],maximum_scaled_difference=difference))
            requests.append(dict(logT=t,input_sha256=hashlib.sha256((probe_sha+request).encode()).hexdigest(),
                                 stderr=result.stderr))
            rows.extend(retained+added[1:])
        assert len(rows) == len(ts)*len(new_q)
        for it in range(len(ts)):
            assert rows[it*len(new_q):it*len(new_q)+len(qs)] == old['data'][it*len(qs):(it+1)*len(qs)]
        raw = dict(old, logQ=new_q, data=rows,
                   core_density_extension=dict(retained_source_sha256=record['source_sha256'],
                       source_probe_sha256=probe_sha,new_logQ=added_q[1:],independent_overlap=overlaps,
                       requests=requests))
        exclusions=[]
        for exclusion in old.get('source_consistency_exclusions',[]):
            it,iq=divmod(exclusion['index'],len(qs))
            exclusions.append(dict(exclusion,index=it*len(new_q)+iq))
        # Preserve source flags and values. Only measured failures of the
        # unchanged first-law criterion may mask an entire derivative stencil.
        for it in range(len(ts)):
            for iq in range(len(qs),len(new_q)):
                row=rows[it*len(new_q)+iq]
                if row is None or row[0] != 0:
                    continue
                rho,T,P=row[2:5];chit,Er,Et,Sr,St=row[8:13]
                defect=max(abs(rho*Er/P+chit-1),abs(T*St/Et-1),abs(rho*T*Sr/P+chit))
                if defect>1e-7:
                    exclusions.append(dict(index=it*len(new_q)+iq,criterion=1e-7,
                        reason='first_law_defect',maximum_absolute_defect=defect))
        raw['source_consistency_exclusions']=exclusions
        absent_source_rows(raw)
        inconsistent_source_rows(raw)
        source=work/'source.json.gz'
        source.write_bytes(gzip.compress(json.dumps(raw,separators=(',',':'),allow_nan=False).encode(),mtime=0))
        potential=work/'potential.dat'
        proc=subprocess.run([sys.executable,str(Path(__file__).with_name('import_freeeos_potential.py')),
                             str(source),str(potential)],capture_output=True,text=True)
        (work/'import.log').write_text(proc.stdout+proc.stderr)
        if proc.returncode:
            raise ValueError('Potential import failed; see retained log')
        checks=[r['maximum_scaled_difference'] for r in overlaps if r['maximum_scaled_difference'] is not None]
        if not checks:
            raise ValueError('No valid independent overlap source')
        return dict(index=index,previous=record,source=str(source),source_sha256=digest(source),
                    potential=str(potential),potential_sha256=digest(potential),
                    old_rows_retained_exactly=True,maximum_overlap_difference=max(checks),
                    excluded_source_rows=len(exclusions),
                    flagged_source_rows=sum(r is not None and r[0]!=0 for r in rows))

    results,failures=[],[]
    with ThreadPoolExecutor(args.jobs) as pool:
        jobs={pool.submit(extend,item):item[0] for item in enumerate(source_records)}
        for future in as_completed(jobs):
            index=jobs[future]
            try:
                results.append(future.result())
                print('completed',len(results),'/',len(unique),'plane',index,flush=True)
            except Exception as error:
                failures.append(dict(index=index,error=repr(error)))
                print('failed',index,repr(error),flush=True)
    manifest=None
    if not failures:
        paths={r['previous']['potential']:r['potential'] for r in results}
        manifest=args.work/'variable_metal.dat'
        manifest.write_text('\n'.join(lines[:4]+[json.dumps(paths[str(p)]) for p in potentials])+'\n')
    for p,h in plan['input_sha256'].items():
        assert digest(p)==h,p
    save(args.report,dict(created_utc=datetime.now(timezone.utc).isoformat(),
         outcome='completed' if not failures else 'completed_with_failures',selected_for_evolution=False,
         planes=results,failures=failures,manifest=None if manifest is None else str(manifest),
         manifest_sha256=None if manifest is None else digest(manifest),input_sha256=plan['input_sha256'],
         plan_sha256=digest(args.work/'plan.json'),
         limitations=['Only actual hot dense source states appended; old source values unchanged.',
                     'Independent runtime, derivative and finite stellar controls required before selection.',
                     'Cold dense states remain outside the declared source coverage.']))


if __name__ == '__main__':
    main()
