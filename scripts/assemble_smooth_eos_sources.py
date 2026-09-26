#!/usr/bin/env python3
"""Assemble the complete helium-3 midpoint and reuse the two existing planes.

Existing potentials and sources are referenced without copying or alteration.
The new raw planes preserve every returned row and explicit absence/failure.
Per-isotherm hashes identify the actual precision variant. Numerical import
does not constitute physical acceptance or select an EOS for evolution.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from metal_eos_composition import mixture
from fetch_tops_composition import fraction_label
from eos_source_coverage import absent_source_rows,inconsistent_source_rows


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(2**20),b''):h.update(chunk)
    return h.hexdigest()


def write(path,data):
    with path.open('x') as stream:json.dump(data,stream,indent=2,allow_nan=False);stream.write('\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source_manifest',type=Path);p.add_argument('parent_manifest',type=Path)
    p.add_argument('output',type=Path);p.add_argument('--jobs',type=int,default=2)
    a=p.parse_args()
    if not 1<=a.jobs<=4:raise ValueError('bounded import requires one through four workers')
    source=json.loads(a.source_manifest.read_text());parent=json.loads(a.parent_manifest.read_text())
    spec_path=Path(source['specification'])
    if sha(spec_path)!=source['specification_sha256']:raise ValueError('collection specification changed')
    spec=json.loads(spec_path.read_text())
    if (not source['all_requested_sources_returned'] or source['failed_requests'] or
            source['unstarted_task_ids']):raise ValueError('source collection is incomplete')
    if source['source']!=spec['source']:raise ValueError('source identity differs from specification')
    for key in ['hydrogen','logT','logQ','source_coverage']:
        if spec[key]!=parent[key]:raise ValueError(f'parent and midpoint geometry differs: {key}')
    if (spec['source']['source_options']!=parent['source_options'] or
            spec['source']['source_archive_sha256']!=parent['source_archive_sha256'] or
            parent['helium3']!=[0.,.12] or spec['helium3']!=.06):
        raise ValueError('unexpected source physics or composition family')
    if source['source_nonpositive_pressure_or_zero_heat_response']:
        raise ValueError('nonpositive converged source requires separate physical review')
    xs,ts,qs=spec['hydrogen'],spec['logT'],spec['logQ'];nt,nq=len(ts),len(qs)
    requests={v['task_id']:v for v in source['requests']}
    if len(requests)!=len(source['requests']) or set(requests)!=set(range(len(xs)*nt)):
        raise ValueError('missing or duplicate source isotherms')
    a.output.mkdir();(a.output/'sources').mkdir();start=time.monotonic()
    inputs={str(p.resolve()):sha(p) for p in [a.source_manifest,a.parent_manifest,spec_path,
            Path(__file__),Path('scripts/import_freeeos_potential.py'),Path('scripts/metal_eos_composition.py'),
            Path('scripts/eos_source_coverage.py')]}
    # The old source manifest includes any tighter-quadrature isotherms used
    # for the density addition. Retain that manifest's identity separately.
    old=[];parent_root=a.parent_manifest.parent.parent
    for v in parent['planes']:
        row=dict(v)
        for key in ['potential','source']:
            path=(parent_root/v[key]).resolve()
            if sha(path)!=v[key+'_sha256']:raise ValueError('changed parent '+key)
            row[key]=os.path.relpath(path,a.output.resolve())
        row['origin']='existing_dense_family_unchanged'
        row['source_identity_manifest']=str(a.parent_manifest.resolve())
        row['source_identity_manifest_sha256']=sha(a.parent_manifest)
        old.append(row)
    write(a.output/'assembly_specification.json',dict(scope=__doc__,accepted_for_evolution=False,
          created_utc=datetime.now(timezone.utc).isoformat(),source_identity=spec['source'],
          parent_identity_manifest=str(a.parent_manifest.resolve()),inputs_sha256=inputs,
          interpolation='C2 cubic hydrogen residual; quadratic helium3 residual; exact ideal mixing'))
    def build(ix):
        x=xs[ix];expected_mixture=mixture(x,.06);rows=[];records=[];exclusions=[];flags=[]
        coverage=spec['source_coverage']
        for it,t in enumerate(ts):
            saved=requests[ix*nt+it];path=Path(saved['source_file'])
            if sha(path)!=saved['source_sha256']:raise ValueError('changed source isotherm')
            r=json.loads(gzip.decompress(path.read_bytes()))
            expected_q=[q for q in qs if q<=coverage['original_logQ_max'] or
                        (t>=coverage['minimum_added_logT'] and x<=coverage['maximum_added_hydrogen'])]
            if (r['mixture']!=expected_mixture or r['actual_probe_sha256']!=spec['source']['probe_sha256'] or
                r['source_options']!=spec['source']['source_options'] or
                r['relative_electron_integral_target']!=spec['source']['relative_electron_integral_target'] or
                r['request']['logQ']!=expected_q or r['request']['logT']!=t or
                r['request']['hydrogen']!=x or r['returncode'] or
                len(r['data'])!=len(expected_q) or
                hashlib.sha256(r['input'].encode()).hexdigest()!=r['input_sha256']):
                raise ValueError('isotherm mixture, coordinates or source identity differ')
            if any(len(v)!=22 for v in r['data']):raise ValueError('malformed source row')
            rows.extend(r['data']);rows.extend([None]*(nq-len(expected_q)))
            records.append(dict(temperature_index=it,source_file=str(path.resolve()),source_sha256=saved['source_sha256']))
            flags.extend(dict(index=it*nq+v['density_index'],info=v['info']) for v in saved['nonconvergence_flags'])
            exclusions.extend(dict(index=it*nq+v['density_index'],criterion=1e-7,reason='first_law_defect',
                              maximum_absolute_defect=v['maximum_absolute_defect']) for v in saved['first_law_failures'])
        raw=dict(version='FreeEOS 3.0.0',options=spec['source']['source_options'],**expected_mixture,
                 logT=ts,logQ=qs,data=rows,source_coverage=coverage,
                 probe_sha256=spec['source']['probe_sha256'],
                 source_archive_sha256=spec['source']['source_archive_sha256'],
                 source_identity=spec['source'],source_radiation_included=False,
                 source_consistency_exclusions=exclusions,source_isotherms=records)
        absent=absent_source_rows(raw);excluded=inconsistent_source_rows(raw)
        actual_flags=[dict(index=i,info=v[0]) for i,v in enumerate(rows) if v is not None and v[0]!=0]
        if actual_flags!=flags:raise ValueError('source failure flags do not match receipt')
        stem=f'freeeos300_gs98_x{fraction_label(x,1000)}_he3060'
        raw_path=a.output/'sources'/(stem+'_raw.json.gz');potential=a.output/(stem+'_potential.dat')
        with raw_path.open('xb') as stream:
            stream.write(gzip.compress((json.dumps(raw,separators=(',',':'),allow_nan=False)+'\n').encode(),mtime=0))
        del rows,raw
        log=a.output/(stem+'_import.log')
        with log.open('x') as stream:
            result=subprocess.run([sys.executable,'-B','scripts/import_freeeos_potential.py',str(raw_path),str(potential)],
                                  stdout=stream,stderr=subprocess.STDOUT,timeout=120)
        if result.returncode:raise RuntimeError('potential import failed; inspect '+str(log))
        return dict(hydrogen=x,helium3=.06,potential=potential.name,potential_sha256=sha(potential),
                    source=str(raw_path.relative_to(a.output)),source_sha256=sha(raw_path),
                    failed_source_states=len(flags),absent_source_states=sum(absent),
                    inconsistent_source_states=len(excluded),origin='new_midpoint_source',
                    source_identity=spec['source'],import_log=log.name,import_log_sha256=sha(log))
    new=[];failures=[]
    with (a.output/'completed_planes.jsonl').open('x') as journal,ThreadPoolExecutor(a.jobs) as pool:
        pending={pool.submit(build,i):i for i in range(len(xs))}
        for future in as_completed(pending):
            try:
                row=future.result();new.append(row);journal.write(json.dumps(row)+'\n');journal.flush()
                print(json.dumps({'completed_planes':len(new),'hydrogen':row['hydrogen']}),flush=True)
            except Exception as error:
                failures.append(dict(hydrogen=xs[pending[future]],error=repr(error)))
    if failures:
        write(a.output/'failed_assembly.json',dict(failures=failures,completed_planes=len(new),accepted_for_evolution=False))
        raise RuntimeError('assembly failed; every completed plane is retained')
    planes=sorted(old+new,key=lambda v:(v['hydrogen'],v['helium3']))
    if len(planes)!=len(xs)*3:raise ValueError('incomplete combined family')
    family=a.output/'freeeos300_gs98_z020.dat'
    family.write_text('EMBER_METAL_HELMHOLTZ 1\nhydrogen '+str(len(xs))+' '+' '.join(map(str,xs))+
        '\nhelium3 3 0.0 0.06 0.12\n'+'\n'.join(json.dumps(v['potential']) for v in planes)+'\n')
    report=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),
                outcome='assembled_validation_pending',accepted_for_evolution=False,
                family=family.name,family_sha256=sha(family),hydrogen=xs,helium3=[0.,.06,.12],
                logT=ts,logQ=qs,source_coverage=spec['source_coverage'],planes=planes,
                inputs_sha256=inputs,elapsed_seconds=time.monotonic()-start,
                retained_bytes=sum(p.stat().st_size for p in a.output.rglob('*') if p.is_file()),
                remaining_work=['Independent physical EOS/composition comparisons','Thermodynamic derivatives and mask coverage',
                                'Stellar/atmosphere queries','Runtime timing and memory','Acceptance before evolution'])
    write(a.output/'sources/freeeos300_gs98_manifest.json',report)
    print(json.dumps({k:report[k] for k in ['outcome','elapsed_seconds','retained_bytes']}),flush=True)


if __name__=='__main__':main()
