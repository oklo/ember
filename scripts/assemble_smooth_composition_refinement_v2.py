#!/usr/bin/env python3
"""Assemble the refined H/He3 family, reusing every existing composition plane.

Existing potentials and sources are referenced without copying or alteration.
The new raw planes preserve every returned row and explicit absence/failure.
Per-isotherm hashes identify the actual precision variant. Explicit validated
accuracy-correction reports may replace flagged isotherms while retaining
the original sources and their hashes. Numerical import
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
    p.add_argument('source_manifests',type=Path,nargs='+')
    p.add_argument('--parent-manifest',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--jobs',type=int,default=2)
    p.add_argument('--accuracy-corrections',type=Path,nargs='*',default=[])
    a=p.parse_args()
    if not 1<=a.jobs<=4:raise ValueError('bounded import requires one through four workers')
    collections=[json.loads(path.read_text()) for path in a.source_manifests]
    source=collections[0];parent=json.loads(a.parent_manifest.read_text())
    spec_path=Path(source['specification'])
    if sha(spec_path)!=source['specification_sha256']:raise ValueError('collection specification changed')
    spec=json.loads(spec_path.read_text())
    if (a.parent_manifest.resolve()!=Path(spec['parent_manifest']).resolve() or
            sha(a.parent_manifest)!=spec['parent_manifest_sha256']):
        raise ValueError('parent differs from the source specification')
    for path,digest in spec['inputs_sha256'].items():
        if sha(path)!=digest:raise ValueError('changed refinement input: '+path)
    for key in ['logT','logQ','source_coverage']:
        if spec[key]!=parent[key]:raise ValueError('parent and refinement geometry differs: '+key)
    if parent['helium3']!=[0.,.06,.12] or spec['target_helium3']!=[0.,.03,.06,.12]:
        raise ValueError('unexpected helium-3 family')
    xs,ys,ts,qs=spec['target_hydrogen'],spec['target_helium3'],spec['logT'],spec['logQ']
    nt,nq=len(ts),len(qs);new_planes=spec['planes']
    existing={(v['hydrogen'],v['helium3']) for v in parent['planes']}
    missing={(v['hydrogen'],v['helium3']) for v in new_planes}
    if (len(existing)!=len(parent['planes']) or len(missing)!=len(new_planes) or existing&missing or
            existing|missing!={(x,y) for x in xs for y in ys}):
        raise ValueError('existing and missing planes do not form the specified family')
    requests={}
    # Bounded allocations can finish partially. Combine their saved results
    # only when every requested isotherm is present, with no conflicting value.
    for collection in collections:
        if (collection['specification_sha256']!=source['specification_sha256'] or
                Path(collection['specification']).resolve()!=spec_path.resolve() or
                collection['source']!=spec['source']):
            raise ValueError('source identity or specification differs between allocations')
        if collection['source_nonpositive_pressure_or_zero_heat_response']:
            raise ValueError('nonpositive converged source requires separate physical review')
        seen=set()
        for record in collection['requests']:
            key=record['task_id']
            if key in seen:raise ValueError('duplicate isotherm within a collection')
            seen.add(key)
            if key in requests and record['source_sha256']!=requests[key]['source_sha256']:
                raise ValueError('conflicting source values across allocations')
            requests.setdefault(key,record)
    if set(requests)!=set(range(len(new_planes)*nt)):
        raise ValueError('source collection is incomplete; retain completed sources and finish only missing isotherms')
    corrections={}
    for report_path in a.accuracy_corrections:
        report=json.loads(report_path.read_text());identity=report['source']
        if (report['outcome']!='precision_correction_candidates_passed' or
                report['criteria']['first_law']>1e-7 or
                report['criteria']['maximum_relative_source_change']>1e-5):
            raise ValueError('accuracy correction did not pass unchanged criteria')
        for path,digest in report['inputs_sha256'].items():
            if sha(path)!=digest:raise ValueError('accuracy correction input changed')
        build=json.loads(Path(identity['build_receipt']).read_text())
        if (sha(identity['build_receipt'])!=identity['build_receipt_sha256'] or
                identity['source_options']!=spec['source']['source_options'] or
                identity['source_archive_sha256']!=spec['source']['source_archive_sha256'] or
                build['source_archive_sha256']!=identity['source_archive_sha256'] or
                build['electron_quadrature_error']!=identity['relative_electron_integral_target'] or
                build['probe_sha256']!=identity['probe_sha256'] or
                build['library_sha256']!=identity['library_sha256'] or
                identity['relative_electron_integral_target']>=spec['source']['relative_electron_integral_target']):
            raise ValueError('correction must change only source accuracy under the same physical model')
        for row in report['corrections']:
            key=row['task_id']
            if (key in corrections or key not in requests or not row['passed'] or
                    row['new_source_info_flags'] or row['new_first_law_failures'] or
                    row['original_source_sha256']!=requests[key]['source_sha256'] or
                    Path(row['original_source_file']).resolve()!=Path(requests[key]['source_file']).resolve() or
                    row['states']!=requests[key]['states']):
                raise ValueError('correction does not match original isotherm')
            corrections[key]=dict(record=row,source_identity=identity,report=str(report_path.resolve()),
                                  report_sha256=sha(report_path))
    a.output.mkdir();(a.output/'sources').mkdir();start=time.monotonic()
    inputs={str(p.resolve()):sha(p) for p in [*a.source_manifests,*a.accuracy_corrections,a.parent_manifest,spec_path,
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
        if 'import_log' in row:
            log=(parent_root/row['import_log']).resolve()
            if sha(log)!=row['import_log_sha256']:raise ValueError('changed parent import log')
            row['import_log']=os.path.relpath(log,a.output.resolve())
        # Preserve original per-plane source identity, including tighter
        # source-quadrature variants in the old dense family.
        row['reused_from_manifest']=str(a.parent_manifest.resolve())
        row['reused_from_manifest_sha256']=sha(a.parent_manifest)
        old.append(row)
    write(a.output/'assembly_specification.json',dict(scope=__doc__,accepted_for_evolution=False,
          created_utc=datetime.now(timezone.utc).isoformat(),source_identity=spec['source'],
          parent_identity_manifest=str(a.parent_manifest.resolve()),inputs_sha256=inputs,
          accuracy_corrections=corrections,
          interpolation='C2 cubic hydrogen residual; four-node cubic helium3 residual; exact ideal mixing'))
    def build(ip):
        x,y=new_planes[ip]['hydrogen'],new_planes[ip]['helium3']
        expected_mixture=mixture(x,y);rows=[];records=[];exclusions=[];flags=[];plane_corrections=[]
        coverage=spec['source_coverage']
        for it,t in enumerate(ts):
            key=ip*nt+it;saved=requests[key];actual_identity=spec['source']
            correction=corrections.get(key)
            if correction:
                candidate=correction['record'];actual_identity=correction['source_identity']
                saved={**saved,'source_file':candidate['source_file'],'source_sha256':candidate['source_sha256'],
                       'nonconvergence_flags':[],'first_law_failures':[]}
                plane_corrections.append(correction)
            path=Path(saved['source_file'])
            if sha(path)!=saved['source_sha256']:raise ValueError('changed source isotherm')
            r=json.loads(gzip.decompress(path.read_bytes()))
            expected_q=[q for q in qs if q<=coverage['original_logQ_max'] or
                        (t>=coverage['minimum_added_logT'] and x<=coverage['maximum_added_hydrogen'])]
            if (r['mixture']!=expected_mixture or r['actual_probe_sha256']!=actual_identity['probe_sha256'] or
                r['source_options']!=spec['source']['source_options'] or
                r['relative_electron_integral_target']!=actual_identity['relative_electron_integral_target'] or
                r['request']['logQ']!=expected_q or r['request']['logT']!=t or
                r['request']['hydrogen']!=x or r['request']['helium3']!=y or
                r['request']['plane_index']!=ip or r['request']['temperature_index']!=it or
                r['request']['request_index']!=ip*nt+it or
                saved['hydrogen']!=x or saved['helium3']!=y or saved['temperature_index']!=it or
                saved['plane_index']!=ip or saved['outcome']!='source_saved' or r['returncode'] or
                len(r['data'])!=len(expected_q) or
                r['input_sha256']!=saved['input_sha256'] or
                hashlib.sha256(r['input'].encode()).hexdigest()!=r['input_sha256']):
                raise ValueError('isotherm mixture, coordinates or source identity differ')
            if any(len(v)!=22 for v in r['data']):raise ValueError('malformed source row')
            rows.extend(r['data']);rows.extend([None]*(nq-len(expected_q)))
            records.append(dict(temperature_index=it,source_file=str(path.resolve()),source_sha256=saved['source_sha256'],
                                source_identity=actual_identity,accuracy_correction=correction))
            flags.extend(dict(index=it*nq+v['density_index'],info=v['info']) for v in saved['nonconvergence_flags'])
            exclusions.extend(dict(index=it*nq+v['density_index'],criterion=1e-7,reason='first_law_defect',
                              maximum_absolute_defect=v['maximum_absolute_defect']) for v in saved['first_law_failures'])
        raw=dict(version='FreeEOS 3.0.0',options=spec['source']['source_options'],**expected_mixture,
                 logT=ts,logQ=qs,data=rows,source_coverage=coverage,
                 probe_sha256=spec['source']['probe_sha256'],
                 source_archive_sha256=spec['source']['source_archive_sha256'],
                 source_identity=spec['source'],source_radiation_included=False,
                 source_consistency_exclusions=exclusions,source_isotherms=records)
        if plane_corrections:
            raw['default_probe_sha256']=raw.pop('probe_sha256')
            raw['default_source_identity']=raw.pop('source_identity')
            raw['source_identity_mode']='per_isotherm_with_explicit_accuracy_corrections'
            raw['accuracy_corrections']=plane_corrections
        absent=absent_source_rows(raw);excluded=inconsistent_source_rows(raw)
        actual_flags=[dict(index=i,info=v[0]) for i,v in enumerate(rows) if v is not None and v[0]!=0]
        if actual_flags!=flags:raise ValueError('source failure flags do not match receipt')
        stem=f'freeeos300_gs98_x{fraction_label(x,1000)}_he3{fraction_label(y,1000)}'
        raw_path=a.output/'sources'/(stem+'_raw.json.gz');potential=a.output/(stem+'_potential.dat')
        with raw_path.open('xb') as stream:
            stream.write(gzip.compress((json.dumps(raw,separators=(',',':'),allow_nan=False)+'\n').encode(),mtime=0))
        del rows,raw
        log=a.output/(stem+'_import.log')
        with log.open('x') as stream:
            result=subprocess.run([sys.executable,'-B','scripts/import_freeeos_potential.py',str(raw_path),str(potential)],
                                  stdout=stream,stderr=subprocess.STDOUT,timeout=120)
        if result.returncode:raise RuntimeError('potential import failed; inspect '+str(log))
        return dict(hydrogen=x,helium3=y,potential=potential.name,potential_sha256=sha(potential),
                    source=str(raw_path.relative_to(a.output)),source_sha256=sha(raw_path),
                    failed_source_states=len(flags),absent_source_states=sum(absent),
                    inconsistent_source_states=len(excluded),origin='new_composition_refinement_source',
                    default_source_identity=spec['source'],source_identity_mode='per_isotherm',
                    accuracy_corrections=plane_corrections,import_log=log.name,import_log_sha256=sha(log))
    new=[];failures=[]
    with (a.output/'completed_planes.jsonl').open('x') as journal,ThreadPoolExecutor(a.jobs) as pool:
        pending={pool.submit(build,i):i for i in range(len(new_planes))}
        for future in as_completed(pending):
            try:
                row=future.result();new.append(row);journal.write(json.dumps(row)+'\n');journal.flush()
                print(json.dumps({'completed_planes':len(new),'hydrogen':row['hydrogen'],'helium3':row['helium3']}),flush=True)
            except Exception as error:
                failures.append(dict(**new_planes[pending[future]],error=repr(error)))
    if failures:
        write(a.output/'failed_assembly.json',dict(failures=failures,completed_planes=len(new),accepted_for_evolution=False))
        raise RuntimeError('assembly failed; every completed plane is retained')
    planes=sorted(old+new,key=lambda v:(v['hydrogen'],v['helium3']))
    if {(v['hydrogen'],v['helium3']) for v in planes}!={(x,y) for x in xs for y in ys}:
        raise ValueError('incomplete combined family')
    family=a.output/'freeeos300_gs98_z020.dat'
    family.write_text('EMBER_METAL_HELMHOLTZ 1\nhydrogen '+str(len(xs))+' '+' '.join(map(str,xs))+
        '\nhelium3 '+str(len(ys))+' '+' '.join(map(str,ys))+'\n'+'\n'.join(json.dumps(v['potential']) for v in planes)+'\n')
    report=dict(scope=__doc__,created_utc=datetime.now(timezone.utc).isoformat(),
                outcome='assembled_validation_pending',accepted_for_evolution=False,
                family=family.name,family_sha256=sha(family),hydrogen=xs,helium3=ys,
                logT=ts,logQ=qs,source_coverage=spec['source_coverage'],planes=planes,
                inputs_sha256=inputs,accuracy_corrections=corrections,elapsed_seconds=time.monotonic()-start,
                retained_bytes=sum(p.stat().st_size for p in a.output.rglob('*') if p.is_file()),
                remaining_work=['Independent physical EOS/composition comparisons','Thermodynamic derivatives and mask coverage',
                                'Stellar/atmosphere queries','Runtime timing and memory','Acceptance before evolution'])
    write(a.output/'sources/freeeos300_gs98_manifest.json',report)
    print(json.dumps({k:report[k] for k in ['outcome','elapsed_seconds','retained_bytes']}),flush=True)


if __name__=='__main__':main()
