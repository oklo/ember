#!/usr/bin/env python3
"""Native EOS and face controls for radiation carried by species redistribution."""
import argparse
import bisect
import csv
import gzip
import json
import math
import os
import selectors
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from audit_full_star_diffusion import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--scratch', type=Path, required=True)
    parser.add_argument('--reuse', action='append', type=Path, default=[],
                        help='Reuse verified native replies from a previous scratch directory.')
    args = parser.parse_args()
    assert not args.output.exists() and not args.scratch.exists()
    args.scratch.mkdir()
    root = Path(__file__).resolve().parents[1]
    eos = Path('/tmp/ember-smooth-eos-refined-family-v1/freeeos300_gs98_z020.dat')
    collision = Path('/tmp/ember-collision-transport-table-v1.dat')
    probe = Path('/tmp/ember-parcel-heat-probe-v1')
    conduction = root / 'data/conduction/condtab21wd_metals.dat'
    profile = root / 'docs/reports/2026-09-11/evolution_latest_profile.csv'
    source_path = root / 'docs/results/native_composition_heat_v2.json'
    source = json.loads(source_path.read_text())
    for name, expected in source['artifacts_sha256'].items(): assert digest(name) == expected, name
    previous = Path('/tmp/ember-native-heat-audit-v2')
    old_queries = gzip.decompress((previous / 'queries.txt.gz').read_bytes()).decode().splitlines()
    old_replies = [json.loads(s) for s in gzip.decompress((previous / 'responses.jsonl.gz').read_bytes()).decode().splitlines()]
    identities = {str(p): digest(p) for p in [Path(__file__).resolve(), probe, eos, collision, conduction,
                    profile, source_path, root/'scripts/parcel_heat_probe.cpp', root/'scripts/conditional_envelope_heat.hpp',
                    root/'src/eos_smooth_mixture.cpp',root/'include/ember/eos_smooth_mixture.hpp',
                    root/'src/screened_microscopic_transport.cpp',root/'include/ember/screened_microscopic_transport.hpp',
                    Path('/tmp/ember-parcel-radiation-build-v1/src/libember.a')]}
    identities.update(source['artifacts_sha256'])
    stamps = {p:(Path(p).stat().st_size,Path(p).stat().st_mtime_ns) for p in identities}
    with profile.open() as stream: rows = [{k:float(v) for k,v in r.items()} for r in csv.DictReader(stream)]
    command = [str(probe),str(eos),str(collision),str(conduction)]
    (args.scratch/'inputs.json').write_text(json.dumps(dict(command=command,input_sha256=identities),indent=2)+'\n')
    qfile=(args.scratch/'queries.txt').open('x');raw=(args.scratch/'responses.jsonl').open('x');err=(args.scratch/'stderr.txt').open('x')
    checks=[];failures=[];centers=[];query_count=0;cache={};reused_count=0
    for directory in args.reuse:
        if directory.name=='ember-parcel-face-width-v1':
            reuse_report=root/'docs/results/parcel_heat_width_diagnostic_v1.json'
        else:
            assert directory.name=='ember-parcel-heat-audit-v1'
            reuse_report=root/'docs/results/native_parcel_heat_v1.json'
        retained=json.loads(reuse_report.read_text())
        for p,h in retained['artifacts_sha256'].items(): assert digest(p)==h,p
        if (directory/'inputs.json').exists():
            old_inputs=json.loads((directory/'inputs.json').read_text())['input_sha256']
        else:old_inputs=retained['input_sha256']
        for p,h in old_inputs.items():
            actual=Path('/tmp/ember-parcel-heat-audit-source-v1.py') if Path(p)==Path(__file__).resolve() else Path(p)
            assert digest(actual)==h,p
        for line,reply in zip((directory/'queries.txt').read_text().splitlines(),(directory/'responses.jsonl').read_text().splitlines(),strict=True):
            cache[line+'\n']=json.loads(reply)
        identities[str(reuse_report)]=digest(reuse_report)
        identities.update(retained['artifacts_sha256'])
    # A derivative of the interpolant is local to its EOS cell. The quintic
    # potential has continuous second derivatives, but not necessarily third
    # derivatives used by these heat responses. Keep the FD stencil in-cell.
    plane=root/'data/eos/numerical_electron_hot_dense_exhaustion_v1/freeeos300_gs98_x000_he3000_potential.dat'
    axes={}
    with plane.open() as stream:
        for line in stream:
            if line.startswith('data'):break
            if line.startswith(('log_t ','log_q ')):
                key,_,*values=line.split();axes[key]=[float(x)*math.log(10) for x in values]
    def eos_cell(q):
        lt=.5*(q[7]+q[13]);lr=math.log(.5*(math.exp(q[6])+math.exp(q[12])))
        return (bisect.bisect_right(axes['log_t'],lt),bisect.bisect_right(axes['log_q'],lr-1.5*(lt-6*math.log(10))))
    (args.scratch/'inputs.json').write_text(json.dumps(dict(command=command,input_sha256=identities),indent=2)+'\n')
    start=time.monotonic();proc=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=err,text=True,bufsize=1)
    (args.scratch/'native_process.json').write_text(json.dumps(dict(pid=proc.pid,controller_pid=os.getpid(),command=command,start_utc=datetime.now(timezone.utc).isoformat(),maximum_elapsed_seconds=300),indent=2)+'\n')
    selector=selectors.DefaultSelector();selector.register(proc.stdout,selectors.EVENT_READ)
    def ask(mode,values):
        nonlocal query_count,reused_count
        line=mode+' '+' '.join(format(float(v),'.17g') for v in values)+'\n'
        if line in cache:
            reused_count+=1
            return cache[line]
        qfile.write(line);qfile.flush();proc.stdin.write(line);proc.stdin.flush()
        if not selector.select(120):raise TimeoutError('No native heat response within 120 seconds.')
        reply=proc.stdout.readline()
        if not reply:raise RuntimeError('Native heat process ended early.')
        raw.write(reply);raw.flush();record=json.loads(reply);query_count+=1
        if 'error' in record:raise RuntimeError(record['error']+'; '+line)
        cache[line]=record;return record
    def check(kind,error,tolerance,context):
        value=float(error);item=dict(kind=kind,error=value,tolerance=tolerance,context=context,
                                     passed=math.isfinite(value) and value<=tolerance)
        checks.append(item)
        if not item['passed']:failures.append(item)
    def relative(a,b):
        a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float);valid=np.isfinite(b)
        assert np.array_equal(np.isfinite(a),valid)
        return float(np.max(abs(a[valid]-b[valid])/np.maximum(abs(b[valid]),1))) if valid.any() else 0.
    try:
        # Preserve material values while extending the response with radiation.
        for i,(q,old) in enumerate(zip(old_queries,old_replies)):
            new=ask('eos',list(map(float,q.split())))
            check('material_value_preserved',relative(new['material'],old['enthalpy']),2e-13,i)
            check('material_partials_preserved',relative(new['d_material'],old['enthalpy_partials']),2e-12,i)
        print(json.dumps(dict(retained_EOS_replies=len(old_replies),queries=query_count,failed_checks=len(failures))),flush=True)
        # Resolve composition derivatives without trace-species subtraction noise.
        for node in (0,160,300,397,440,480,511):
            row=rows[node];q=[row['X'],max(row['Y3'],.001),row['temperature_K'],row['density_g_cm3'],1,1]
            base=ask('eos',q);analytic=np.asarray(base['d_radiation']);hr=np.asarray(base['radiation'])
            centers.append(dict(kind='eos',node=node,query=q,response=base))
            for coordinate in range(4):
                for width in ((1e-4,5e-5) if coordinate<2 else (1e-5,5e-6)):
                    values=[]
                    for delta in (-2,-1,1,2):
                        trial=q.copy()
                        if coordinate<2:trial[coordinate+2]*=math.exp(delta*width)
                        else:trial[coordinate-2]+=delta*width
                        values.append(np.asarray(ask('eos',trial)['radiation']))
                    numerical=(values[0]-8*values[1]+8*values[2]-values[3])/(12*width)
                    scale=np.maximum(np.max(abs(analytic),axis=1),abs(hr))
                    check('EOS_radiation_derivative',np.max(abs(numerical-analytic[:,coordinate])/scale),3e-6,[node,coordinate,width])
            # Independent isobaric density solves differentiate total and radiation enthalpy.
            for species in (0,1):
                values=[];width=1e-4
                for delta in (-2,-1,1,2):
                    trial=q[:2];trial[species]+=delta*width
                    reply=ask('bulk',[*trial,q[2],base['P']])
                    check('isobaric_pressure',abs(reply['P']/base['P']-1),1e-10,[node,species,delta])
                    values.append(np.array([reply['h_rad'],reply['h_bulk']]))
                numerical=(values[0]-8*values[1]+8*values[2]-values[3])/(12*width)
                expected=np.array([hr[species],hr[species]+base['material'][species]])
                check('isobaric_bulk_enthalpy',np.max(abs(numerical-expected)/np.maximum(abs(expected),1)),3e-6,[node,species])
        def point(row):return [math.log(row['radius_cm']),math.log(row['density_g_cm3']),math.log(row['temperature_K']),row['luminosity_erg_s'],row['X'],row['Y3']]
        face_cases=[]
        for combined,faces in [(0,[0,300,396,410,420]),(1,[396,410,420,430,470,510])]:
            for face in faces:
                q=[combined,1,3,rows[face]['mass_g'],rows[face+1]['mass_g'],*point(rows[face]),*point(rows[face+1]),5e12,-1e11]
                face_cases.append((face,q))
        for combined,face in [(0,300),(1,470)]:
            for mask in (0,1,2):
                q=[combined,1,mask,rows[face]['mass_g'],rows[face+1]['mass_g'],*point(rows[face]),*point(rows[face+1]),5e12,-1e11]
                for k in (0,1):
                    if not mask&(1<<k):q[9+k]=q[15+k]=q[17+k]=0
                face_cases.append((face,q))
        def gate(lt):
            z=(lt-math.log(2e6))/math.log(1.5)
            return 0 if z<=0 else 1 if z>=1 else z**3*(10+z*(-15+6*z))
        for face,q in face_cases:
            base=ask('face',q);total=np.asarray(q[17:19]);weight=gate(q[7])*gate(q[13]) if q[0] else 1.
            if q[0] and weight>0:
                hot=q.copy();hot[0]=0;micro=np.asarray(ask('face',hot)['micro'])
            else:micro=np.asarray(base['micro'])
            mid=[.5*(q[9]+q[15]),.5*(q[10]+q[16]),math.exp(.5*(q[7]+q[13])),.5*(math.exp(q[6])+math.exp(q[12])),bool(q[2]&1),bool(q[2]&2)]
            h=ask('eos',mid);hr=np.nan_to_num(np.asarray(h['radiation'],dtype=float),nan=0.)
            expected=float(hr@(total-weight*micro));delta=base['corrected']['Q']-base['old']['Q']
            scale=max(abs(expected),sum(abs(hr*(total-weight*micro))),1.)
            check('parcel_heat_split',abs(delta-expected)/scale,3e-9,[q[0],face,q[2]])
            check('total_rate_coefficient',relative(np.asarray(base['corrected']['h_total'])-base['old']['h_total'],hr),3e-9,[q[0],face,q[2]])
            check('conductivity_unchanged',relative([base['corrected']['K'],*base['corrected']['dKlo'],*base['corrected']['dKhi']],
                                                       [base['old']['K'],*base['old']['dKlo'],*base['old']['dKhi']]),0,[q[0],face,q[2]])
            centers.append(dict(kind='face',face=face,query=q,extra_heat=delta,expected_extra_heat=expected,weight=weight,response=base))
            matrix=np.array([base['corrected']['dQlo']+base['corrected']['dQhi'],
                             base['old']['dQlo']+base['old']['dQhi']])
            added=matrix[0]-matrix[1]
            for end in (0,1):
                for coordinate in range(3):
                    width=1e-4;cell=eos_cell(q)
                    for refinement in range(20):
                        trials=[]
                        for step in (-2,-1,1,2):
                            trial=q.copy();trial[1]=0;trial[5+6*end+coordinate]+=step*width;trials.append(trial)
                        if all(eos_cell(t)==cell for t in trials):break
                        width*=.5
                    else:raise ValueError('No finite in-cell derivative stencil.')
                    widths=[width,width/2] if width<1e-4 else [width]
                    for width in widths:
                        values=[]
                        for step in (-2,-1,1,2):
                            trial=q.copy();trial[1]=0;trial[5+6*end+coordinate]+=step*width
                            value=ask('face',trial);values.append(np.array([value['corrected']['Q'],value['corrected']['Q']-value['old']['Q']]))
                        numerical=(values[0]-8*values[1]+8*values[2]-values[3])/(12*width)
                        target=np.array([matrix[0,4*end+coordinate],added[4*end+coordinate]])
                        scale=np.array([max(abs(matrix[0]).max(),1.),max(abs(added).max(),1.)])
                        check('face_heat_derivative',np.max(abs(numerical-target)/scale),3e-5,[q[0],face,q[2],end,coordinate,width])
        proc.stdin.close();proc.wait(timeout=10);assert proc.returncode==0
    except Exception as exc:failures.append(dict(exception=type(exc).__name__,message=str(exc)))
    finally:
        if proc.poll() is None:
            proc.terminate()
            try:proc.wait(timeout=3)
            except subprocess.TimeoutExpired:proc.kill();proc.wait(timeout=3)
        selector.close();qfile.close();raw.close();err.close()
    for p,stamp in stamps.items():
        if (Path(p).stat().st_size,Path(p).stat().st_mtime_ns)!=stamp and digest(p)!=identities[p]:failures.append(dict(changed_input=p))
    maxima={}
    for c in checks:maxima[c['kind']]=max(maxima.get(c['kind'],0),c['error'])
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed_native_parcel_heat_controls' if not failures else 'failed_native_parcel_heat_controls',
                accepted_for_stellar_evolution=False,queries=query_count,reused_queries=reused_count,checks=len(checks),maxima=maxima,failures=failures,centers=centers,
                refined_face_derivatives=[c for c in checks if c['kind']=='face_heat_derivative' and c['context'][-1]<1e-4],
                elapsed_seconds=time.monotonic()-start,native_exit=proc.returncode,
                artifacts_sha256={str(p):digest(p) for p in args.scratch.iterdir()},
                limitations=['The extra term assumes LTE, isobaric parcels and linear enthalpy contrasts.',
                             'Cool microscopic drift remains omitted in the conditional wrapper.',
                             'This audit does not validate moving convective boundaries or extend the selected trajectory.'])
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('outcome','queries','checks','maxima','failures','elapsed_seconds')},indent=2),flush=True)
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
