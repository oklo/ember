#!/usr/bin/env python3
"""Calculate reusable gas-opacity planes and atmosphere columns at specified Z.

The first batch supplies low-metal source comparisons at two hydrogen
abundances. It does not by itself define an interpolated atmosphere family
or select a composition approximation for stellar evolution.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import generate_nongrey_grid as gen
from assemble_nongrey_grid import load_continuation
from nongrey_opacity import validate_table, merge_isotherms
from prepare_nongrey_sources import digest, data_digest
from write_scientific_result import write_result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('plan',type=Path);p.add_argument('work',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--jobs',type=int,default=8)
    a=p.parse_args()
    if a.output.exists() or a.work.exists():raise FileExistsError('use a new batch directory')
    if not 1<=a.jobs<=8:raise ValueError('one to eight workers required')
    a.work.mkdir(parents=True);root=Path(__file__).resolve().parents[1]
    plan=json.loads(a.plan.read_text());prepared_path=Path(plan['prepared']);prepared=json.loads(prepared_path.read_text())
    base=json.loads(Path(plan['reference_provenance']).read_text())['specification']
    for k in ['synspec','tlusty']:
        if digest(prepared[k])!=prepared['executables'][k]:raise ValueError('changed source executable')
    for path,h in zip(prepared['line_lists'],prepared['line_list_sha256'],strict=True):
        if digest(path)!=h:raise ValueError('changed line list')
    if data_digest(Path(prepared['synple'])/'data')!=prepared['data_sha256']:raise ValueError('changed source data')
    sources=[Path(__file__),a.plan,prepared_path,Path(plan['reference_provenance'])]+[
        root/'scripts'/v for v in ['generate_nongrey_grid.py','run_nongrey_continuation.py','assemble_nongrey_grid.py','import_nongrey_grid.py']]
    hashes={str(v.resolve()):digest(v) for v in sources}
    (a.work/'scripts').mkdir()
    for v in sources:
        if v.suffix=='.py':shutil.copy2(v,a.work/'scripts'/v.name)
    planes=[]
    for i,entry in enumerate(plan['planes']):
        x,y,z=entry['XH'],entry['X3'],entry['Z']
        if not 0<x<1 or y<0 or z<=0 or x+y+z>=1:raise ValueError('invalid physical plane')
        seed=Path(entry['seed']);record=json.loads((seed/'validated.json').read_text())
        spec=dict(base);spec['metals']=[v*z/sum(base['metals']) for v in base['metals']]
        spec.update(hydrogen=[x],helium3=[y],teff_K=plan['teff_K'],log_g=plan['log_g'],
                    source='LTE gas atmospheres with specified GS98 metal mass; no grains')
        spec.pop('composition_control',None)
        # Revalidate the seed before spending time on a new opacity plane.
        initial=gen.completed_initial_structure(seed,spec,record['teff_K'],record['log_g'],spec['depths'])
        if initial is None:raise ValueError('seed exceeds target depth count')
        d=a.work/f'plane-{i:03d}';d.mkdir();sp=d/'specification.json';write_result(sp,spec)
        abundance,_=gen.composition(x,y,spec['metals']);T=gen.temperatures(spec);rho=gen.sequence(spec['log_density'])
        planes.append(dict(index=i,entry=entry,work=d,specification=spec,path=sp,abundance=abundance,temperature=T,density=rho))
        hashes[str(sp)]=digest(sp);hashes[str(seed/'validated.json')]=digest(seed/'validated.json')
    rp=root/'docs/research/fable/coordination/reservations'/plan['reservation']
    reservation=dict(owner='primary',state='running',threads=a.jobs,memory_bytes=6000000000,output_bytes=plan['maximum_output_bytes'],
        controller_pid=os.getpid(),started_utc=datetime.now(timezone.utc).isoformat(),scratch_paths=[str(a.work.resolve())],
        command=[sys.executable,str(Path(__file__).resolve()),str(a.plan.resolve()),str(a.work.resolve()),str(a.output.resolve()),'--jobs',str(a.jobs)],
        task_ids=['E-METAL-ATMOSPHERES-001'])
    write_result(rp,reservation);write_result(a.work/'plan.json',dict(plan=plan,input_sha256=hashes))
    def check_size():
        total=sum(v.stat().st_size for v in a.work.rglob('*') if v.is_file() and not v.is_symlink())
        if total>plan['maximum_output_bytes']:raise RuntimeError('batch output cap exceeded')
    def source(job):
        plane,j=job;e=plane['entry'];T=plane['temperature'][j]
        d=plane['work']/'opacity'/f'temperature-{j:03d}'
        gen.opacity_inputs(d,prepared,plane['specification'],e['XH'],e['X3'],T)
        gen.execute(prepared['synspec'],d,['fort.63','fort.29'])
        validate_table(d/'fort.63',plane['abundance'],[T],plane['density'])
        return dict(plane=plane['index'],isotherm=j,path=str(d/'fort.63'),sha256=digest(d/'fort.63'),seconds=json.loads((d/'completed.json').read_text())['seconds'])
    start=time.monotonic();completed=[];failures=[];states=[]
    try:
        jobs=[(plane,j) for j in reversed(range(len(planes[0]['temperature']))) for plane in planes]
        with ThreadPoolExecutor(a.jobs) as pool:
            futures={pool.submit(source,job):job for job in jobs}
            for future in as_completed(futures):
                plane,j=futures[future]
                try:
                    result=future.result();completed.append(result)
                    write_result(plane['work']/f'isotherm-{j:03d}.json',result)
                    print('opacity',len(completed),'/',len(jobs),flush=True);check_size()
                except Exception as e:failures.append(dict(phase='opacity',plane=plane['index'],isotherm=j,error=repr(e)))
        good=[]
        for plane in planes:
            if any(v['plane']==plane['index'] for v in failures):continue
            table=plane['work']/'opacity/fort.63'
            merge_isotherms([plane['work']/'opacity'/f'temperature-{j:03d}'/'fort.63' for j in range(len(plane['temperature']))],table)
            validate_table(table,plane['abundance'],plane['temperature'],plane['density'])
            write_result(plane['work']/'opacity_complete.json',dict(table=str(table),sha256=digest(table),entry=plane['entry']))
            good.append(plane)
        def column(job):
            plane,t,g=job;e=plane['entry'];d=plane['work']/f't{t:g}-g{g:g}'
            command=[sys.executable,str(root/'scripts/run_nongrey_continuation.py'),str(prepared_path),str(plane['path']),str(plane['work']/'opacity/fort.63'),e['seed'],str(d),
                '--initializer',plan['initializer'],'--hydrogen',str(e['XH']),'--helium3',str(e['X3']),'--teff',str(t),'--logg',str(g),
                '--initializer-depths','100','--initializer-frequencies','5000','--convective-iterations','10','--composition-continuation',
                '--initial-column-factor',str(plan['initial_column_factor'])]
            if plan.get('initial_bottom_tau') is not None:
                command += ['--initial-bottom-tau',str(plan['initial_bottom_tau'])]
            with (plane['work']/f't{t:g}-g{g:g}.log').open('x') as log:
                run=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,cwd=root,check=True)
            key,state,spec,record,deps=load_continuation(d)
            return dict(plane=plane['index'],coordinates=list(key),Z=e['Z'],work=str(d),state=state,input_sha256=deps)
        jobs=[(plane,t,g) for plane in good for t in plan['teff_K'] for g in plan['log_g']]
        with ThreadPoolExecutor(a.jobs) as pool:
            futures={pool.submit(column,job):job for job in jobs}
            for future in as_completed(futures):
                plane,t,g=futures[future]
                try:
                    result=future.result();states.append(result);write_result(plane['work']/f't{t:g}-g{g:g}.json',result)
                    print('atmosphere',len(states),'/',len(jobs),flush=True);check_size()
                except Exception as e:failures.append(dict(phase='atmosphere',plane=plane['index'],teff_K=t,logg=g,error=repr(e)))
    except Exception as e:failures.append(dict(phase='batch',error=repr(e)))
    for path,h in hashes.items():
        if digest(path)!=h:failures.append(dict(changed_input=path))
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='completed' if not failures else 'completed_with_failures',
        selected_for_evolution=False,scope=__doc__,plan=plan,input_sha256=hashes,isotherms=completed,states=states,failures=failures,
        elapsed_seconds=time.monotonic()-start,limitations=['Two H compositions and one reduced Z in the first batch do not define a complete variable-metal atmosphere grid.',
        'He3=0 sources isolate the metal response; any borrowed helium-isotope dependence needs a separate check.',
        'Fixed GS98 relative abundances, gas only; no independent selective metal depletion or condensation.'])
    write_result(a.output,report);reservation.update(state='released',ended_utc=datetime.now(timezone.utc).isoformat(),outcome=report['outcome'])
    rp.write_text(json.dumps(reservation,indent=2)+'\n');print(report['outcome'],len(states),'states',len(failures),'failures',flush=True)


if __name__=='__main__':main()
