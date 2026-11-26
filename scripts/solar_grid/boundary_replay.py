"""Paired native experiments with a withheld measured atmosphere knot."""
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

from release import thermal_turns
from cold_restart_check import audited_steps, same_restart
from cooling_interval import BASE, LIMITS, probe, crossing, digest, read, rows, save


def inputs(root, study):
    band=root/'m06/cooling-intervals/band-1'
    if study=='m05':
        case=root/'warm-marcs/extended-start-m05/native-continuation/below-1100K'
        return case,case,1050,1100,1000,None
    warmer=(root/BASE).parent
    if study=='m06':
        return warmer.parent/'direct-envelope-opacity-check',warmer,900,1000,850,None
    if study=='early850':
        return warmer,warmer,850,900,850,root/'atmosphere-validation/early850/fixed-profile.json'
    d=read(band/'decision.json')
    return band/'native',band/'native',d['midpoint_K'],d['edge_K'],d['target_K'],None


def change_family(family, directory, knot, scale=None):
    directory.mkdir(parents=True,exist_ok=True)
    manifest=[]
    for line in family.read_text().splitlines():
        if line.startswith('file '):
            src=Path(shlex.split(line)[1])
            data=src.read_text().splitlines()
            axis=data[4].split()[2:]
            j=next(i for i,v in enumerate(axis) if math.isclose(10**float(v),knot,abs_tol=1e-6))
            ng=int(data[5].split()[1])
            off=7+j*ng
            assert len(data)>=7+ng*len(axis)
            if scale is None:
                assert 0<j<len(axis)-1
                del data[off:off+ng]
                del axis[j]
                data[4]=' '.join(['log_teff',str(len(axis)),*axis])
            else:
                for i in range(off,off+ng):
                    t,p,h=map(float,data[i].split())
                    data[i]=f'{t+math.log10(scale[0]):.17g} {p+math.log10(scale[1]):.17g} {h*scale[2]:.17g}'
            member=directory/src.name
            member.write_text('\n'.join(data)+'\n')
            line=f'file "{member}"'
        manifest.append(line)
    output=directory/'family.dat'
    output.write_text('\n'.join(manifest)+'\n')
    return output


def prepare(root,study):
    from atmosphere_repair import members
    case,table_case,knot,warm,cold,early=inputs(root,study)
    work=root/'atmosphere-validation'/study
    cfg=(case/'config.txt').read_text()
    old=Path(dict(shlex.split(l) for l in cfg.splitlines() if l.strip())['atmosphere_pressure_matched'])
    measured=Path(dict(shlex.split(l) for l in (table_case/'config.txt').read_text().splitlines() if l.strip())['atmosphere_pressure_matched'])
    column_hash=None
    if study in ('m06','early850'):
        local=root/f'cool-m06-{knot}-birth-depth10'
        spec=read(local/'spec.json')
        receipt=Path(read(root/'recovery'/f'm06-{knot}-birth-depth10-column'/'selected-column.json')['receipt'])
        native=read(receipt)
        assert native['quality']['usable'] and native['teff_used']==knot
        sys.path.insert(0,str(root/'thermal-holdout/m06'))
        from archive_matching_points import matching_points
        pt=matching_points(native)['10']
        truth=[pt['T'],pt['P_bar']*1e6,pt['depth_cm_mu_He3corr']]
        predicted=probe(root,measured,knot,spec)
        measured=change_family(measured,work/'measured-table',knot,[v/p for v,p in zip(truth,predicted,strict=True)])
        assert max(abs(v/p-1) for v,p in zip(truth,probe(root,measured,knot,spec),strict=True))<1e-8
        column_hash=digest(receipt)
    scale=None if early is None else read(early)['early_over_final']
    changed=change_family(measured,work/'variant-table',knot,scale)
    ck=case/'start.checkpoint'
    seed=ck.read_bytes().splitlines(keepends=True)
    off=9+int(seed[8])
    ids=dict(shlex.split(line.decode()) for line in seed[9:off])
    for k,v in members(old).items():
        assert ids.pop(k)==v
    initial_dir=case/('segment-0' if study=='band1' else 'evolution')
    initial=rows(initial_dir/'history.jsonl')[0]
    for arm,family in [('measured',measured),('trial',changed)]:
        directory=work/arm
        directory.mkdir(exist_ok=True)
        assert cfg.count(str(old))==1
        (directory/'config.txt').write_text(cfg.replace(str(old),str(family)))
        keys=dict(ids,**members(family))
        import json
        header=''.join(f'{json.dumps(k)} {json.dumps(v)}\n' for k,v in sorted(keys.items()))
        start=directory/'start.checkpoint'
        start.write_bytes(b''.join(seed[:8])+f'{len(keys)}\n'.encode()+header.encode()+b''.join(seed[off:]))
        shutil.copyfile(str(ck)+'.lithium',str(start)+'.lithium')
    save(work/'setup.json',dict(knot=knot,warm=warm,cold=cold,initial=initial,source=str(ck),
        source_sha256=digest(ck),lithium_sha256=digest(Path(str(ck)+'.lithium')),
        column_sha256=column_hash,measured_members=members(measured),trial_members=members(changed),
        scope='Identical initial interior, fuel, microphysics and native solver. Only atmosphere knot changes.'))


def run(root,study,arm):
    from atmosphere_repair import members
    from coordinate import solver
    work=root/'atmosphere-validation'/study
    d=read(work/'setup.json')
    directory=work/arm
    seed=directory/'start.checkpoint'
    family=Path(dict(shlex.split(l) for l in (directory/'config.txt').read_text().splitlines() if l.strip())['atmosphere_pressure_matched'])
    assert members(family)==d[arm+'_members']
    source=Path(d['source'])
    assert digest(source)==d['source_sha256'] and digest(Path(str(seed)+'.lithium'))==d['lithium_sha256']
    def material(p):
        lines=p.read_bytes().splitlines(keepends=True)
        return lines[:8],lines[9+int(lines[8]):]
    assert material(seed)==material(source)
    previous=d['initial']
    segments=[]
    # Compare strictly inside the table; no need to converge to its edge.
    end=d['cold']+1
    for n in range(5):
        seg=directory/f'segment-{n}'
        cmd=[str(solver),'--lifetime',str(directory/'config.txt'),str(5*previous['years']),str(seg),'--restart',str(seed),'--cpu-seconds','900']
        if n==0:cmd.extend(['--restart-step-years','10000'])
        if not (seg/'report.json').exists():
            with (directory/f'segment-{n}.log').open('w') as log:
                subprocess.run(cmd,stdout=log,stderr=log,check=False)
        report=read(seg/'report.json')
        states=rows(seg/'history.jsonl')
        assert same_restart(previous,states[0]) and audited_steps(states,rows(seg/'attempts.jsonl'),report)
        segments.append(str(seg))
        if min(r['Teff_K'] for r in states)<=end:break
        if report['stop_reason'].strip().lower()!='cpu budget':
            raise ValueError(report['stop_reason'])
        previous=states[-1]
        seed=seg/'final.checkpoint'
    else:raise ValueError('Replay CPU budgets exhausted')
    save(directory/'finished.json',dict(segments=segments,source_sha256=d['source_sha256'],
         setup_sha256=digest(work/'setup.json'),hashes={str(Path(s)/f):digest(Path(s)/f) for s in segments for f in ('history.jsonl','attempts.jsonl','report.json')}))


def compare_series(a,b,warm,knot,cold):
    start=warm-1
    fields=('years','central_T_K','radius_Rsun','nuclear_fraction')
    result={}
    def times(series,temp):
        for l,r in zip(series,series[1:]):
            if l['Teff_K']>=temp>=r['Teff_K'] and l['Teff_K']>r['Teff_K']:
                if l['Teff_K']==temp:return l['years'],l['years']
                if r['Teff_K']==temp:return r['years'],r['years']
                return l['years'],r['years']
        raise ValueError(f'Missing crossing {temp}')
    for target in sorted({knot,cold+1},reverse=True):
        if target<cold+1:continue
        va,vb=[crossing(s,target,fields) for s in (a,b)]
        dt=[v[0]-crossing(s,start,fields)[0] for s,v in ((a,va),(b,vb))]
        err=dict(elapsed_relative=dt[1]/dt[0]-1,core_relative=vb[1]/va[1]-1,
                 radius_relative=vb[2]/va[2]-1,nuclear_fraction_difference=vb[3]-va[3])
        bounds=[]
        for s in (a,b):
            lo,hi=times(s,start);x,y=times(s,target)
            bounds.append((x-hi,y-lo))
        (al,ah),(bl,bh)=bounds
        interval=[bl/ah-1,bh/al-1] if min(al,bl)>0 else None
        result[str(target)]=dict(elapsed_years=dt,errors=err,elapsed_bracket=interval)
    return result


def score(root,study):
    work=root/'atmosphere-validation'/study
    d=read(work/'setup.json')
    series=[]
    for arm in ('measured','trial'):
        receipt=read(work/arm/'finished.json')
        assert digest(work/'setup.json')==receipt['setup_sha256']
        assert all(digest(Path(p))==v for p,v in receipt['hashes'].items())
        samples=[]
        for s in receipt['segments']:
            block=rows(Path(s)/'history.jsonl')
            samples.extend(block[1:] if samples else block)
        series.append(samples)
    scores=compare_series(*series,d['warm'],d['knot'],d['cold'])
    # An interval with a resolved temperature reversal is ineligible for the smooth cooling shortcut.
    reversals=[thermal_turns(samples) for samples in series]
    ok=not any(reversals) and all(all(abs(v)<=LIMITS[k] for k,v in point['errors'].items()) and
         point['elapsed_bracket'] is not None and max(map(abs,point['elapsed_bracket']))<=LIMITS['elapsed_relative'] for point in scores.values())
    calibration={}
    if study=='early850':
        fixed=read(work/'fixed-profile.json')
        calibration=dict(fixed_profile_sha256=digest(work/'fixed-profile.json'),native_clock_passes=ok,
            physical_checks_passed=fixed['physical_checks_passed'],matching_within_limits=fixed['matching_within_limits'])
        ok=ok and fixed['physical_checks_passed'] and fixed['matching_within_limits']
    save(work/'comparison.json',dict(calibration=calibration,passes=ok,scores=scores,limits=LIMITS,reversals=reversals,warm=d['warm'],cold=d['cold'],
         setup_sha256=digest(work/'setup.json'),scope='Local native interval; a change in reversal behavior excludes interpolation.'))


if __name__=='__main__':
    root=Path(sys.argv[1]).resolve()
    os.environ['EMBER_GRID_ROOT']=str(root)
    study,action=sys.argv[2:]
    if action=='prepare':prepare(root,study)
    elif action=='score':score(root,study)
    else:run(root,study,action)
