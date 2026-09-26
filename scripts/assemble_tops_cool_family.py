#!/usr/bin/env python3
"""Assemble corrected cold and transition opacities, preserving hot inputs.

Combine actual native uncut means with plasma ratios evaluated through Ember.
Reuse grey means accompanying groups, completed grey refinements and the hot
source archive. Preserve original values and source exclusions; never fill
missing cells. Assembly is not acceptance of interpolation or stellar physics.
"""
import argparse
import json
import math
from pathlib import Path
import shlex
import shutil

from audit_refractive_opacity_family_v2 import read_table
from audit_tops_electron_dispersion import KEV, KB
from audit_tops_group_factor_tables import query
from assemble_tops_refractive_family_v2 import ratio_corners
from fetch_tops_composition import digest
from import_tops_composition import read
from reduce_tops_group_factors import add_inputs, verify


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('plan','hot_assembly','hot_family','probe','probe_reference','scratch','output','report'):
        p.add_argument('--'+name.replace('_','-'),required=True,type=Path)
    a=p.parse_args()
    if any(path.exists() for path in (a.scratch,a.output,a.report)):
        raise FileExistsError('preserve candidate assembly outputs')
    inputs={}
    def pin(path):
        path=Path(path);add_inputs(inputs,{str(path.resolve()):digest(path)});return path
    def report(path):
        r=json.loads(pin(path).read_text());add_inputs(inputs,r['input_sha256']);add_inputs(inputs,r.get('output_sha256',{}));return r
    spec=json.loads(pin(a.plan).read_text())
    hot=report(a.hot_assembly)
    report(a.probe_reference)
    if inputs.get(str(a.probe.resolve()))!=digest(a.probe):raise ValueError('unidentified C++ opacity probe')
    for name in ('assemble_tops_cool_family.py','audit_refractive_opacity_family_v2.py',
                 'audit_tops_group_factor_tables.py','assemble_tops_refractive_family_v2.py',
                 'audit_tops_factor_hydrogen.py','import_tops_composition.py','reduce_tops_group_factors.py'):
        pin(Path(__file__).with_name(name))
    ratios={}
    for plane in spec['ratio_planes']:
        if digest(Path(plane['plan']))!=plane['sha256']:raise ValueError('ratio plane specification changed')
        r=report(plane['reduction_report']);key=plane['X'],plane['Z']
        if (r['X'],r['Z'])!=key or r['temperatures_keV']!=spec['temperatures_keV']:
            raise ValueError('ratio plane differs')
        if any(n!=len(r['densities_atomic_g_cm3']) for n in r['density_prefix_sizes']):
            raise ValueError('incomplete requested ratio density rectangle')
        ratios[key]=r
    if set(ratios)!={(x,z) for x in spec['ratio_X'] for z in spec['ratio_Z']}:
        raise ValueError('missing ratio composition')
    tt=spec['temperatures_keV'];cold=set(spec['cold_temperatures_keV']);overlap=set(tt)-cold
    ratio_rr=next(iter(ratios.values()))['densities_atomic_g_cm3']
    if any(r['densities_atomic_g_cm3']!=ratio_rr for r in ratios.values()):raise ValueError('ratio density axes differ')
    grey_plan=json.loads(pin(spec['grey_baseline_plan']).read_text())
    mf=pin(grey_plan['baseline_manifest'])
    if digest(mf)!=grey_plan['baseline_manifest_sha256']:raise ValueError('native mixture manifest changed')
    native=json.loads(mf.read_text())['planes']
    composition={(r['X'],r['Z']) for r in native}
    xs=sorted({x for x,z in composition});zs=sorted({z for x,z in composition})
    if composition!={(x,z) for x in xs for z in zs}:raise ValueError('native compositions are not rectangular')
    native_rr=None
    for entry in native:
        path=pin(mf.parent/entry['file']);native_t,nr,_,_=read(path,entry)
        if any(t not in native_t for t in tt):raise ValueError('non-native temperature in assembly')
        nr=[r for r in nr if r<=ratio_rr[-1]]
        if native_rr is not None and native_rr!=nr:raise ValueError('native density axes differ')
        native_rr=nr
    rr=sorted(set(native_rr)|set(ratio_rr))
    required={(x,z,t,rho) for x,z in composition for t in tt for rho in rr}
    cells={};excluded=set();sources=[];overlaps=[];retained_logs={}
    def add(key,value,flag,origin):
        if key not in required:return
        if not math.isfinite(value) or value<=0:raise ValueError('invalid opacity source value')
        if key in cells:
            if flag!=(key in excluded):raise ValueError('source exclusion differs: '+str(key))
            error=abs(cells[key]/value-1)
            if not flag and error>1e-12:raise ValueError('overlapping grey sources differ: '+str(key))
            overlaps.append(error)
        else:
            cells[key]=value
            if flag:excluded.add(key)
    def grey_source(job,temperatures):
        if not set(job['temperatures_keV'])&temperatures:return
        directory=Path(job['work'])
        receipt=json.loads((directory/'receipt.json').read_text())
        request=json.loads((directory/'request.json').read_text())
        recipe=json.loads((directory/'recipe.json').read_text())
        if ((receipt['X'],receipt['Z'])!=(job['X'],job['Z']) or recipe['job']!=job or
            request['datype']!='gray' or request['plasnu']!='off' or request['lib']!='new' or
            digest(directory/'request.json')!=receipt['request_sha256']):
            raise ValueError('wrong grey scientific request')
        nt,nr,part,flags=read(directory/'source.txt',receipt,dimensions=tuple(receipt['dimensions']))
        if nt!=job['temperatures_keV'] or len(nr)!=len(job['densities_g_cm3']) or any(
            abs(r/q-1)>5e-5 for r,q in zip(nr,job['densities_g_cm3'],strict=True)):
            raise ValueError('grey coordinates differ from request')
        count=0
        for (t,rho),value in part.items():
            if t in temperatures:
                key=job['X'],job['Z'],t,rho
                if key in required:count+=1
                add(key,value,(t,rho) in flags,str(directory))
        if count:
            for name in ('receipt.json','request.json','recipe.json','source.txt'):pin(directory/name)
            add_inputs(inputs,recipe.get('input_sha256',{}))
            sources.append({'work':str(directory),'kind':'grey','states_in_required_rectangle':count})
    jobs=list(grey_plan['requests'])
    for item in grey_plan.get('retained_requests',[]):
        path=pin(item['plan'])
        if digest(path)!=item['plan_sha256']:raise ValueError('retained grey plan changed')
        jobs.append(json.loads(path.read_text())['requests'][item['request_index']])
    if len({j['work'] for j in jobs})!=len(jobs):raise ValueError('duplicate native grey source ownership')
    for job in jobs:grey_source(job,cold)
    original_native={key:value for key,value in cells.items() if key not in excluded}
    # The accepted hot assembly identifies all source batches, including its
    # completed pilot refinements. Read only batches touching the five overlap T.
    for item in hot['sources']:
        directory=Path(item['work']);recipe=json.loads((directory/'recipe.json').read_text())
        grey_source(recipe['job'],overlap)
    supplemental=json.loads(pin(spec['grey_refinement_plan']).read_text())
    if supplemental['baseline_manifest_sha256']!=digest(mf):raise ValueError('supplemental mixture source differs')
    for job in supplemental['requests']:grey_source(job,cold)
    for item in spec['retained_grey']:
        path=pin(item['report'])
        if digest(path)!=item['sha256']:raise ValueError('retained grey reduction changed')
        r=report(path);table_path=Path(next(iter(r['output_sha256'])));table=read_table(table_path)
        if table['X']!=[r['X']] or table['Z']!=r['Z']:raise ValueError('retained grey table composition differs')
        for it,t in enumerate(r['temperatures_keV']):
            for ir,logvalue in enumerate(table['planes'][r['X']][it]):
                rho=r['densities_atomic_g_cm3'][ir];key=r['X'],r['Z'],t,rho
                add(key,10**logvalue,False,str(table_path));retained_logs[key]=logvalue
    group_grey_count=0
    for (x,z),r in ratios.items():
        if (x,z) not in composition:continue
        for row in r['records']:
            key=x,z,row['temperature_keV'],row['density_atomic_g_cm3']
            if key in required:group_grey_count+=1
            add(key,row['uncut_source_rosseland'],False,row['source'])
    if required-set(cells):raise ValueError('missing grey source cells: '+str(sorted(required-set(cells))[:10]))
    if excluded:raise ValueError('source substitutions occur inside the requested cold rectangle: '+str(sorted(excluded)[:10]))
    if any(cells[k]!=v for k,v in original_native.items()):raise ValueError('original native mean changed')
    if any(abs(math.log10(cells[k])-v)>2e-15 for k,v in retained_logs.items()):
        raise ValueError('retained grey logarithm differs')
    verify(inputs)
    a.scratch.mkdir(parents=True)
    points=[(t,rho) for t in tt for rho in rr]
    ratio_at={}
    for key,r in ratios.items():
        table=next(Path(p) for p in r['output_sha256'] if Path(p).name=='factor.dat')
        values=query(a.probe,table,points,key[0],a.scratch,f'x{key[0]:g}-z{key[1]:g}')
        if any(v[0]<1-1e-12 for v in values):raise ValueError('interpolated plasma ratio below unity')
        ratio_at[key]={point:v[0] for point,v in zip(points,values,strict=True)}
    corners={key:ratio_corners(*key,spec['ratio_X'],spec['ratio_Z']) for key in composition}
    a.output.mkdir(parents=True)
    copies=[]
    for filename in ('aesopus21_gs98_mixture.dat','tops_gs98_mixture_high.dat'):
        manifest_path=pin(a.hot_family/filename);lines=manifest_path.read_text().splitlines()
        fields=lines[0].split()
        if fields[:2]!=['EMBER_OPACITY_MIXTURE','1'] or int(fields[2])!=len(lines)-1:
            raise ValueError('invalid preserved opacity manifest')
        rewritten=[lines[0]]
        for line in lines[1:]:
            z,name=shlex.split(line);source=(manifest_path.parent/name).resolve(strict=True)
            if inputs.get(str(source))!=digest(source):raise ValueError('unidentified preserved opacity table')
            target=a.output/source.name;shutil.copyfile(source,target)
            if digest(source)!=digest(target):raise ValueError('preserved table changed')
            rewritten.append(f'{z} "{source.name}"');copies.append(str(target.resolve()))
        (a.output/filename).write_text('\n'.join(rewritten)+'\n')
    manifest=[f'EMBER_OPACITY_MIXTURE 1 {len(zs)} logRho ATOMIC_cold_native_uncut_times_refractive_ratio']
    for z in zs:
        lines=[f'EMBER_OPACITY_TABLE 2 {len(xs)} {len(tt)} {len(rr)} Native uncut grey times plasma ratio; complete source rectangle',
               ' '.join(format(math.log10(r),'.17g') for r in rr),
               ' '.join(format(math.log10(t*KEV/KB),'.17g') for t in tt)]
        for x in xs:
            lines.append(f'{x:.17g} {z:.17g}')
            for t in tt:
                logs=[]
                for rho in rr:
                    lf=sum(weight*math.log(ratio_at[corner][t,rho]) for corner,weight in corners[x,z])
                    value=(math.log(cells[x,z,t,rho])+lf)/math.log(10)
                    if not math.isfinite(value):raise ValueError('nonfinite assembled opacity')
                    logs.append(value)
                lines.append(str(len(rr))+' '+' '.join(format(v,'.17g') for v in logs))
        name=f'tops_gs98_mixture_z{round(z*1000):03d}_low.dat'
        (a.output/name).write_text('\n'.join(lines)+'\n');manifest.append(f'{z:.17g} "{name}"')
    (a.output/'tops_gs98_mixture_low.dat').write_text('\n'.join(manifest)+'\n')
    verify(inputs)
    result={'scope':__doc__,'accepted_for_stellar_opacity':False,'complete_native_composition_rectangle':True,
        'X':xs,'Z':zs,'temperature_keV':tt,'density_atomic_g_cm3':rr,'cold_and_transition_states_written':len(required),
        'native_cold_grey_states_retained':len(original_native),'group_accompanying_grey_states_used':group_grey_count,
        'maximum_overlapping_grey_difference':max(overlaps,default=0.),'overlapping_grey_comparisons':len(overlaps),
        'unique_ratio_runtime_queries':len(ratio_at)*len(points),'preserved_opacity_files':copies,
        'source_exclusions':[],'sources':sources,'input_sha256':inputs,
        'output_sha256':{str(p.resolve()):digest(p) for p in a.output.iterdir() if p.is_file()},
        'scratch_sha256':{p.name:digest(p) for p in a.scratch.iterdir() if p.is_file()},
        'limitations':['Assembly checks source completeness and retention, not independent interpolation accuracy.',
                       'Temperature joins, isotope mapping, independent spectra and actual stellar queries remain to be checked.',
                       'Dense cold conduction and absolute plasma-opacity physics are not validated by this table assembly.']}
    a.report.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('cold_and_transition_states_written','native_cold_grey_states_retained',
        'group_accompanying_grey_states_used','maximum_overlapping_grey_difference','unique_ratio_runtime_queries')}),flush=True)


if __name__=='__main__':main()
