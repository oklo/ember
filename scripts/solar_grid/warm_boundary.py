"""Measure the evolved mixture with warm opacities at the stellar matching depth."""
import gzip
import hashlib
import itertools
import tarfile
import math
from pathlib import Path
import subprocess
import sys

from release import digest, read, save
from cooling_interval import MATCH_LIMITS, probe


def paths(root, mass):
    case = root/f'warm-marcs/extended-start-m{mass:02}/hydrogen-edge-column'
    return case, case/'warm-source'


def triage(root,mass):
    case,work=paths(root,mass)
    spec=read(case/'spec.json')
    t,g=spec['teff'],spec['logg']
    seed=Path(spec['seed_receipt'])
    if not seed.is_absolute():seed=root.parent/'picaso-pilot-oct11-fable'/seed
    limit=read(seed)['fresh']['coverage']['CK_T'][1]
    out=dict(spec_sha256=digest(case/'spec.json'),temperature_ceiling_K=limit,warm=False)
    archive=root.parent/'lowmass-metallicity-oct10-fable/sources/marcs/marcs_st_mod.tar'
    # Outside this archive a missing estimate is recorded explicitly.
    if 2500<=t<4000 and 0<=g<5.5:
        low_t=100*math.floor(t/100);low_g=math.floor(g*2)/2
        points={}; hashes={};ratio={}
        with tarfile.open(archive) as ar:
            for dt,dg,dz in itertools.product((0,1),repeat=3):
                name=f'./p{low_t+dt*100}_g{low_g+dg*.5:+.1f}_m0.0_t02_st_z{dz*.25:+.2f}_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod.gz'
                try:raw=ar.extractfile(name).read()
                except KeyError:continue
                hashes[name]=hashlib.sha256(raw).hexdigest()
                lines=gzip.decompress(raw).decode().splitlines()
                x,_,z=map(float,next(l for l in lines if 'are X, Y and Z' in l).split()[:3])
                ratio[dz]=z/x
                i=next(i for i,l in enumerate(lines) if l.split()[:3]==['k','lgTauR','lgTau5'])
                size=int(next(l for l in lines if 'Number of depth points' in l).split()[0])
                rows=[list(map(float,l.split())) for l in lines[i+1:i+1+size]]
                vals=[]
                for tau,field in [(3,4),(3,6),(30,4)]:
                    a,b=next((a,b) for a,b in zip(rows,rows[1:]) if a[1]<=math.log10(tau)<=b[1])
                    w=(math.log10(tau)-a[1])/(b[1]-a[1])
                    vals.append(math.log(a[field])+w*math.log(b[field]/a[field]))
                points[dt,dg,dz]=vals
        baselines={}
        for label,metal in [('stellar',spec['Z_over_X']),('reference',.02/.7)]:
            if len(points)!=8 or not ratio[0]<=metal<=ratio[1]:continue
            weights=((t-low_t)/100,(g-low_g)/.5,math.log(metal/ratio[0])/math.log(ratio[1]/ratio[0]))
            vals=[math.exp(sum(v[k]*math.prod(w if bit else 1-w for w,bit in zip(weights,key))
                               for key,v in points.items())) for k in range(3)]
            baselines[label]=dict(Z_over_X=metal,MARCS_T_and_Pgas=vals[:2],bottom_K=vals[2])
        if len(baselines)==2:
            bottom=baselines['stellar']['bottom_K']
            out.update(warm=bottom>limit,baseline=baselines,estimated_bottom_K=bottom,source_members=hashes)
    save(work/'triage.json',out)


def prepare(root, mass):
    case, work = paths(root, mass)
    spec = read(case/'spec.json')
    triage_file=work/'triage.json'
    pre=read(triage_file)
    assert pre['warm'] and pre['spec_sha256']==digest(case/'spec.json')
    template=read(root.parent/'warm-reference-atmospheres-oct2-v1/h700_y000_t4400_g515/specification.json')
    template.update(tau=3,teff_K=[spec['teff']],log_g=[spec['logg']],newton_relaxation=.3)
    save(work/'reference-spec.json',template)
    mixture = spec['arms']['anchor']
    assert math.isclose(sum(mixture.values()),1,abs_tol=1e-8)
    assert math.isclose(pre['baseline']['stellar']['Z_over_X'],mixture['Z']/mixture['X'])
    target = dict(template,hydrogen=[mixture['X']],helium3=[mixture['He3']],
                  metals=[m*mixture['Z']/sum(template['metals']) for m in template['metals']])
    assert target['tau']==3 and target['teff_K']==[spec['teff']] and target['log_g']==[spec['logg']]
    save(work/'target-spec.json',target)
    save(work/'route.json',dict(source='MARCS baseline and matched-depth TLUSTY/SYNSPEC differential',
         cold_opacity_max_K=pre['temperature_ceiling_K'],estimated_bottom_K=pre['estimated_bottom_K'],
         baseline=pre['baseline'],tau=3,mixture=mixture,
         input_hashes={str(p):digest(p) for p in [case/'spec.json',triage_file, work/'reference-spec.json',case/'candidate-family.dat',work/'target-spec.json']}))


def opacity_plane(root,work,prepared):
    from generate_nongrey_grid import composition, execute, opacity_inputs, temperatures, sequence
    from prepare_nongrey_sources import data_digest, digest as stream_digest
    # The retained warm source writes double-precision binary SYNSPEC tables.
    sys.path.insert(0,str(root.parent/'dense-capacity-sept29-v1'))
    from nongrey_opacity_r64 import merge_isotherms,validate_table
    spec=read(work/'target-spec.json')
    x,y=spec['hydrogen'][0],spec['helium3'][0]
    abundance,_=composition(x,y,spec['metals'])
    exe=Path(prepared['synspec'])
    assert stream_digest(exe)==prepared['executables']['synspec']
    assert all(stream_digest(p)==h for p,h in zip(prepared['line_lists'],prepared['line_list_sha256'],strict=True))
    assert data_digest(Path(prepared['synple'])/'data')==prepared['data_sha256']
    temps=temperatures(spec)
    densities=sequence(spec['log_density'])
    base=work/'opacity-grid'
    paths=[base/f'temperature-{i:03d}'/'fort.63' for i in range(len(temps))]
    middle=len(temps)//2
    for i in [middle,*[j for j in range(len(temps)) if j!=middle]]:
        folder=paths[i].parent
        opacity_inputs(folder,prepared,spec,x,y,temps[i])
        execute(exe,folder,['fort.63'])
        validate_table(paths[i],abundance,[temps[i]],densities)
    table=base/'opacity.bin'
    temp=table.with_suffix('.tmp')
    merge_isotherms(paths,temp)
    result=validate_table(temp,abundance,temps,densities)
    temp.replace(table)
    save(work/'opacity.json',dict(table=str(table),sha256=digest(table),spec_sha256=digest(work/'target-spec.json'),
         shape=result['shape'],isotherm_sha256=[digest(p) for p in paths]))


def solve(root, mass, action):
    from generate_nongrey_grid import atmosphere_inputs, composition, continuation_structure
    case, work = paths(root,mass)
    recipe_dir = root.parent/'production-source-oct2-v1/scripts'
    archive = root.parent/'warm-reference-atmospheres-oct2-v1'
    donor = archive/'h700_y000_t4400_g515'
    prepared_file = donor/'provenance.json'
    prepared = read(prepared_file)
    route = read(work/'route.json')
    assert all(digest(Path(p))==v for p,v in route['input_hashes'].items())
    if action=='opacity':
        opacity_plane(root,work,prepared)
        return
    target = action=='target'
    if target:
        spec = read(work/'target-spec.json')
        receipt = read(work/'opacity.json')
        table = Path(receipt['table'])
        assert digest(table)==receipt['sha256'] and receipt['spec_sha256']==digest(work/'target-spec.json')
    else:
        spec = read(work/'reference-spec.json')
        table=archive/'h700_y000/opacity.bin'
        assert digest(table)==read(table.parent/'READY.json')['sha256']
    assert digest(Path(prepared['tlusty']))==prepared['executables']['tlusty']
    record=read(donor/'result.json')
    assert record['accepted'] and digest(donor/'fort.7')==record['retained_sha256']['fort.7']
    col=work/action
    col.mkdir(exist_ok=True)
    if not (col/'result.json').exists():
        seed=donor/'fort.7'
        if target and (work/'reference/boundary.json').is_file() and read(work/'reference/result.json')['accepted']:
            seed=work/'reference/fort.7'
            assert digest(seed)==read(seed.parent/'result.json')['retained_sha256']['fort.7']
        temp,gravity=spec['teff_K'][0],spec['log_g'][0]
        initial=continuation_structure(seed.read_text(),4400 if seed.parent==donor else temp,
                                       5.15 if seed.parent==donor else gravity,temp,gravity)
        abund,masses=composition(spec['hydrogen'][0],spec['helium3'][0],spec['metals'])
        atmosphere_inputs(col,prepared,spec,table,abund,masses,temp,gravity,initial)
        save(col/'specification.json',spec)
        save(col/'provenance.json',dict(prepared,initial_guess_sha256=digest(seed)))
        subprocess.run([sys.executable,recipe_dir/'solve_column.py',col,'--executable',prepared['tlusty'],
                        '--executable-sha256',prepared['executables']['tlusty'],'--attempt-cpu','1800',
                        '--total-cpu','5400','--strategy','hybrid','--chmax','1e-6'],check=False)
    result=read(col/'result.json')
    if not result['accepted']:
        raise ValueError(f'{action} warm solve failed its independent quality checks: {result["attempts"][-1].get("failure")}')
    log=next(p for p in (col/'run.log',col/'run.log.gz') if p.is_file())
    assert digest(log)==result['retained_sha256'][log.name]
    source=gzip.decompress(log.read_bytes()).decode() if log.suffix=='.gz' else log.read_text()
    state=result['diagnostics']
    save(col/'boundary.json',dict(T=state['T'],P=state['Pgas'],depth=thickness(source),
        tau=spec['tau'],opacity_sha256=digest(table),receipt_sha256=digest(col/'result.json')))


def thickness(log):
    # Integrate dz = dm/rho with a local power law in column mass across each layer.
    layers=[]
    for line in log.rsplit('FINAL MODEL ATMOSPHERE',1)[1].splitlines():
        fields=line.replace('D','E').split()
        if len(fields)==11 and fields[0].isdigit():
            m,tau,rho=map(float,(fields[1],fields[2],fields[5]))
            layers.append((m,tau,rho))
    assert layers[0][1]<2/3<3<layers[-1][1]
    total=0.
    for (m1,t1,r1),(m2,t2,r2) in zip(layers,layers[1:]):
        assert m2>m1>0 and t2>t1>0 and min(r1,r2)>0
        lo,hi=max(t1,2/3),min(t2,3)
        if lo>=hi:continue
        logm=math.log(m2/m1)
        a=logm*math.log(lo/t1)/math.log(t2/t1)
        b=logm*math.log(hi/t1)/math.log(t2/t1)
        exponent=1-math.log(r2/r1)/logm
        total+=m1/r1*math.exp(exponent*a)*(math.expm1(exponent*(b-a))/exponent if exponent else b-a)
    assert math.isfinite(total) and total>0
    return total


def assess(root,mass):
    case,work=paths(root,mass)
    route=read(work/'route.json')
    assert all(digest(Path(k))==v for k,v in route['input_hashes'].items())
    ref,star=[read(work/k/'boundary.json') for k in ('reference','target')]
    for n,d in [('reference',ref),('target',star)]:
        assert d['receipt_sha256']==digest(work/n/'result.json') and d['tau']==3
    target=read(case/'spec.json')
    predicted=probe(root,case/'candidate-family.dat',target['teff'],target)
    assert all(math.isclose(a,b,rel_tol=1e-10) for a,b in zip(predicted,read(case/'prediction.json'),strict=True))
    measured=[star[k] for k in ('T','P','depth')]
    reference=route['baseline']['reference']['MARCS_T_and_Pgas']
    corrected=[v*star[k]/ref[k] for v,k in zip(reference,('T','P'),strict=True)]
    discrepancies=[a/b-1 for a,b in zip(predicted,measured,strict=True)]
    cross=[a/b-1 for a,b in zip(corrected,measured[:2],strict=True)]
    ok=all(abs(e)<=l for e,l in zip(discrepancies,MATCH_LIMITS,strict=True)) and all(abs(e)<=l for e,l in zip(cross,MATCH_LIMITS[:2],strict=True))
    source=work/'target/boundary.json'
    save(case/'assessment.json',dict(passes=ok,source=str(source),receipt_sha256=digest(source),native=measured,
         predicted=predicted,relative=discrepancies,limits=MATCH_LIMITS,marcs_corrected_T_P=corrected,
         marcs_cross_source_relative=cross,reference_sha256=digest(work/'reference/boundary.json'),
         scope='Warm source at actual X, He3 and Z; GS98 metal proxy and effective helium mass.'))


if __name__=='__main__':
    root=Path(sys.argv[1]).resolve()
    mass=int(sys.argv[2])
    action=sys.argv[3]
    sys.path.insert(0,str(root.parent/'production-source-oct2-v1/scripts'))
    if action=='triage':triage(root,mass)
    elif action=='prepare':prepare(root,mass)
    elif action=='assess':assess(root,mass)
    else:solve(root,mass,action)
