#!/usr/bin/env python3
"""Compute canonical depleted-atmosphere chains at fixed composition/gravity.

Opacity planes are supplied explicitly and must already be complete. A plan
may reference previously validated cells; remaining cells receive full source
solves and independent chemistry audits. Fifty-kelvin continuation resolves
cold condensation fronts; hotter cells start from matching gas atmospheres.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import itertools
import json
from pathlib import Path
import subprocess
import sys
from validate_condensate_model import validate
from prepare_nongrey_sources import digest

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['prepared','specification','plan','work']:p.add_argument(name,type=Path)
    p.add_argument('--jobs',type=int,default=4)
    p.add_argument('--plane',type=int,nargs='+',help='selected composition planes; separate work directories for concurrent subsets')
    p.add_argument('--gravity',type=float,nargs='+')
    p.add_argument('--minimum-teff',type=int,choices=[2600,2650,2700,2750,2800,3000],default=2600,
                   help='explicit partial chain when colder source states are unsupported; final family still requires a complete rectangle')
    p.add_argument('--convective-initializer',type=Path,
                   help='optional verified CONREF source receipt for intermediate cold-temperature guesses')
    p.add_argument('--convective-iterations',type=int,default=20,
                   help='initial CONREF iterations before ordinary equations resume in the auxiliary solve; canonical replay remains mandatory')
    p.add_argument('--initialize-first',action='store_true',
                   help='also precondition the first cold model when beginning a partial chain from a gas atmosphere')
    p.add_argument('--initial-model',type=Path,
                   help='validated colder condensate seed for a single composition/gravity chain')
    a=p.parse_args();spec=json.loads(a.specification.read_text());plan=json.loads(a.plan.read_text())
    prepared=json.loads(a.prepared.read_text())
    initializer=a.convective_initializer or (Path(plan['convective_initializer']) if plan.get('convective_initializer') else None)
    if initializer:
        init=json.loads(initializer.read_text());method=init.get('initialization_only',{})
        if (method.get('canonical_prepared')!=prepared or method.get('method')!='native CONREF with explicit bounds repairs' or
                digest(init['tlusty'])!=init['executables']['tlusty'] or
                digest(Path(prepared['tlusty']).parent/'tlusty208.f')!=method['base_source_sha256']):
            raise ValueError('convective initializer does not match canonical source')
    elif a.initialize_first:raise ValueError('initializing the first cell requires a verified initializer')
    if not 1<=a.jobs<=12:raise ValueError('invalid worker count')
    if not 1<=a.convective_iterations<=200:raise ValueError('invalid convective initializer iterations')
    if spec['teff_K']!=[2600,2800,3000,3200]:raise ValueError('temperature chain needs explicit adaptation for different axes')
    planes=list(itertools.product(spec['hydrogen'],spec['helium3']))
    selected=a.plane if a.plane is not None else list(range(len(planes)))
    gravity=a.gravity if a.gravity is not None else spec['log_g']
    if len(set(selected))!=len(selected) or any(n<0 or n>=len(planes) for n in selected):raise ValueError('invalid planes')
    if len(set(gravity))!=len(gravity) or any(g not in spec['log_g'] for g in gravity):raise ValueError('invalid gravities')
    if a.initial_model:
        if len(selected)!=1 or len(gravity)!=1:raise ValueError('a supplied seed requires exactly one chain')
        seed=validate(a.initial_model,spec,prepared)
        x,y,t,g=seed['coordinates']
        if (x,y)!=planes[selected[0]] or g!=gravity[0] or t>=a.minimum_teff:
            raise ValueError('supplied seed must be colder and have the same composition and gravity')
    a.work.mkdir(parents=True,exist_ok=True)
    inputs={'prepared':prepared,'specification':spec,'plan':plan,'planes':selected,'gravity':gravity}
    if a.minimum_teff!=2600:inputs['minimum_teff']=a.minimum_teff
    if a.initialize_first:inputs['initialize_first']=True
    if initializer:inputs['convective_initializer']=init
    if a.initial_model:inputs['initial_model']={'directory':str(a.initial_model.resolve()),
                                              'receipt_sha256':digest(a.initial_model/'completed.json')}
    receipt=a.work/'provenance.json'
    if receipt.exists() and json.loads(receipt.read_text())!=inputs:raise ValueError('grid inputs changed')
    receipt.write_text(json.dumps(inputs,indent=2)+'\n')
    def known(x,y,t,g):
        records=[r for r in plan.get('models',[]) if [r[k] for k in ['XH','X3','teff_K','log_g']]==[x,y,t,g]]
        if len(records)>1:raise ValueError('duplicate precomputed model')
        return Path(records[0]['directory']) if records else None
    def chain(job):
        n,g=job;x,y=planes[n];ig=spec['log_g'].index(g)
        table_records=[r for r in plan['opacity_planes'] if (r['XH'],r['X3'])==(x,y)]
        if len(table_records)!=1:raise ValueError('missing opacity plane')
        table=Path(table_records[0]['directory'])/'opacity.bin'
        if not table.exists():raise ValueError('opacity plane is not complete')
        previous=a.initial_model;models=[];continuation_models=[]
        for t in [t for t in [2600,2650,2700,2750,2800,3000,3200] if t>=a.minimum_teff]:
            target=known(x,y,t,g)
            if target is None:
                target=a.work/f'plane-{n:03d}'/f'teff-{t}-gravity-{ig:03d}'
                if not (target/'validated.json').exists():
                    if previous is None or t>=3000:
                        gas_spec=json.loads((Path(plan['gas_grid'])/'specification.json').read_text())
                        it=min(range(len(gas_spec['teff_K'])),key=lambda i:abs(gas_spec['teff_K'][i]-t))
                        if (gas_spec['hydrogen']!=spec['hydrogen'] or gas_spec['helium3']!=spec['helium3']
                                or gas_spec['log_g']!=spec['log_g']):raise ValueError('gas initializer axes differ')
                        initial=Path(plan['gas_grid'])/f'plane-{n:03d}'/f'model-{it:03d}-{ig:03d}'
                    else:initial=previous
                    if initializer and (t in [2750,2800] or (a.initialize_first and previous is None and t<3000)):
                        precondition=target.with_name(target.name+'-initializer')
                        command=[sys.executable,str(ROOT/'scripts/run_condensate_atmosphere.py'),
                            str(initializer.resolve()),str(a.specification.resolve()),str(table.resolve()),
                            str(initial.resolve()),str(precondition.resolve()),'--hydrogen',str(x),'--helium3',str(y),
                            '--teff',str(t),'--logg',str(g),'--mode','equilibrium','--convective-iterations',str(a.convective_iterations)]
                        precondition.parent.mkdir(parents=True,exist_ok=True)
                        with Path(str(precondition)+'.log').open('w') as log:
                            subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
                        initial=precondition
                    cmd=[sys.executable,str(ROOT/'scripts/run_condensate_atmosphere.py'),
                        str(a.prepared.resolve()),str(a.specification.resolve()),str(table.resolve()),
                        str(initial.resolve()),str(target.resolve()),'--hydrogen',str(x),'--helium3',str(y),
                        '--teff',str(t),'--logg',str(g),'--mode','equilibrium']
                    target.parent.mkdir(parents=True,exist_ok=True)
                    with Path(str(target)+'.log').open('w') as log:
                        subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
            if not (target/'chemistry-audit.json').exists():
                with Path(str(target)+'-audit.log').open('w') as log:
                    subprocess.run([sys.executable,str(ROOT/'scripts/audit_condensate_atmosphere.py'),str(target)],
                        stdout=log,stderr=subprocess.STDOUT,check=True)
            result=validate(target,spec,prepared,[x,y,t,g])
            previous=target
            record={'XH':x,'X3':y,'teff_K':t,'log_g':g,'directory':str(target.resolve()),
                    'diagnostics':result['diagnostics'],'chemistry':result['chemistry']}
            continuation_models.append(record)
            if t in spec['teff_K']:models.append(record)
            print('Validated',x,y,t,g,flush=True)
        path=a.work/f'chain-{n:03d}-{ig:03d}.json'
        path.write_text(json.dumps(models,indent=2)+'\n')
        (a.work/f'continuation-chain-{n:03d}-{ig:03d}.json').write_text(json.dumps(continuation_models,indent=2)+'\n')
        return models
    jobs=list(itertools.product(selected,gravity));records=[]
    def checked(job):
        try:return chain(job)
        except Exception as error:
            n,g=job;failure={'plane':n,'gravity':g,'error':str(error)}
            (a.work/f'failure-{n:03d}-{spec["log_g"].index(g):03d}.json').write_text(json.dumps(failure,indent=2)+'\n')
            print('REJECTED',failure,flush=True);return None
    with ThreadPoolExecutor(a.jobs) as pool:
        for result in pool.map(checked,jobs):
            if result is not None:records.extend(result)
    (a.work/'models.json').write_text(json.dumps(records,indent=2)+'\n')
    if len(records)!=len(jobs)*sum(t>=a.minimum_teff for t in spec['teff_K']):raise RuntimeError('incomplete condensate chains')


if __name__=='__main__':main()
