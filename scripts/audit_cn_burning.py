#!/usr/bin/env python3
"""Check coupled pp/CN burning and convection on three saved stellar structures."""
import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--native',type=Path,required=True)
    ap.add_argument('--scratch',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();a.scratch.mkdir()
    if a.output.exists():raise ValueError('preserve existing result')
    identities={}
    def pin(p):
        p=Path(p);identities[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();return p
    for p in [a.native,Path(__file__),'scripts/cn_burning_probe.cpp','include/ember/cn_burning.hpp','src/cn_burning.cpp']:
        pin(p)
    C0=.02*.171836/12;N0=.02*.050335/14
    cases=[('fully_convective_2500gyr','out/evolution-metal-512-2500gyr-gas-checkpoint.json',0.,1e9),
           ('fully_convective_3400gyr','out/evolution-cold-remnant-forward-512-3400gyr-v1.json',.5,1e8),
           ('late_800myr','/tmp/ember-adaptive-diffusion-v11/checkpoint.json',1.,1e6)]
    queries=[];replies=[];summaries=[]
    def physical(c):
        out=np.array(c,copy=True)
        out[:,3:]*=[12,13,14]
        extra=out[:,3:].sum(axis=1)-(12*C0+14*N0)
        out[:,2]-=extra
        return out
    for name,path,fraction,years in cases:
        data=json.loads(pin(path).read_text())
        if 'profile' in data:
            model=np.array(data['profile'])[:,:7];regions=[[0,len(model)]]
        else:
            model=np.array(data['model_record']['model']);regions=data['model_record']['mixing_regions']
        m=model[:,0];weights=np.empty(len(m));weights[0]=(m[0]+m[1])/2;weights[-1]=(m[-1]-m[-2])/2;weights[1:-1]=(m[2:]-m[:-2])/2
        old=np.column_stack([model[:,5:7],.98-model[:,5:7].sum(axis=1),np.tile([C0*(1-fraction),0,N0+C0*fraction],(len(m),1))])
        def step(composition,dt,part):
            rows=np.column_stack([model[:,:5],composition[:,:2],composition[:,3:]])
            values=[len(m),len(regions),dt,1e-13,*rows.flat,*np.array(regions).flat]
            query=' '.join(format(v,'.17g') for v in values)+'\n';queries.append(query)
            run=subprocess.run([str(a.native)],input=query,text=True,capture_output=True,check=True,timeout=120)
            (a.scratch/(name+'-'+part+'.stderr')).write_text(run.stderr)
            result=json.loads(run.stdout);replies.append(result)
            if 'error' in result:raise RuntimeError(result['error'])
            new=np.array(result['composition']);source=np.array(result['sources'])
            assert np.isfinite(new).all() and np.min(new)>=0
            assert abs(new[:,3:].sum(axis=1)/(C0+N0)-1).max()<2e-12
            assert abs(result['nuclear_mass_balance'])<2e-6
            residual=physical(new)-physical(composition)-dt*source
            region_balance=[]
            for begin,end in regions:
                averaged=np.dot(weights[begin:end],residual[begin:end])/m[-1]
                region_balance.append(float(abs(averaged).max()))
                if end>begin+1:assert np.all(new[begin:end]==new[begin])
            maximum=max(region_balance)
            assert maximum<2e-12,maximum
            return new,{k:v for k,v in result.items() if k not in ['composition','sources']}|{'maximum_region_balance':maximum}
        seconds=years*365.25*86400
        coarse,c1=step(old,seconds,'full')
        first,c2=step(old,seconds/2,'half1')
        fine,c3=step(first,seconds/2,'half2')
        (a.scratch/(name+'-fine.json')).write_text(json.dumps({'composition':fine.tolist(),'model':model.tolist(),'regions':regions})+'\n')
        summaries.append({'case':name,'nodes':len(m),'mixed_regions':sum(end>begin+1 for begin,end in regions),
                          'initial_carbon_converted_fraction':fraction,'interval_years':years,'steps':[c1,c2,c3],
                          'full_vs_two_halves_maximum_physical_abundance_difference':float(abs(physical(coarse)-physical(fine)).max()),
                          'fine_central_lookup_H_He3_He4_and_CN_molalities':fine[0].tolist(),
                          'hydrogen_consumed_g':float(weights@(old[:,0]-fine[:,0]))})
    (a.scratch/'queries.txt').write_text(''.join(queries));(a.scratch/'responses.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in replies))
    for p in a.scratch.iterdir():
        if p.is_file():pin(p)
    report={'created_utc':datetime.now(timezone.utc).isoformat(),'outcome':'passed_fixed_structure_burning_controls',
            'cases':summaries,'input_sha256':identities,'limitations':[
                'Temperature, density, mass mesh and mixing regions are held at each saved model.',
                'Chosen initial carbon fractions are controls, not reconstructed physical abundances or restart selections.',
                'No microscopic/secular transport, atmosphere response or stellar thermal evolution is calculated.',
                'Full-versus-half-step differences are reported; these intervals are not asserted to satisfy production time-accuracy targets.']}
    a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'outcome':report['outcome'],'cases':summaries},indent=2))


if __name__=='__main__':main()
