#!/usr/bin/env python3
"""Measure H/He composition curvature before using EOS-derived diffusion forces.

All controls hold baryonic rho, T, Z=.02 and Y3=0 fixed while H1 replaces He4.
The temperatures and densities come from three saved stellar layers; X is
the midpoint of its source-table cell, so every difference remains inside one
cell. Additional controls span larger hydrogen fractions at central T/rho.
FreeEOS options match the selected family. Two independently built electron
quadrature tolerances and two finite-difference spacings test numerical error.
This diagnostic changes neither the EOS nor any stellar abundance.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import subprocess

import numpy as np

from metal_eos_composition import mixture

RGAS = 6.02214076e23*1.380649e-16


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def differences(f, h):
    # Symmetric differences reduce the composition-independent entropy offset.
    second = ((f[2]-f[0])+(f[-2]-f[0]))/(4*h*h)
    fine = ((f[1]-f[0])+(f[-1]-f[0]))/(h*h)
    finer = ((f[.5]-f[0])+(f[-.5]-f[0]))/(.25*h*h)
    return {'coarse':second, 'fine':fine, 'finer':finer,
            'richardson':(4*fine-second)/3,
            'richardson_fine':(4*finer-fine)/3,
            'slope':(f[1]-f[-1])/(2*h)}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('ember_probe',type=Path)
    ap.add_argument('work',type=Path)
    ap.add_argument('output',type=Path)
    a=ap.parse_args();root=Path(__file__).resolve().parents[1]
    family=root/'data/eos/numerical_electron_base_v1/freeeos300_gs98_z020.dat'
    profile=root/'docs/reports/2026-09-11/evolution_latest_profile.csv'
    a.work.mkdir(exist_ok=False)
    axis=np.array([float(x) for x in family.read_text().splitlines()[1].split()[2:]])
    rows=np.genfromtxt(profile,delimiter=',',names=True)
    cases=[]
    for label,i,target in [('center',0,None),('interior',250,None),
                           ('outer_nonconvective',396,None),
                           ('hydrogen_01625',0,.1625),('hydrogen_065',0,.65)]:
        row=rows[i];x=float(row['X']) if target is None else target
        j=int(np.searchsorted(axis,x,side='right')-1)
        lo,hi=map(float,axis[j:j+2]);mid=(lo+hi)/2
        cases.append({'label':label,'zone':i,'saved_X':float(row['X']),
                      'T':float(row['temperature_K']),'rho':float(row['density_g_cm3']),
                      'X':mid,'Y3':0.,'source_interval':[lo,hi],'h':(hi-lo)/8})
    builds=[]
    for name in ['ember-freeeos-precision-builder-check-v1',
                 'ember-freeeos-quadrature13-builder-check-v1']:
        directory=Path('/tmp')/name
        receipt=json.loads((directory/'build_receipt.json').read_text())
        probe=directory/'probe'
        assert digest(probe)==receipt['probe_sha256']
        assert digest(receipt['library'])==receipt['library_sha256']
        builds.append({'probe':str(probe),'receipt':receipt,
                       'receipt_sha256':digest(directory/'build_receipt.json')})
    points=[-2.,-1.,-.5,0.,.5,1.,2.]
    queries=[(ci,z,c['X']+z*c['h']) for ci,c in enumerate(cases) for z in points]
    query_text=''.join(f'{x:.17g} 0 {cases[i]["T"]:.17g} {cases[i]["rho"]:.17g}\n'
                       for i,z,x in queries)
    (a.work/'ember_queries.txt').write_text(query_text)
    p=subprocess.run([str(a.ember_probe.resolve()),str(family)],input=query_text,
                     text=True,capture_output=True,timeout=90)
    (a.work/'ember_stdout.txt').write_text(p.stdout)
    (a.work/'ember_stderr.txt').write_text(p.stderr)
    if p.returncode: raise RuntimeError(p.stderr)
    values=[list(map(float,line.split())) for line in p.stdout.splitlines()]
    assert len(values)==len(queries) and all(len(v)==10 and all(map(math.isfinite,v)) for v in values)
    runtime={}
    for (ci,z,x),v in zip(queries,values,strict=True):
        assert abs(v[0]-x)<1e-15
        runtime[ci,z]=v

    def direct(task):
        bi,qi=task;ci,z,x=queries[qi];c=cases[ci];m=mixture(x,0.)
        mass=m['source_mass_scale']
        inp=' '.join(format(v,'.17g') for v in m['eps'])+'\n3 223 -2\n'
        inp+=f'{math.log(c["rho"]*mass):.17g} {math.log(c["T"]):.17g}\n'
        p=subprocess.run([builds[bi]['probe']],input=inp,text=True,capture_output=True,timeout=20)
        v=list(map(float,p.stdout.split()))
        record={'build':bi,'case':ci,'offset':z,'X':x,'input':inp,
                'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr,
                'source_mass_scale':mass}
        if p.returncode or len(v)!=22 or v[0]!=0 or not all(map(math.isfinite,v)):
            raise RuntimeError(record)
        assert abs(v[2]/(c['rho']*mass)-1)<1e-9
        assert abs(v[3]/c['T']-1)<1e-10
        record.update(phi=(v[5]/c['T']-v[6])*mass,P=v[4],E=v[5]*mass)
        return record
    tasks=[(bi,qi) for bi in range(len(builds)) for qi in range(len(queries))]
    with ThreadPoolExecutor(max_workers=2) as pool: source=list(pool.map(direct,tasks))
    (a.work/'direct_source.json').write_text(json.dumps(source,indent=2)+'\n')
    checked=[]
    for ci,c in enumerate(cases):
        table=differences({z:runtime[ci,z][4] for z in points},c['h'])
        ds=[differences({v['offset']:v['phi'] for v in source if v['build']==bi and v['case']==ci},c['h'])
            for bi in range(len(builds))]
        ref=ds[-1]['richardson_fine']
        step_error=abs(ds[-1]['richardson']/ref-1)
        precision_error=abs(ds[0]['richardson_fine']/ref-1)
        ideal=RGAS*(1/c['X']+1/(4*(.98-c['X'])))
        # These criteria accept the measurement, not the table for diffusion.
        assert ref>0 and step_error<1e-4 and precision_error<1e-4
        assert abs(table['richardson_fine']/ref)<1e-5
        checked.append({**c,'table':table,'source':ds,'source_curvature_over_Rgas':ref/RGAS,
                        'ideal_ion_curvature_over_Rgas':ideal/RGAS,
                        'table_over_source_curvature':table['richardson_fine']/ref,
                        'source_step_relative_change':step_error,
                        'source_quadrature_relative_change':precision_error})
    result={'scope':__doc__,'created_utc':datetime.now(timezone.utc).isoformat(),
            'outcome':'completed_diagnostic_selected_table_lacks_H_composition_curvature',
            'accepted_for_diffusion_forces':False,'new_stellar_evolution':False,
            'freeeos_options':[3,223,-2],'source_builds':builds,'direct_queries':len(source),
            'runtime_queries':len(queries),'controls':checked,
            'interpretation':'At Y3=0, linear interpolation of F/T in X gives zero within-cell second derivative. Direct FreeEOS has positive curvature. Existing P/E derivative and heat-capacity checks do not validate chemical-potential derivatives for diffusion. This diagnostic does not change or invalidate the recorded integration of the implemented stellar equations.',
            'inputs_sha256':{str(p):digest(p) for p in [Path(__file__),root/'scripts/diffusion_eos_potential_probe.cpp',root/'src/eos_mixture.cpp',root/'scripts/metal_eos_composition.py',a.ember_probe,family,profile]},
            'raw_source_file':str(a.work/'direct_source.json'),
            'raw_source_sha256':digest(a.work/'direct_source.json')}
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'outcome':result['outcome'],'cases':len(checked),
                      'max_step_error':max(c['source_step_relative_change'] for c in checked),
                      'max_quadrature_error':max(c['source_quadrature_relative_change'] for c in checked),
                      'max_table_source_curvature':max(abs(c['table_over_source_curvature']) for c in checked)}))


if __name__=='__main__':main()
