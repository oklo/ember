#!/usr/bin/env python3
"""Compare selected conductive material data with retained microscopic faces.

Query the original conductive prescription at nodes and microscopic face states.
This is a material-coefficient comparison, not a new stellar structure solution.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess

from fetch_tops_composition import digest
from reduce_tops_group_factors import verify


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args()
    if args.work.exists() or args.report.exists():raise FileExistsError('preserve previous extract')
    args.work.mkdir()
    profile=Path('docs/reports/2026-09-11/evolution_latest_profile.csv')
    probe=Path('/private/tmp/ember-conduction-mixture-probe-v1')
    table=Path('data/conduction/condtab21wd_metals.dat')
    old=Path('docs/results/tops_density_fine_transport_v1.json')
    old_data=json.loads(old.read_text())
    for p in (probe,table):
        if old_data['input_sha256'][str(p.resolve())]!=digest(p):raise ValueError('selected conduction input changed')
    microscopic=Path('docs/results/screened_microscopic_transport_v3.json')
    micro=json.loads(microscopic.read_text())
    if micro['input_sha256'][str(profile.resolve())]!=digest(profile):raise ValueError('profile changed')
    inputs={str(p.resolve()):digest(p) for p in (Path(__file__),profile,probe,table,old,microscopic,
        Path('scripts/conduction_mixture_probe.cpp'),Path('src/conduction_table.cpp'),
        Path('include/ember/conduction_table.hpp'),Path('src/screened_microscopic_transport.cpp'))}
    if digest(Path('scripts/conduction_mixture_probe.cpp'))!=old_data['probe_source_sha256']:
        raise ValueError('conduction probe source identity changed')
    with profile.open() as f:rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
    def state(r):return [r['X'],1-r['X']-r['Y3']-r['Y4'],r['Y3'],r['temperature_K'],r['density_g_cm3']]
    queries=[state(rows[i]) for i in range(380,422)]
    for i in range(380,421):
        a,b=state(rows[i]),state(rows[i+1])
        q=[.5*(u+v) for u,v in zip(a,b)]
        q[3]=math.sqrt(a[3]*b[3])
        queries.append(q)
    text=''.join(' '.join(format(v,'.17g') for v in q)+'\n' for q in queries)
    (args.work/'queries.txt').write_text(text)
    reply=subprocess.run([str(probe),str(table)],input=text,text=True,capture_output=True,timeout=30)
    (args.work/'stdout.txt').write_text(reply.stdout);(args.work/'stderr.txt').write_text(reply.stderr)
    reply.check_returncode()
    kappas=[float(line) for line in reply.stdout.splitlines()]
    if len(kappas)!=len(queries) or any(not math.isfinite(k) or k<=0 for k in kappas):raise ValueError('invalid conduction reply')
    arad=7.565733250033928e-15;clight=2.99792458e10
    coefficients=[4*arad*clight*q[3]**3/(3*q[4]*k) for q,k in zip(queries,kappas)]
    saved={r['face']:r for r in micro['faces'] if r['kind']=='saved'}
    faces=[]
    for i in range(380,421):
        j=42+i-380;a,b=rows[i],rows[i+1];r=saved[i]['response'];q=queries[j]
        radius=.5*(a['radius_cm']+b['radius_cm']);area=4*math.pi*radius**2
        geom=area**2*q[4]/(b['mass_g']-a['mass_g'])
        nominal=-coefficients[j]*geom*(b['temperature_K']-a['temperature_K'])
        lum=.5*(a['luminosity_erg_s']+b['luminosity_erg_s'])
        faces.append(dict(face=i,query_X_Z_Y3_T_rho=q,selected_conductive_opacity=kappas[j],
            selected_conductivity_cgs=coefficients[j],microscopic_conductivity_cgs=r['conductivity'],
            microscopic_relative_to_selected=r['conductivity']/coefficients[j]-1,
            selected_local_conductive_luminosity=nominal,microscopic_conductive_luminosity=r['conductive'],
            conductive_luminosity_difference_over_stellar=(r['conductive']-nominal)/lum,
            microscopic_carried_luminosity=r['carried'],microscopic_carried_over_stellar=r['carried']/lum,
            microscopic_species_rate=r['rate'],stellar_luminosity=lum))
    nodes=[dict(zone=i,query=queries[j],kappa=kappas[j],conductivity=coefficients[j]) for j,i in enumerate(range(380,422))]
    verify(inputs)
    report=dict(scope=__doc__,outcome='completed_material_comparison',accepted_for_stellar_evolution=False,
        basis='baryonic isotope fractions and density; GS98 metal inventory',nodes=nodes,faces=faces,
        face_state_convention='arithmetic radius, density and abundances; geometric temperature, matching the microscopic evaluator',
        limitations=['The selected conductivity is evaluated at the microscopic face state; this is not the nodal opacity average in the selected structure equation.',
            'The conductive luminosities use the saved temperature gradient; convection and structure have not readjusted.',
            'The microscopic prescription assumes fully ionized collision partners and is not accepted in the cool envelope.'],
        input_sha256=inputs,output_sha256={str(p.resolve()):digest(p) for p in args.work.iterdir() if p.is_file()})
    args.report.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    base=next(r for r in faces if r['face']==396)
    print(json.dumps(dict(nodes=len(nodes),faces=len(faces),base=base)),flush=True)


if __name__=='__main__':main()
