#!/usr/bin/env python3
"""Exercise coupled physical hot-interior evolution with a held outer envelope.

This is a solver control on a truncated saved interior, not a whole-star
trajectory. The pressure and temperature of its sealed outer edge are fixed.
"""
import argparse,csv,gzip,hashlib,json,subprocess,time
from datetime import datetime,timezone
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('probe',type=Path);p.add_argument('output',type=Path);p.add_argument('--scratch',type=Path,required=True)
    a=p.parse_args();assert not a.output.exists() and not a.scratch.exists();a.scratch.mkdir()
    root=Path(__file__).resolve().parents[1];sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    family=Path('/private/tmp/ember-smooth-eos-refined-family-v1/freeeos300_gs98_z020.dat')
    collision=Path('/tmp/ember-collision-transport-table-v1.dat')
    opacity=Path('/private/tmp/ember-refractive-hot-family-refined-v1')
    profile=root/'docs/reports/2026-09-11/evolution_latest_profile.csv'
    with profile.open() as f:rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)][:397]
    fields=['mass_g','radius_cm','density_g_cm3','temperature_K','luminosity_erg_s','X','Y3']
    model=' '.join(format(r[k],'.17g') for r in rows for k in fields)
    cases=[dict(ions=ions,absent_He3=absent,years=1e6) for ions in [0,1] for absent in [0,1]]
    queries=[f"{c['ions']} {c['absent_He3']} {c['years']*365.25*86400:.17g} {len(rows)} "+model for c in cases]
    query=a.scratch/'queries.txt';query.write_text('\n'.join(queries)+'\n')
    sources=[Path(__file__),a.probe,family,collision,profile,root/'scripts/physical_core_step_probe.cpp',
       root/'include/ember/microscopic_transport.hpp',root/'src/microscopic_transport.cpp',
       root/'src/screened_microscopic_transport.cpp',root/'include/ember/screened_microscopic_transport.hpp',
       root/'src/species_transport.cpp',root/'src/evolution.cpp',root/'src/mixing.cpp',root/'src/structure.cpp',
       root/'src/henyey.cpp',root/'src/relaxation.cpp',root/'src/eos_smooth_mixture.cpp',root/'src/collision_transport.cpp',
       root/'src/opacity_mixture.cpp',root/'src/opacity_table.cpp',root/'src/nuclear_pp.cpp',
       opacity/'tops_gs98_mixture_high.dat']+[opacity/f'tops_gs98_mixture_z{z:03d}_high.dat' for z in [10,20,30]]
    identities={str(p):sha(p) for p in sources}
    assert identities[str(family)]=='7d6bfa0b26264834307dad651852bc5ff693f91f3861b4934e582dee047fd7f5'
    assert identities[str(collision)]=='48f7328da79a6163a997ab524c623ecb334d964f47bb7fbdefc084e63b46e90f'
    stdout=a.scratch/'responses.jsonl';stderr=a.scratch/'stderr.txt'
    start=time.monotonic();error=None;code=None
    with query.open() as inp,stdout.open('x') as out,stderr.open('x') as err:
        try:code=subprocess.run([str(a.probe),str(family),str(collision),str(opacity)],stdin=inp,stdout=out,stderr=err,timeout=160,check=False).returncode
        except subprocess.TimeoutExpired:error='native probe exceeded 160 seconds; partial replies retained'
    elapsed=time.monotonic()-start
    replies=[json.loads(s) for s in stdout.read_text().splitlines()]
    failures=[];summaries=[]
    if error:failures.append(error)
    if code not in [0,None]:failures.append(f'process exit {code}')
    if len(replies)!=len(cases):failures.append(f'{len(replies)} replies for {len(cases)} queries')
    for i,(c,r) in enumerate(zip(cases,replies)):
        summary={**c,**{k:v for k,v in r.items() if k!='model'}};summaries.append(summary)
        if 'error' in r or not r.get('converged'):failures.append(dict(case=i,message=r.get('error',r.get('message'))));continue
        if not r['input_preserved'] or r['age']!=c['years']*365.25*86400:failures.append(dict(case=i,message='age or old model changed incorrectly'))
        if abs(r['luminosity_balance'])>2e-8 or abs(r['nuclear_mass_balance'])>2e-6 or r['abundance_residual']>1e-12:
            failures.append(dict(case=i,message='conservation or coupling criterion failed'))
        if not (r['max_dlnT']>0 and r['max_dlnrho']>0 and r['max_dX']>0):
            failures.append(dict(case=i,message='thermal, mechanical or composition variables did not advance'))
        if c['absent_He3'] and any(row[-1]<=0 for row in r['model']):
            failures.append(dict(case=i,message='nuclear He3 production failed'))
    for name,h in identities.items():assert sha(name)==h,name
    for f in [query,stdout]:(f.with_suffix(f.suffix+'.gz')).write_bytes(gzip.compress(f.read_bytes(),mtime=0))
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed_held_envelope_physical_core_controls' if not failures else 'failed',
       accepted_for_stellar_evolution=False,points=len(rows),enclosed_mass_fraction=rows[-1]['mass_g']/(.1*1.98847e33),
       queries=len(cases),completed_cases=len(replies),elapsed_seconds=elapsed,failures=failures,summaries=summaries,
       input_sha256=identities,artifacts_sha256={str(p):sha(p) for p in a.scratch.iterdir()},
       limitations=['Hot interior truncated at point 397, with fixed outer pressure and temperature and sealed species flux.',
         'Conditional fully stripped collision model; corrected hot radiative opacity and native EOS.',
         'No cool envelope, photosphere, moving truncation boundary, or continuation of the selected stellar trajectory.',
         'This checks coupled equations; it does not establish timestep accuracy or omitted physical effects.'])
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ['outcome','elapsed_seconds','failures','summaries']},indent=2))
    if failures:raise SystemExit(1)
if __name__=='__main__':main()
