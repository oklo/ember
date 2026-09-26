#!/usr/bin/env python3
"""Check absolute composition counts and direct variable-Z FreeEOS sources."""
import argparse,json,math,subprocess
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from metal_eos_composition import mixture
from prepare_nongrey_sources import digest
from write_scientific_result import write_result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('output',type=Path);a=ap.parse_args()
    if a.output.exists():raise FileExistsError('preserve results')
    probe=Path('/tmp/ember-freeeos-precision-builder-check-v1/probe')
    assert digest(probe)=='44c284b28594ca6d19c0d69a1402686f8dbb758746d220bced89927f0eb15a4f'
    cp=Path('/tmp/ember-cn-main-sequence-v22/checkpoint.json')
    original=Path('/tmp/ember-variable-metal-source-controls-v1/original_mixtures.json')
    paths=[Path(__file__),Path('scripts/metal_eos_composition.py'),Path('scripts/generate_metal_eos.py'),
           Path('scripts/stellar_composition.py'),probe,cp,original,
           Path('data/atmosphere/sources/synple-elements.json'),Path('data/opacity/sources/tops_gs98_x070_z020.request.json')]
    inputs={str(p.resolve()):digest(p) for p in paths}
    for r in json.loads(original.read_text()):assert mixture(r['X'],r['Y3'])==r['mixture']
    for x,y,z in [(1,0,0),(0,0,0),(0,.06,.003),(.999,0,.001),(.98,.015,.003)]:
        m=mixture(x,y,metallicity=z);assert abs(sum(m['composition'])-1)<1e-14
        if x==0:assert m['eps'][0]==0
        if x+z==1 and y==0:assert m['eps'][1]==0
        if z==0:assert all(e==0 for e in m['eps'][2:])
    rejected=0
    for x,y,z in [(.99,.02,.001),(.9,0,-.001),(.9,0,float('nan'))]:
        try:mixture(x,y,metallicity=z)
        except ValueError:rejected+=1
        else:raise ValueError('invalid source mixture accepted')
    model=np.array(json.loads(cp.read_text())['model_record']['model'])
    cases=[dict(zone=i,X=p[5],Y3=p[6],T=p[3],rho=p[2]) for i,p in enumerate(model) if i in [0,150,300,387,430,480,511]]
    cases += [dict(zone=None,X=x,Y3=y,T=t,rho=r) for x,y,t,r in
              [(1.,0.,6000.,1e-5),(0.,0.,1e7,1e4),(.999,0.,1e5,.01)]]
    responses=[]
    def source(m,T,rho):
        scale=m['source_mass_scale']
        request=' '.join(format(v,'.17g') for v in m['eps'])+'\n3 223 -2\n'+f'{math.log(rho*scale):.17g} {math.log(T):.17g}\n'
        p=subprocess.run([str(probe)],input=request,text=True,capture_output=True,timeout=60,check=True)
        v=np.array([float(s) for s in p.stdout.split()])
        assert len(v)==22 and np.isfinite(v).all() and v[0]==0,(m['hydrogen'],T,rho,p.stdout,p.stderr)
        assert abs(v[2]/scale/rho-1)<1e-8 and abs(v[3]/T-1)<1e-8
        return dict(request=request,raw_source=v.tolist(),pressure=v[4],energy=v[5]*scale,
                    cv=v[12]*scale,cp=v[13]*scale,adiabatic_gradient=v[14])
    fixed_errors=[]
    for case in cases:
        zs=[.001,.003,.01,.02] if case['zone'] is not None else ([0.] if case['X'] in [0.,1.] else [.001])
        states=[]
        for z in zs:
            m=mixture(case['X'],case['Y3'],metallicity=z)
            s=source(m,case['T'],case['rho']);states.append(dict(Z=z,mixture=m,**s))
        if case['zone'] is not None:
            old=source(mixture(case['X'],case['Y3']),case['T'],case['rho'])
            error=max(abs(states[-1][k]/old[k]-1) for k in ['pressure','energy','cv','cp'])
            assert error<1e-7
            fixed_errors.append(error)
            for s in states:
                s['relative_to_Z020']={k:float(s[k]/states[-1][k]-1) for k in ['pressure','energy','cv','cp']}
        responses.append(dict(**case,states=states))
    assert all(digest(p)==v for p,v in inputs.items())
    result=dict(utc=datetime.now(timezone.utc).isoformat(),outcome='variable_metal_source_controls_passed',
        selected=False,default_mixtures_unchanged=6,invalid_compositions_rejected=rejected,
        maximum_explicit_Z020_source_difference=max(fixed_errors),queries=sum(len(r['states']) for r in responses)+len(fixed_errors),
        records=responses,input_sha256=inputs,
        limitations=['Direct FreeEOS source values, not a variable-Z interpolation family or selected stellar EOS.',
                     'Helium uses source He4 electronic physics; analytic isotope entropy is applied during table import.',
                     'Uniformly scaled GS98 pattern; selective elemental changes need their own material treatment.'])
    write_result(a.output,result)
    print(json.dumps({k:result[k] for k in ['outcome','queries','maximum_explicit_Z020_source_difference']},indent=2))


if __name__=='__main__':main()
