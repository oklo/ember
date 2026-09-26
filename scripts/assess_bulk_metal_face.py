#!/usr/bin/env python3
"""Combine variable-metal material forces and the moving-group collision law.

Uses the actual current convective-envelope base, retaining its composition
and screening length. A local reservoir time is a scale, not an evolved
surface abundance. The accepted stellar model remains unchanged.
"""
import argparse
from datetime import datetime,timezone
import gzip,hashlib,json,math
from pathlib import Path
import subprocess
import numpy as np
from trace_metal_transport import mass_weights
from write_scientific_result import write_result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    stage=Path('/tmp/ember-bulk-metal-face-v1');stage.mkdir()
    checkpoint=Path('/tmp/ember-cn-main-sequence-v22/checkpoint.json')
    reference=Path('docs/results/thoul_diffusion_comparison_v1.json')
    eos=Path('/tmp/ember-variable-metal-eos-v1/probe');family=Path('/tmp/ember-variable-metal-eos-v1/family/variable_metal.dat')
    collision=Path('/tmp/ember-bulk-metal-collision-v1/probe');table=Path('/tmp/ember-collision-eta32-lowb-192-v1.dat')
    record=next(r for r in json.loads(reference.read_text())['records'] if r['version']==22 and r['face']==386)
    model=np.array(json.loads(checkpoint.read_text())['model_record']['model']);left,right=model[386:388]
    x=.5*(left[5]+right[5]);y=.5*(left[6]+right[6]);T=math.sqrt(left[3]*right[3]);rho=.5*(left[2]+right[2]);Z=.02
    radius=.5*(left[1]+right[1]);area=4*math.pi*radius**2;geometry=area**2*rho/(right[0]-left[0])
    queries=[[p[5],p[6],Z,p[3],p[2],1,1,1]for p in [left,right]]+[[x,y,Z,T,rho,1,1,1]]
    text=''.join(' '.join(format(float(v),'.17g')for v in q)+'\n'for q in queries)
    run=subprocess.run([str(eos),str(family)],input=text,text=True,capture_output=True,check=True)
    (stage/'material.jsonl').write_text(run.stdout);material=[json.loads(line) for line in run.stdout.splitlines()]
    if any('error'in r for r in material):raise ValueError(material)
    query=[T,rho,x,y,Z,record['native']['screening_length']]
    run=subprocess.run([str(collision),str(table)],input=' '.join(format(float(v),'.17g')for v in query)+'\n',text=True,capture_output=True,check=True)
    (stage/'collision.json').write_text(run.stdout);kinetic=json.loads(run.stdout)
    if 'error'in kinetic:raise ValueError(kinetic)
    hmaterial=np.array(material[2]['enthalpy']);hkinetic=np.array(kinetic['bulk']['enthalpy']);h=hmaterial+hkinetic
    chemical=np.array(material[1]['gradient'])-np.array(material[0]['gradient'])
    inverse_temperature=(right[3]-left[3])/(left[3]*right[3])
    force=chemical+h*inverse_temperature
    mobility=np.array(kinetic['bulk']['mobility'])[:3,:3];rate=-geometry*mobility@force
    fixedh=hmaterial[:2]+np.array(kinetic['fixed']['enthalpy'])
    fixedforce=chemical[:2]+fixedh*inverse_temperature
    fixedrate=-geometry*np.array(kinetic['fixed']['mobility'])[:2,:2]@fixedforce
    mass=mass_weights(model[:,0]);envelope_mass=float(sum(mass[387:]));yr=365.25*86400
    seconds=envelope_mass*Z/-rate[2] if rate[2]<0 else None
    # Opposing a concentration gradient must diffuse a positive metal peak.
    hessian=np.array(material[2]['hessian']);diffusivity=float((mobility@hessian[:,2])[2]/rho)
    if not np.isfinite(rate).all() or diffusivity<=0:raise ValueError('invalid group face response')
    result=dict(utc=datetime.now(timezone.utc).isoformat(),outcome='local_bulk_metal_face_evaluated',selected=False,
                age_years=3703913421275.8203,face=386,temperature_K=T,density_g_cm3=rho,
                composition=dict(X=x,Y3=y,Z=Z),species_order=['H1','He3','GS98_metal_group'],
                mass_rates_g_per_s=rate.tolist(),mass_velocities_cm_per_s=(rate/(area*rho*np.array([x,y,Z]))).tolist(),
                local_metal_reservoir_years=seconds/yr if seconds else None,
                local_metal_gradient_diffusivity_cm2_per_s=diffusivity,
                metal_chemical_force=force.tolist(),material_exchange_enthalpy_erg_per_g=hmaterial.tolist(),
                kinetic_transport_enthalpy_erg_per_g=hkinetic.tolist(),
                fixed_metal_rate_with_same_EOS_g_per_s=fixedrate.tolist(),
                mobile_over_fixed_H_rate=float(rate[0]/fixedrate[0]),
                same_EOS_fixed_over_selected_H_rate=float(fixedrate[0]/record['native']['native_rate'][0]),
                material_carried_luminosity_erg_per_s=float(hmaterial@rate),
                kinetic_carried_luminosity_erg_per_s=float(hkinetic@rate),
                input_sha256={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest()for p in
                              [Path(__file__),checkpoint,reference,eos,family,collision,table]},
                limitations=['One local interface, not a coupled evolutionary step or integrated depletion history.',
                             'Fully stripped GS98 metals share one velocity; distinct metal heat variables are retained.',
                             'No atmosphere or structural feedback, changed nuclear ledger, or global material coverage is supplied here.'])
    write_result(a.output,result)
    print(json.dumps({k:result[k]for k in ['outcome','mass_rates_g_per_s','local_metal_reservoir_years','local_metal_gradient_diffusivity_cm2_per_s','mobile_over_fixed_H_rate','same_EOS_fixed_over_selected_H_rate']},indent=2))


if __name__=='__main__':main()
