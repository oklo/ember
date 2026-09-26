#!/usr/bin/env python3
"""Finite passive metal redistribution on the saved evolving H/He structures.

This controls material-coverage needs; it is not a coupled stellar track.
Uses trace, fully stripped ions and either classical or suppressed heat flow.
Carbon and nitrogen are passive transport controls here, without reactions.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import time

import numpy as np

from trace_metal_transport import ELEMENTS, coefficients, mass_weights, bernoulli, step
from prepare_nongrey_sources import digest
from write_scientific_result import write_result

YEAR=31557600.


def transport_frame(item):
    index,model=item
    value=np.zeros((2,len(ELEMENTS),len(model)-1,2))
    valid=np.minimum(model[:-1,3],model[1:,3])>=2e6
    for face in np.flatnonzero(valid):
        for k,(_,a,z) in enumerate(ELEMENTS):
            for h,heat in enumerate(["classical","suppressed"]):
                value[h,k,face]=coefficients(model[face],model[face+1],a,z,heat)
    return index,value,valid


def partition(count,outer,inner):
    regions=[]
    if inner>1: regions.append((0,inner))
    regions.extend((i,i+1) for i in range(max(inner,0),outer))
    regions.append((outer,count))
    return regions


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--workers",type=int,default=6)
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError("preserve completed report")
    a.stage.mkdir(exist_ok=False)
    started=time.monotonic();root=Path(__file__).resolve().parents[1]
    inputs={str(p.resolve()):digest(p) for p in [Path(__file__),root/"scripts/trace_metal_transport.py",
            root/"scripts/diffusion_burgers.py",root/"scripts/ion_collision_integrals.py"]}
    checkpoints=[];ages=[];models=[];outer=[];inner=[]
    for version in range(5,23):
        p=Path(f"/tmp/ember-cn-main-sequence-v{version}/checkpoint.json")
        inputs[str(p)]=digest(p);d=json.loads(p.read_text())
        checkpoints.append(p);ages.append(d["age_seconds"])
        models.append(np.array(d["model_record"]["model"])[:,:7])
        mixed=[(i,j) for i,j in d["model_record"]["mixing_regions"] if j-i>1]
        assert mixed[-1][1]==len(models[-1])
        outer.append(mixed[-1][0]);inner.append(mixed[0][1] if len(mixed)>1 and mixed[0][0]==0 else 0)
    times=np.array(ages)-ages[0];mass=models[0][:,0]
    assert np.all(np.diff(times)>0) and all(np.array_equal(m[:,0],mass) for m in models)
    w=mass_weights(mass);total=mass[-1];count=len(mass)
    def at(t):
        j=min(int(np.searchsorted(times,t,side="right"))-1,len(times)-2);j=max(j,0)
        f=(t-times[j])/(times[j+1]-times[j]);lo,hi=models[j:j+2]
        model=lo.copy();model[:,1:4]=np.exp((1-f)*np.log(lo[:,1:4])+f*np.log(hi[:,1:4]))
        model[:,4:]=(1-f)*lo[:,4:]+f*hi[:,4:]
        b=int(np.rint((1-f)*outer[j]+f*outer[j+1]))
        c=int(np.rint((1-f)*inner[j]+f*inner[j+1]))
        # A one-cell convective core is just a singleton.
        if c==1:c=0
        return model,b,c
    frame_times=[]
    for t0,t1 in zip(times[:-1],times[1:]):
        n=max(1,int(math.ceil((t1-t0)/(4e9*YEAR))))
        frame_times.extend(np.linspace(t0,t1,n+1)[:-1])
    frame_times.append(times[-1]);frame_times=np.array(frame_times)
    frames=[at(t)[0] for t in frame_times]
    coefficients_all=np.zeros((len(frames),2,len(ELEMENTS),count-1,2))
    masks=np.zeros((len(frames),count-1),dtype=bool)
    write_result(a.stage/"plan.json",dict(scope=__doc__,started_utc=datetime.now(timezone.utc).isoformat(),
        input_sha256=inputs,frame_count=len(frames),workers=a.workers,
        initial_age_years=ages[0]/YEAR,final_age_years=ages[-1]/YEAR))
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        pending=[pool.submit(transport_frame,(i,m)) for i,m in enumerate(frames)]
        for done,future in enumerate(as_completed(pending),1):
            index,value,valid=future.result();coefficients_all[index]=value;masks[index]=valid
            print(f"transport frames {done}/{len(frames)}",flush=True)
    cache=a.stage/"coefficients.npz"
    np.savez_compressed(cache,times=frame_times,coefficients=coefficients_all,valid=masks)
    # Check the finite trace-control population at four actual envelope bases.
    trace_checks=[]
    for version in [6,10,19,22]:
        j=version-5;m=models[j];face=outer[j]-1
        for name,A,Z in ELEMENTS:
            for heat in ["classical","suppressed"]:
                small=np.array(coefficients(m[face],m[face+1],A,Z,heat,1e-10))
                larger=np.array(coefficients(m[face],m[face+1],A,Z,heat,1e-8))
                error=float(max(abs(larger/small-1)))
                assert error<1e-5
                trace_checks.append(dict(version=version,element=name,heat_flow=heat,relative_error=error))
    def integrate(dt_years):
        values=np.ones((2,len(ELEMENTS),count));t=0.;history=[];error=0.;steps=0
        for frame in range(len(frame_times)-1):
            t0,t1=frame_times[frame:frame+2]
            while t<t1:
                dt=min(dt_years*YEAR,t1-t);mid=t+.5*dt
                f=(mid-t0)/(t1-t0);model,b,c=at(mid)
                regions=partition(count,b,c)
                interfaces=np.array([end-1 for begin,end in regions[:-1]])
                valid=masks[frame]&masks[frame+1]
                if not valid[interfaces].all():raise ValueError("missing hot trace coefficients on active interface")
                v=(1-f)*coefficients_all[frame,:,:,:,0]+f*coefficients_all[frame+1,:,:,:,0]
                D=np.ones_like(v)
                D[:,:,valid]=np.exp((1-f)*np.log(coefficients_all[frame,...,1][:,:,valid])
                                   +f*np.log(coefficients_all[frame+1,...,1][:,:,valid]))
                radius=.5*(model[:-1,1]+model[1:,1]);rho=.5*(model[:-1,2]+model[1:,2])
                area=4*np.pi*radius**2;dr=np.diff(mass)/(area*rho)
                pe=v*dr/D;rate=area*rho*D/dr
                forward=rate*bernoulli(-pe);backward=rate*bernoulli(pe)
                # Unused interfaces inside convective reservoirs have no
                # transport coefficients; mixing is handled by the partition.
                forward[:,:,~valid]=0;backward[:,:,~valid]=0
                for h in range(2):
                    for k in range(len(ELEMENTS)):
                        values[h,k],e=step(values[h,k],w,forward[h,k],backward[h,k],regions,dt)
                        error=max(error,e)
                t+=dt;steps+=1
            history.append(dict(age_years=(ages[0]+t)/YEAR,surface_ratio=values[:,:,-1].tolist(),
                                maximum_ratio=np.max(values,axis=2).tolist(),outer_base_index=b))
        assert max(abs(np.sum(values*w,axis=2)/total-1).flat)<1e-8
        return values,history,error,steps
    coarse,chistory,cerror,csteps=integrate(2e8)
    fine,fhistory,ferror,fsteps=integrate(1e8)
    endpoint_path=a.stage/"endpoint.npz"
    np.savez_compressed(endpoint_path,mass=mass,weights=w,coarse_ratios=coarse,fine_ratios=fine)
    surface_difference=abs(coarse[:,:,-1]/fine[:,:,-1]-1)
    summary=[]
    for h,heat in enumerate(["classical","suppressed"]):
        summary.append(dict(heat_flow=heat,elements=[dict(element=name,
            surface_over_initial=float(fine[h,k,-1]),maximum_over_initial=float(max(fine[h,k])),
            minimum_over_initial=float(min(fine[h,k])),
            surface_time_refinement=float(surface_difference[h,k]),
            final_global_inventory_ratio=float(w@fine[h,k]/total)) for k,(name,A,Z) in enumerate(ELEMENTS)]))
    assert all(digest(p)==value for p,value in inputs.items())
    result=dict(utc=datetime.now(timezone.utc).isoformat(),outcome="passive_redistribution_completed",selected=False,
        scope=__doc__,initial_age_years=ages[0]/YEAR,final_age_years=ages[-1]/YEAR,
        elapsed_evolution_years=(ages[-1]-ages[0])/YEAR,wall_seconds=time.monotonic()-started,
        frame_count=len(frames),coarse_steps=csteps,fine_steps=fsteps,
        maximum_step_inventory_error=max(cerror,ferror),
        maximum_surface_time_refinement=float(max(surface_difference.flat)),trace_limit_checks=trace_checks,
        summary=summary,history=fhistory,input_sha256=inputs,
        artifacts_sha256={str(p):digest(p) for p in [cache,endpoint_path,a.stage/"plan.json"]},
        limitations=["Prescribed H/He structure and gravity; no composition, energy or opacity feedback on the star.",
          "Uniform passive tracers start at the 3.551 Tyr saved structure; prior short radiative interval is omitted.",
          "Classical versus suppressed heat flow are separate approximations, not uncertainty bounds.",
          "Fully stripped trace metals; H/He background renormalized after removing its fixed 2% metal mass.",
          "Core electron degeneracy and nonideal ion forces are not fully described by the classical heat-flow control.",
          "No C/N nuclear conversion; C and N are passive transport controls.",
          "Intermediate structures and convective boundaries are interpolated between retained checkpoints.",
          "Timestep comparison does not assess spatial mesh or saved-background interpolation error."])
    write_result(a.output,result)
    print(json.dumps({k:result[k] for k in ["outcome","wall_seconds","maximum_surface_time_refinement","summary"]},indent=2))


if __name__=="__main__":main()
