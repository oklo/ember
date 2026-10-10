"""Inject existing cluster endpoints into Gaia background and apply CP cuts.

Conditional recovery test: the Gaia photometry/RUWE eligibility of faint mock
objects is not calibrated. Manual cluster isolation and SCCP are omitted.
"""
from pathlib import Path
from contextlib import closing
import sys,json,sqlite3,time
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
W=Path(__file__).resolve().parent;V3=W.parent/'cluster-stellar-fit-oct9-v1'
sys.path.insert(0,str(V3))
from common import V1,V2,C,old,CAL
from build_templates import bound_pairs
from selection import m10,completeness
from cp_selection import coordinates,window,neighbors

def inverse_mass(q,lum,rng):
    lm=np.log10(q.mass_Msun.to_numpy());ll=q.logL_Lsun.to_numpy()
    lim=np.r_[0,np.flatnonzero(np.sign(np.diff(ll))[1:]!=np.sign(np.diff(ll))[:-1])+1,len(ll)-1]
    vals=[];weights=[]
    for lo,hi in zip(lim[:-1],lim[1:]):
        x=ll[lo:hi+1];y=lm[lo:hi+1]
        if x[-1]<x[0]:x=x[::-1];y=y[::-1]
        f=PchipInterpolator(x,y,extrapolate=False);v=f(lum)
        vals.append(v);weights.append(np.where(np.isfinite(v),abs(f.derivative()(lum)),0))
    vals=np.array(vals);weights=np.array(weights);total=weights.sum(0)
    ix=np.argmax(np.cumsum(weights,axis=0)>rng.random(len(lum))*total,axis=0)
    return np.where(total>0,10**vals[ix,np.arange(len(lum))],np.nan)

def run():
    with closing(sqlite3.connect(C/'cluster_catalogues.sqlite')) as con:
        # Remove every candidate source, including uncertain entries, so that
        # the observed cluster is not counted again underneath its simulation.
        known=set(pd.read_sql('SELECT gaia_dr3_id FROM objects WHERE gaia_dr3_id IS NOT NULL',con).gaia_dr3_id)
    iso=pd.read_csv(C/'cluster_isochrones_provisional.csv');end=pd.read_csv(V2/'structural_endpoints.csv')
    rows=[];receipts=[]
    for name in ['hyades_r2.5_assort_s43','pleiades_r3_assort_s43','pleiades_r4_fieldlow_s43']:
        start=time.monotonic();r=end[end.run==name].iloc[0];cl=r.cluster
        meta=json.loads((V2/'runs'/name/'initial.json').read_text());h,p=old.read_snapshot(r.snapshot)
        summary,p,rel,rad=old.structure(h,p,meta);pairs=bound_pairs(p,rel)
        vel=p.vel.copy()
        for a,b in pairs:vel[a]=vel[b]=(p.mass[a]*p.vel[a]+p.mass[b]*p.vel[b])/(p.mass[a]+p.mass[b])
        core=(rad<5)&(p.mass>.1)&(p.mass<2);vel-=np.median(vel[core],axis=0)
        ra,de,dist,angle,pm=old.sky_positions(rel,meta,vel)
        # Center-of-mass velocities avoid counting binary orbital speed as
        # the star's several-year astrometric proper motion.
        data=np.load(W/f'{cl}_cp_parent.npz');background=data['points'][~np.isin(data['source_id'],list(known))]
        err=[];mid=np.arange(8.5,22,.5)
        for mag in mid:
            ix=abs(data['G']-mag)<.35
            err.append([np.median(data[key][ix]) if ix.sum()>20 else np.nan for key in ['parallax_error','pmra_error','pmdec_error']])
        err=np.array(err)
        for col in range(3):
            ok=np.isfinite(err[:,col]);err[:,col]=np.interp(mid,mid[ok],err[ok,col])
        alive=(p.mass>.0099)&(p.star.type<10);suppress=~alive
        merged=[]
        for a,b in pairs:
            dra=((ra[a]-ra[b]+180)%360-180)*np.cos(np.deg2rad((de[a]+de[b])/2))
            if np.hypot(dra,de[a]-de[b])*3600>=1:continue
            if not (alive[a] and alive[b]):continue
            if p.mass[b]>p.mass[a]:a,b=b,a
            suppress[b]=True;merged.append((a,b))
        local_m10=m10(ra,de)
        for cloud in ['nc','f2']:
            q=iso[(iso.cluster==cl)&(iso.cloud==cloud)].sort_values('mass_Msun')
            lum=np.full(len(p.mass),-np.inf);low=alive&(p.mass<=.25)&(p.mass>=q.mass_Msun.min())
            lum[low]=PchipInterpolator(np.log10(q.mass_Msun),q.logL_Lsun)(np.log10(p.mass[low]))
            high=alive&(p.mass>.25)
            lum[high]=np.log10(np.maximum(p.star.lum[high],1e-30))
            gg=old.bandmag(lum,CAL[cl,'Gmag'])+5*np.log10(dist/10)
            for a,b in merged:
                lum[a]=np.logaddexp(lum[a]*np.log(10),lum[b]*np.log(10))/np.log(10)
                gg[a]=-2.5/np.log(10)*np.logaddexp(-.4*np.log(10)*gg[a],-.4*np.log(10)*gg[b])
            for noisy in [False,True]:
                rng=np.random.default_rng(583+(cloud=='f2'))
                for trial in range(24):
                    gm=gg+rng.normal(0,CAL[cl,'Gmag']['scatter'],len(gg))
                    detected=(rng.random(len(gg))<completeness(gm,local_m10))&~suppress
                    errors=np.array([np.interp(gm,mid,err[:,i]) for i in range(3)]).T
                    par=1000/dist;motions=pm.copy()
                    if noisy:
                        offset=rng.normal(size=(len(gm),3))*errors;par=par+offset[:,0];motions=motions+offset[:,1:]
                    quality=(par>2)&(1000*errors[:,0]/np.maximum(par,1e-8)**2<10)
                    pos=coordinates(ra,de,np.maximum(par,1e-8),motions[:,0],motions[:,1],cl)
                    eligible=detected&quality&window(pos,cl)
                    ix=np.flatnonzero(eligible);allpoints=np.vstack([background,pos[ix]])
                    k=neighbors(allpoints,pos[ix],cl,np.arange(len(ix))+len(background))
                    recovered=np.zeros(len(gm),bool);recovered[ix[k>=6]]=True
                    mass=inverse_mass(q,lum+rng.normal(0,.05,len(lum)),rng)
                    for population,mm in [('BD',(mass>=.03)&(mass<.08)),('reference',(mass>=.1)&(mass<.2))]:
                        for label,rr in [('all',np.ones(len(rad),bool)),('r0_10',rad<10),('r10_30',(rad>=10)&(rad<30)),('r30_100',(rad>=30)&(rad<100)),('r100plus',rad>=100)]:
                            target=mm&rr&~suppress
                            rows.append(dict(run=name,cluster=cl,birth_R=r.expected_birth_ratio,cloud=cloud,astrometric_noise=noisy,trial=trial,population=population,region=label,
                                             intrinsic_sources=int(target.sum()),gaia_detected=int((target&detected).sum()),automatic_window=int((target&eligible).sum()),CP_recovered=int((target&recovered).sum())))
        receipts.append(dict(run=name,background_count=len(background),removed_known_candidates=int(len(data['points'])-len(background)),wall_seconds=time.monotonic()-start,
                             core_mean_pm_masyr=np.median(pm[core],axis=0).tolist()))
        pd.DataFrame(rows).to_csv(W/'cp_injection_draws.csv',index=False)
        (W/'cp_injection_execution.json').write_text(json.dumps(receipts,indent=2)+'\n')
        print(receipts[-1],flush=True)
    frame=pd.DataFrame(rows)
    means=frame.groupby(['run','cluster','birth_R','cloud','astrometric_noise','population','region'])[['intrinsic_sources','gaia_detected','automatic_window','CP_recovered']].mean().reset_index()
    means['CP_recovery_given_Gaia']=means.CP_recovered/means.gaia_detected.replace(0,np.nan)
    means.to_csv(W/'cp_injection_summary.csv',index=False)

if __name__=='__main__':run()
