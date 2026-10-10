"""Matched stellar training counts and separate brown-dwarf predictions."""
from contextlib import closing
import json,sqlite3,hashlib
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.interpolate import PchipInterpolator
from common import W,V1,V2,C,old,CAL,ANGLES,MASS
from selection import m10,completeness

def observed():
    receipts=[]
    for cl in ['hyades','pleiades']:
        raw=np.load(C/f'{cl}_mass_sensitivity_draws.npz');mass=raw['mass_Msun'];ids=raw['object_id']
        d=pd.read_csv(V1/f'{cl}_spatial_selection.csv').set_index('object_id').loc[ids]
        ang=d.angle_to_center_deg.to_numpy();sel=np.ones(len(d),bool)
        if cl=='pleiades':
            with closing(sqlite3.connect(C/'cluster_catalogues.sqlite')) as con:
                chosen=pd.read_sql("SELECT DISTINCT object_id FROM source_records WHERE source='olivares2018_bouy2015' AND selected=1",con).object_id
            sel=(ang>=3)|(d.index.isin(chosen)&(d.Kmag.to_numpy()<17))
        counts=np.array([np.histogram2d(ang[sel],m[sel],bins=[ANGLES[cl],MASS])[0] for m in mass])
        allcounts=np.array([np.histogram2d(ang,m,bins=[ANGLES[cl],MASS])[0] for m in mass])
        np.savez_compressed(W/f'observed_{cl}.npz',counts=counts,all_associated_counts=allcounts,angles=ANGLES[cl],mass_edges=MASS)
        receipts.append(dict(cluster=cl,all_associated_reference=float(allcounts[:,:,2:].sum((1,2)).mean()),
                             training_reference=float(counts[:,:,2:].sum((1,2)).mean()),
                             all_associated_BD=float(allcounts[:,:,0].sum(1).mean()),
                             comparison_BD=float(counts[:,:,0].sum(1).mean())))
    (W/'observed_scope.json').write_text(json.dumps(receipts,indent=2)+'\n')

def bound_pairs(p,rel):
    pairs=[]
    for i,j in cKDTree(rel).query_pairs(.01):
        r=np.linalg.norm(rel[i]-rel[j]);dv=np.linalg.norm(p.vel[i]-p.vel[j])
        if .5*dv*dv<old.G*(p.mass[i]+p.mass[j])/r:pairs.append((r,i,j))
    used=set();answer=[]
    for _,i,j in sorted(pairs):
        if i in used or j in used:continue
        used.update([i,j]);answer.append((i,j))
    return np.array(answer,dtype=int).reshape(-1,2)

def counts_for_population(mass,rel,pairs,meta,normalization=1.,replicates=32,types=None):
    cl=meta['cluster'];ra,dec,dist,ang,_=old.sky_positions(rel,meta)
    radius=np.linalg.norm(rel,axis=1)
    use=(mass>=.0099)&(mass<=.25)
    if types is not None:use&=types<10
    ids=np.flatnonzero(use);local=np.full(len(mass),-1,int);local[ids]=np.arange(len(ids))
    sky=(ra[ids]>=46)&(ra[ids]<=85)&(dec[ids]>=-4.5)&(dec[ids]<=38.5)&(radius[ids]<18)
    local_m10=m10(ra[ids],dec[ids]);dm=5*np.log10(dist[ids]/10)
    apparent_pairs=[];suppressed=np.zeros(len(ids),bool)
    if len(pairs):
        i,j=pairs.T
        dra=((ra[i]-ra[j]+180)%360-180)*np.cos(np.deg2rad((dec[i]+dec[j])/2))
        sep=np.hypot(dra,dec[i]-dec[j])*3600
        for ii,jj in pairs[sep<1]:
            a,b=local[ii],local[jj]
            if a<0 and b<0:continue
            if a<0 or b<0:
                other=ii if a<0 else jj;k=b if a<0 else a
                if mass[other]>.25 and (types is None or types[other]<10):suppressed[k]=True
                continue
            if mass[ids[b]]>mass[ids[a]]:a,b=b,a
            apparent_pairs.append((a,b));suppressed[b]=True
    pp=np.array(apparent_pairs,dtype=int).reshape(-1,2)
    iso=pd.read_csv(C/'cluster_isochrones_provisional.csv');result={}
    for cloud in ['nc','f2']:
        q=iso[(iso.cluster==cl)&(iso.cloud==cloud)].sort_values('mass_Msun')
        lm=np.log10(q.mass_Msun.to_numpy());ll=q.logL_Lsun.to_numpy()
        lum=PchipInterpolator(lm,ll,extrapolate=False)(np.log10(mass[ids]))
        lum=np.where(mass[ids]<q.mass_Msun.min(),-np.inf,lum)
        bands={b:old.bandmag(lum,CAL[cl,b])+dm for b in ['Gmag','Jmag','Kmag']}
        if len(pp):
            a,b=pp.T;lum[a]=np.logaddexp(lum[a]*np.log(10),lum[b]*np.log(10))/np.log(10)
            for value in bands.values():value[a]=-2.5/np.log(10)*np.logaddexp(-.4*np.log(10)*value[a],-.4*np.log(10)*value[b])
        limits=np.r_[0,np.flatnonzero(np.sign(np.diff(ll))[1:]!=np.sign(np.diff(ll))[:-1])+1,len(ll)-1]
        inverses=[]
        for lo,hi in zip(limits[:-1],limits[1:]):
            x=ll[lo:hi+1];y=lm[lo:hi+1]
            if x[-1]<x[0]:x=x[::-1];y=y[::-1]
            inverses.append(PchipInterpolator(x,y,extrapolate=False))
        rng=np.random.default_rng(765)
        arrays={(cloud,quality,deep):[] for quality in [.8,1.] for deep in [False,True]}
        for rep in range(replicates):
            gm=bands['Gmag']+rng.normal(0,CAL[cl,'Gmag']['scatter'],len(ids))
            jm=bands['Jmag']+rng.normal(0,CAL[cl,'Jmag']['scatter'],len(ids))
            km=bands['Kmag']+rng.normal(0,CAL[cl,'Kmag']['scatter'],len(ids))
            yy=lum+rng.normal(0,.05,len(ids));roots=np.array([f(yy) for f in inverses])
            weights=np.array([np.where(np.isfinite(root),abs(f.derivative()(yy)),0) for root,f in zip(roots,inverses)])
            total=weights.sum(0);pick=np.argmax(np.cumsum(weights,axis=0)>rng.random(len(ids))*total,axis=0)
            inferred=np.where(total>0,10**roots[pick,np.arange(len(ids))],np.nan)
            pg=completeness(gm,local_m10)
            for quality in [.8,1.]:
                for deep in [False,True]:
                    if cl=='hyades':ir=sky&(jm<(18 if deep else 17.5))
                    else:ir=(ang[ids]<(5 if deep else 3))&(km<(18 if deep else 17))
                    probability=1-(1-pg*quality)*(1-ir.astype(float))
                    if cl=='pleiades':
                        core=ang[ids]<3
                        probability[core]=(km[core]<17)*(.90/(1-.043))
                    probability[suppressed]=0
                    hist=np.histogram2d(ang[ids],inferred,bins=[ANGLES[cl],MASS],weights=probability)[0]/normalization
                    arrays[cloud,quality,deep].append(hist)
        result.update({key:np.array(value) for key,value in arrays.items()})
    return result

def save(name,arrays,meta):
    for (cloud,quality,deep),a in arrays.items():
        tag=f'{name}_{cloud}_q{quality:g}_'+('deep' if deep else 'shallow')
        target=W/'templates'/f'{tag}.npz'
        np.savez_compressed(target,counts=a)
        (W/'templates'/f'{tag}.json').write_text(json.dumps(dict(meta,tag=tag,cloud=cloud,membership_recovery=quality,deep=deep),indent=2)+'\n')

if __name__=='__main__':
    (W/'templates').mkdir(exist_ok=True);observed()
    endpoints=pd.read_csv(V2/'structural_endpoints.csv')
    for r in endpoints.itertuples():
        if r.binary_mode=='random' and r.binary_fraction>0:continue
        marker=W/'templates'/f'{r.run}.complete'
        if marker.exists():continue
        base=V1 if r.workspace==V1.name else V2
        meta=json.loads((base/'runs'/r.run/'initial.json').read_text())
        h,p=old.read_snapshot(r.snapshot);summary,p,rel,rad=old.structure(h,p,meta)
        a=counts_for_population(p.mass,rel,bound_pairs(p,rel),meta,types=p.star.type)
        save(r.run,a,dict(kind='bound',run=r.run,cluster=r.cluster,expected_birth_ratio=r.expected_birth_ratio,
                         birth_mass_Msun=r.birth_mass_Msun,structural_screen_pass=bool(r.structural_screen_pass),
                         aperture_mass_Msun=r.aperture_mass_excluding_BH_NS_Msun,aperture_halfmass_pc=r.aperture_halfmass_pc))
        marker.write_text(hashlib.sha256((W/'build_templates.py').read_bytes()).hexdigest()+'\n');print(r.run,flush=True)
    for path in sorted((W/'diffuse').glob('*.npz')):
        marker=W/'templates'/f'{path.stem}.complete'
        if marker.exists():continue
        info=json.loads(path.with_suffix('.json').read_text());d=np.load(path)
        a=counts_for_population(d['mass'],d['rel'],d['pairs'],info['meta'],normalization=info['N_reference_components'],replicates=16)
        aperture=9 if info['cluster']=='hyades' else 11.6
        fraction=float(d['mass'][np.linalg.norm(d['rel'],axis=1)<aperture].sum()/info['N_reference_components'])
        save(path.stem,a,dict(kind='diffuse',run=path.stem,cluster=info['cluster'],expected_birth_ratio=info['expected_birth_ratio'],
                             size_pc=info['size_pc'],speed_kms=info['speed_kms'],mass_per_reference_star=info['mass_per_reference_star'],
                             aperture_initial_mass_per_reference_upper=fraction))
        marker.write_text('complete\n');print(path.stem,flush=True)
