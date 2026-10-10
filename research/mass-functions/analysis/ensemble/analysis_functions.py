"""Population accounting and deliberately bounded survey approximations.

This is not an injection/recovery completeness model of the catalogue union.
No fitted recovery factor is used. All descendants are retained in the ledger.
"""
from pathlib import Path
from contextlib import closing
import json,sys,sqlite3
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares
from astropy.coordinates import SkyCoord
from astropy import units as u
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014,rtide

W=Path(__file__).resolve().parent
sys.path.insert(0,str(W/'install/include'))
import petar
C=W.parent/'cluster-mass-functions-oct9-v1'
G=.00449830997959438;KMS=1.022712165045695
PROFILE_DRAWS=[]
COUNT_DRAWS=[]
ANGULAR_PROFILE_DRAWS=[]
MASS_FUNCTION_DRAWS=[]
MAGNITUDE_DRAWS=[]
MASS_BINS=np.array([.02,.03,.04,.05,.06,.08,.1,.125,.16,.2,.25])

def ratio(m,mask):
    b=np.count_nonzero(mask&(m>=.03)&(m<.08));s=np.count_nonzero(mask&(m>=.1)&(m<.2))
    return dict(N_BD=b,N_star=s,ratio=b/s if s else np.nan)

def center(pos,mass):
    # Density center, then shrinking spheres; excludes sparse distant escapers.
    ds,_=cKDTree(pos).query(pos,k=24)
    d=1/np.maximum(ds[:,-1],1e-8)**3
    x=np.average(pos,weights=d*d,axis=0)
    for _ in range(12):
        r=np.linalg.norm(pos-x,axis=1);sel=r<=np.quantile(r,.1)
        xnew=np.average(pos[sel],weights=mass[sel],axis=0)
        if np.linalg.norm(xnew-x)<1e-4:break
        x=xnew
    return x

def read_snapshot(path):
    h=petar.PeTarDataHeader(str(path),external_mode='galpy')
    p=petar.Particle(interrupt_mode='bse',external_mode='galpy')
    p.loadtxt(str(path),skiprows=1)
    return h,p

def structure(h,p,meta):
    good=p.mass>0;p=p[good]
    cen=center(p.pos,p.mass);rel=p.pos-cen;r=np.linalg.norm(rel,axis=1)
    absolute=cen+np.array(h.pos_offset);rr=np.hypot(*absolute[:2])/8000;zz=absolute[2]/8000
    Rj=20.
    for _ in range(30):
        mass=p.mass[r<Rj].sum()
        new=float(rtide(MWPotential2014,rr,zz,M=G*mass/(8000*(220*KMS)**2),use_physical=False))*8000
        if abs(new-Rj)<.001:break
        Rj=new
    mask=r<Rj;order=np.argsort(r[mask]);m=p.mass[mask][order]
    rh=float(np.interp(m.sum()/2,np.cumsum(m),r[mask][order])) if len(m) else np.nan
    vv=p.vel[mask];sig_raw=float(np.sqrt(np.mean(np.var(vv,axis=0)))/KMS) if len(vv) else np.nan
    # Internal binary speeds are not the velocity dispersion of the cluster.
    claimed=set();system_v=[];near=[]
    for i,j in cKDTree(rel).query_pairs(.01):
        dr=np.linalg.norm(rel[i]-rel[j]);dv=np.linalg.norm(p.vel[i]-p.vel[j])
        if .5*dv*dv<G*(p.mass[i]+p.mass[j])/dr:near.append((dr,i,j))
    for _,i,j in sorted(near):
        if i in claimed or j in claimed:continue
        claimed.update([i,j])
        if mask[i] and mask[j]:system_v.append((p.mass[i]*p.vel[i]+p.mass[j]*p.vel[j])/(p.mass[i]+p.mass[j]))
    system_v.extend(p.vel[i] for i in np.flatnonzero(mask) if i not in claimed)
    sig=float(np.sqrt(np.mean(np.var(system_v,axis=0)))/KMS) if len(system_v) else np.nan
    row=dict(run=meta['name'],cluster=meta['cluster'],time_Myr=h.time,initial_rh_pc=meta['rh'],seed=meta['seed'],binary_fraction=meta['binary_fraction'],binary_mode=meta.get('binary_mode','random'),N=len(p.mass),total_mass_Msun=p.mass.sum(),jacobi_radius_pc=Rj,mass_within_jacobi_Msun=mass,rhalf_within_jacobi_pc=rh,sigma_1d_kms=sig,sigma_1d_individual_kms=sig_raw,mass_within_9pc_Msun=p.mass[r<9].sum(),mass_within_11p6pc_Msun=p.mass[r<11.6].sum(),N_BH_within_jacobi=int(np.sum(mask&(p.star.type==14))),N_WD_within_jacobi=int(np.sum(mask&(p.star.type>=10)&(p.star.type<=12))),center_offset_from_orbit_pc=np.linalg.norm(absolute-np.array(meta['orbit']['present'][:3])))
    for key,sel in [('all',np.ones(len(r),bool)),('r10',r<10),('r18',r<18),('r30',r<30),('r200',r<200),('jacobi',mask)]:
        row.update({key+'_'+k:v for k,v in ratio(p.mass,sel).items()})
    return row,p,rel,r

def sky_positions(rel,meta,velocity=None):
    # Condition present cluster center on observed sky position; no change of
    # tail orientation, cluster sizes, distances relative to center, or velocities.
    pos=rel+np.array(meta['orbit']['present'][:3]);rho=np.hypot(pos[:,0],pos[:,1]);phi=np.arctan2(pos[:,1],pos[:,0])
    vel=np.broadcast_to(np.array(meta['orbit']['present'][3:]),pos.shape).copy()
    if velocity is not None:vel+=velocity
    vr=(pos[:,0]*vel[:,0]+pos[:,1]*vel[:,1])/rho/(220*KMS)
    vt=(pos[:,0]*vel[:,1]-pos[:,1]*vel[:,0])/rho/(220*KMS)
    o=Orbit(np.column_stack([rho/8000,vr,vt,pos[:,2]/8000,vel[:,2]/(220*KMS),phi]),ro=8,vo=220,zo=.0208,solarmotion='schoenrich')
    ra=o.ra();dec=o.dec();dist=o.dist()*1000
    coord=SkyCoord(ra*u.deg,dec*u.deg);c=SkyCoord(meta['ra']*u.deg,meta['dec']*u.deg)
    return ra,dec,dist,coord.separation(c).deg,np.column_stack([o.pmra(),o.pmdec()])

def calibrations():
    ret={};receipt=[]
    for cl in ['hyades','pleiades']:
        d=pd.read_csv(C/f'catalogues/{cl}_members_and_ember_masses.csv')
        filled={}
        if cl=='hyades' and (W/'hyades_2mass.csv').exists():
            phot=pd.read_csv(W/'hyades_2mass.csv',dtype={'source_id':str}).drop_duplicates('source_id').set_index('source_id')
            with closing(sqlite3.connect(C/'cluster_catalogues.sqlite')) as con:
                exact=pd.read_sql('SELECT object_id,gaia_dr3_id FROM objects',con).set_index('object_id').gaia_dr3_id
            for band,column in [('Jmag','j_m'),('Kmag','ks_m')]:
                extra=d.object_id.map(exact).map(phot[column])
                filled[band]=int((d[band].isna()&extra.notna()).sum())
                d[band]=d[band].fillna(extra)
        for band in ['Gmag','Jmag','Kmag']:
            q=d[['logL_observed',band,'distance_for_luminosity_pc']].dropna()
            q=q[(q.logL_observed< -1.65)&(q.logL_observed> -5.4)&(q.distance_for_luminosity_pc>0)]
            x=q.logL_observed.to_numpy();y=q[band].to_numpy()-5*np.log10(q.distance_for_luminosity_pc.to_numpy()/10)
            order=np.argsort(x);groups=np.array_split(order,min(12,len(order)//4))
            if cl=='hyades' and band in ['Jmag','Kmag']:
                # Keep the sparse, faint calibrators represented even when
                # newly recovered stellar photometry outnumbers them greatly.
                edges=np.arange(-5.5,-1.49,.25)
                groups=[np.flatnonzero((x>=lo)&(x<hi)) for lo,hi in zip(edges[:-1],edges[1:])]
                groups=[g for g in groups if len(g)>=2]
            xx=np.array([np.median(x[g]) for g in groups]);yy=np.array([np.median(y[g]) for g in groups])
            fit=least_squares(lambda a:np.polynomial.polynomial.polyval(xx+3,a)-yy,[12.,-3.,0.],loss='soft_l1',f_scale=.2).x
            scatter=1.4826*np.median(abs(y-np.polynomial.polynomial.polyval(x+3,fit)))
            ret[cl,band]=dict(coefficients=fit,xmin=float(xx.min()),xmax=float(xx.max()),scatter=max(.12,float(scatter)))
            receipt.append(dict(cluster=cl,band=band,n=len(q),new_2mass_measurements=filled.get(band,0),logL_min=xx.min(),logL_max=xx.max(),coefficients=fit.tolist(),robust_scatter_mag=scatter))
    (W/'photometric_calibration.json').write_text(json.dumps(receipt,indent=2))
    return ret

def bandmag(logL,c):
    # Linear tangent continuation outside the empirical interval. In most cases
    # the only long continuation is for Gaia-faint BDs that are rejected anyway.
    x=np.clip(logL,c['xmin'],c['xmax']);y=np.polynomial.polynomial.polyval(x+3,c['coefficients'])
    slope=np.clip(c['coefficients'][1]+2*c['coefficients'][2]*(x+3),-8.,-2.)
    return y+slope*(logL-x)

def mock(p,rel,r,meta,cal):
    cl=meta['cluster'];corev=np.median(p.vel[r<10],axis=0)
    ra,dec,dist,angle,pm=sky_positions(rel,meta,p.vel-corev)
    _,_,_,_,bulkpm=sky_positions(rel,meta)
    cosine=np.sum(pm*bulkpm,axis=1)/(np.linalg.norm(pm,axis=1)*np.linalg.norm(bulkpm,axis=1))
    pmangle=np.degrees(np.arccos(np.clip(cosine,-1,1)))
    iso=pd.read_csv(C/'cluster_isochrones_provisional.csv');rows=[];rng=np.random.default_rng(889)
    use=(p.mass>=.0099)&(p.mass<=.25);mass=p.mass[use];radius=r[use];ang=angle[use]
    phot={}
    variants=[(c,e) for c in ['nc','f2'] for e in (['dark_below_grid','continued_slope'] if cl=='hyades' else ['covered'])]
    for cloud,faint_extension in variants:
        rng=np.random.default_rng(889+(cloud=='f2'))
        q=iso[(iso.cluster==cl)&(iso.cloud==cloud)].sort_values('mass_Msun')
        forward=PchipInterpolator(np.log10(q.mass_Msun),q.logL_Lsun,extrapolate=False)
        logL=forward(np.log10(mass));dm=5*np.log10(dist[use]/10)
        missing=mass<q.mass_Msun.iloc[0]
        if faint_extension=='continued_slope':
            slope=(q.logL_Lsun.iloc[1]-q.logL_Lsun.iloc[0])/np.log10(q.mass_Msun.iloc[1]/q.mass_Msun.iloc[0])
            logL[missing]=q.logL_Lsun.iloc[0]+slope*np.log10(mass[missing]/q.mass_Msun.iloc[0])
        else:logL[missing]=-np.inf
        for band in ['Gmag','Jmag','Kmag']:phot[band]=bandmag(logL,cal[cl,band])+dm
        # Dynamically surviving unresolved pairs, including primordial pairs.
        # Bind-energy test rejects chance projected alignments. Mock angular
        # resolution is 1 arcsec, bracketed separately from membership losses.
        pairs=[]
        for i,j in cKDTree(rel).query_pairs(.01):
            dr=np.linalg.norm(rel[i]-rel[j]);dv=np.linalg.norm(p.vel[i]-p.vel[j])
            if .5*dv*dv-G*(p.mass[i]+p.mass[j])/dr>=0:continue
            a=SkyCoord(ra[i]*u.deg,dec[i]*u.deg).separation(SkyCoord(ra[j]*u.deg,dec[j]*u.deg)).arcsec
            if a<1:pairs.append((i,j,dr))
        pairs.sort(key=lambda q:q[2]);claimed=set();blend={};local={i:j for j,i in enumerate(np.flatnonzero(use))}
        suppress=np.zeros(len(mass),bool)
        brightening=np.zeros(len(mass))
        # Two calibrated low-mass components: sum the Ember and band fluxes.
        # A low-mass component paired to >.25 Msun is hidden from both bins.
        for i,j,_ in pairs:
            if i in claimed or j in claimed:continue
            claimed.update([i,j])
            if i not in local and j not in local:continue
            if i not in local or j not in local:
                k=local.get(i,local.get(j));other=j if i in local else i
                if p.mass[other]>.25 and p.star.type[other]<10:suppress[k]=True
                continue
            a,b=local[i],local[j]
            if mass[b]>mass[a]:a,b=b,a
            primary_lum=logL[a]
            logL[a]=np.logaddexp(logL[a]*np.log(10),logL[b]*np.log(10))/np.log(10)
            brightening[a]=2.5*(logL[a]-primary_lum) if np.isfinite(primary_lum) else 0.
            for band in phot:phot[band][a]=-2.5/np.log(10)*np.logaddexp(-.4*np.log(10)*phot[band][a],-.4*np.log(10)*phot[band][b])
            suppress[b]=True
        # Recover source-equivalent masses through the same branch-preserving
        # PCHIP inversion as the observational pipeline, with 0.05 dex scatter.
        mm=q.mass_Msun.to_numpy();ll=q.logL_Lsun.to_numpy();changes=np.flatnonzero(np.sign(np.diff(ll))[1:]!=np.sign(np.diff(ll))[:-1])+1
        limits=np.r_[0,changes,len(mm)-1]
        def invert(obs):
            roots=[];weights=[]
            for lo,hi in zip(limits[:-1],limits[1:]):
                y=ll[lo:hi+1];x=np.log10(mm[lo:hi+1])
                if y[-1]<y[0]:y=y[::-1];x=x[::-1]
                it=PchipInterpolator(y,x,extrapolate=False);root=it(obs);wt=np.where(np.isfinite(root),abs(it.derivative()(obs)),0.)
                roots.append(root);weights.append(wt)
            roots=np.array(roots);weights=np.array(weights);norm=weights.sum(axis=0)
            ix=np.argmax(np.cumsum(weights,axis=0)>rng.random(len(obs))*norm,axis=0)
            return np.where(norm>0,10**roots[ix,np.arange(len(obs))],np.nan)
        for scenario in ['central_matched','wide_union_shallow','wide_union_deep','ir_matched']:
            samples=[]
            for rep in range(64):
                # Survey thresholds are sensitivity brackets, not calibrated
                # completeness curves for the heterogeneous published union.
                gm=phot['Gmag']+rng.normal(0,cal[cl,'Gmag']['scatter'],len(mass))
                jm=phot['Jmag']+rng.normal(0,cal[cl,'Jmag']['scatter'],len(mass))
                km=phot['Kmag']+rng.normal(0,cal[cl,'Kmag']['scatter'],len(mass))
                if scenario=='ir_matched':
                    # Exact angular/magnitude limits for the DANCe core check.
                    # Hyades uses a common IR aperture, without claiming that
                    # the heterogeneous candidate union is complete there.
                    if cl=='pleiades':selected=(ang<3)&(km<17)
                    else:
                        sky=(ra[use]>=46)&(ra[use]<=85)&(dec[use]>=-4.5)&(dec[use]<=38.5)
                        selected=sky&(radius<18)&(jm<17.5)
                elif scenario=='central_matched':
                    footprint=radius<10 if cl=='hyades' else ang<3
                    detect=(gm<20.3)|((jm<18) if cl=='hyades' else (km<17))
                    selected=footprint&detect
                else:
                    deep=scenario=='wide_union_deep'
                    gaia=(radius<200)&(gm<(20.7 if deep else 20.0))
                    if cl=='hyades':
                        sky=(ra[use]>=46)&(ra[use]<=85)&(dec[use]>=-4.5)&(dec[use]<=38.5)
                        kin=(pmangle[use]<7) if deep else ((pm[use,0]>42)&(pm[use,0]<197)&(pm[use,1]>-92)&(pm[use,1]<43))
                        ir=sky&(radius<18)&kin&(jm<(18 if deep else 17.5))
                    else:
                        ir=(ang<(5 if deep else 3))&(km<(18 if deep else 17))
                    selected=gaia|ir
                selected &= ~suppress
                if cl=='pleiades' and scenario=='ir_matched':
                    edges=np.arange(13.5,17.01,.5)
                    hist=np.histogram(km[selected],edges)[0]
                    for k,n in enumerate(hist):
                        MAGNITUDE_DRAWS.append(dict(run=meta['name'],cloud=cloud,rep=rep,
                            K_lo=edges[k],K_hi=edges[k+1],N=n))
                inferred=invert(logL+rng.normal(0,.05,len(mass)))
                result=ratio(inferred,selected)
                samples.append(result)
                hist=np.histogram(inferred[selected],MASS_BINS)[0]
                for k,n in enumerate(hist):
                    MASS_FUNCTION_DRAWS.append(dict(run=meta['name'],cluster=cl,cloud=cloud,
                        faint_extension=faint_extension,scenario=scenario,rep=rep,bin=k,N=n,
                        mass_lo_Msun=MASS_BINS[k],mass_hi_Msun=MASS_BINS[k+1]))
                COUNT_DRAWS.append(dict(run=meta['name'],cluster=cl,cloud=cloud,
                    faint_extension=faint_extension,scenario=scenario,rep=rep,**result,
                    N_bright_binary_BD=int(np.sum(selected&(inferred>=.03)&(inferred<.08)&(brightening>.5))),
                    N_bright_binary_star=int(np.sum(selected&(inferred>=.1)&(inferred<.2)&(brightening>.5)))))
                rbins=np.array([0,5,10,20,50,200,np.inf])
                for pop,lo,hi in [('BD',.03,.08),('stars',.1,.2)]:
                    hist=np.histogram(radius[selected&(inferred>=lo)&(inferred<hi)],rbins)[0]
                    for k,n in enumerate(hist):
                        PROFILE_DRAWS.append(dict(run=meta['name'],cluster=cl,cloud=cloud,
                            faint_extension=faint_extension,scenario=scenario,rep=rep,
                            population=pop,bin=k,rlo_pc=rbins[k],rhi_pc=rbins[k+1],N=n))
                    abins=np.array([0,5,10,20,40,80,180] if cl=='hyades' else [0,1,2,3,5,10,30,90,180])
                    hist=np.histogram(ang[selected&(inferred>=lo)&(inferred<hi)],abins)[0]
                    for k,n in enumerate(hist):
                        ANGULAR_PROFILE_DRAWS.append(dict(run=meta['name'],cluster=cl,cloud=cloud,
                            faint_extension=faint_extension,scenario=scenario,rep=rep,
                            population=pop,bin=k,angle_hi_deg=abins[k+1],N=n))
            vals=np.array([s['ratio'] for s in samples]);nb=np.array([s['N_BD'] for s in samples]);ns=np.array([s['N_star'] for s in samples])
            rows.append(dict(run=meta['name'],cluster=cl,cloud=cloud,faint_extension=faint_extension,scenario=scenario,unresolved_pairs=len(claimed)//2,suppressed_lowmass_components=int(suppress.sum()),ratio=float(np.nanmedian(vals)),ratio_q16=float(np.nanquantile(vals,.16)),ratio_q84=float(np.nanquantile(vals,.84)),N_BD=np.median(nb),N_star=np.median(ns)))
    ledger=pd.DataFrame(dict(id=p.id,mass_Msun=p.mass,r_pc=r,ra_deg=ra,dec_deg=dec,distance_pc=dist,angle_deg=angle,pmra_masyr=pm[:,0],pmdec_masyr=pm[:,1],pm_direction_offset_deg=pmangle,x_pc=rel[:,0],y_pc=rel[:,1],z_pc=rel[:,2]))
    return rows,ledger

def main():
    cal=calibrations();allrows=[];mocks=[];ends=[]
    for path in sorted((W/'runs').glob('*_r*_s*')):
        if not (path/'initial.json').exists():continue
        meta=json.loads((path/'initial.json').read_text());snapshots=sorted([f for f in path.glob('data.*') if f.name[5:].isdigit()],key=lambda f:int(f.suffix[1:]))
        if (path/'continuation.json').exists():
            snapshots=[f for f in snapshots if int(f.suffix[1:])<=6]
            more=[f for f in (path/'continuation').glob('data.*') if f.name[5:].isdigit()]
            # Each resumed snapshot supersedes the copy at the same file index.
            byid={int(f.suffix[1:]):f for f in snapshots+more};snapshots=[byid[k] for k in sorted(byid)]
        for f in snapshots:
            try:h,p=read_snapshot(f)
            except (ValueError,IndexError):continue # potentially currently writing
            if len(p.mass)!=h.n:continue
            row,p,rel,r=structure(h,p,meta);allrows.append(row)
            if abs(h.time-meta['age_Myr'])<1e-6:
                mr,ledger=mock(p,rel,r,meta,cal);mocks+=mr;ends.append(row);ledger.to_csv(path/'final_descendants.csv',index=False)
    pd.DataFrame(allrows).to_csv(W/'dynamical_history.csv',index=False)
    pd.DataFrame(ends).to_csv(W/'dynamical_endpoints.csv',index=False)
    pd.DataFrame(mocks).to_csv(W/'mock_survey_ratios.csv',index=False)
    print(pd.DataFrame(ends).to_string(index=False,float_format=lambda x:f'{x:.4g}'))
    if mocks:print(pd.DataFrame(mocks).to_string(index=False,float_format=lambda x:f'{x:.4g}'))

if __name__=='__main__':main()
