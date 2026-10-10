"""Published automatic convergent-point selection, with explicit missing stages.

Implements Risbud et al. (2025) equations 2--10 and table 3. Manual rejection
of other overdensities and the Hyades SCCP extension are not reproduced.
"""
from pathlib import Path
from contextlib import closing
import json,sqlite3,time
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from astropy.coordinates import SkyCoord
from astropy import units as u
W=Path(__file__).resolve().parent;C=W.parent/'cluster-mass-functions-oct9-v1'
PARAMS={'hyades':dict(ac=97.3,dc=6.9,ra=66.701,dec=16.081,parallax=21.2,pmra=104.136,pmdec=-28.732,a=2.5,b=.8,r=15,zmax=35),
        'pleiades':dict(ac=89.8,dc=-44.7,ra=56.680,dec=24.108,parallax=7.38,pmra=19.955,pmdec=-45.457,a=1.5,b=.9,r=20,zmax=50)}

def coordinates(ra,dec,parallax,pmra,pmdec,cl):
    p=PARAMS[cl];ra=np.deg2rad(ra);de=np.deg2rad(dec)
    ac=np.deg2rad(p['ac']);dc=np.deg2rad(p['dc'])
    aa=np.sin(ac-ra);bb=np.cos(de)*np.tan(dc)-np.sin(de)*np.cos(ac-ra)
    norm=np.hypot(aa,bb);ss=aa/norm;cc=bb/norm
    vpar=4.74047*(ss*pmra+cc*pmdec)/parallax
    vper=4.74047*(-cc*pmra+ss*pmdec)/parallax
    mean=4.74047*np.hypot(p['pmra'],p['pmdec'])/p['parallax']
    co=SkyCoord(ra=ra*u.rad,dec=de*u.rad).galactic
    l=co.l.rad;b=co.b.rad;d=1000/np.asarray(parallax)
    xyz=np.column_stack([d*np.cos(b)*np.cos(l),d*np.cos(b)*np.sin(l),d*np.sin(b)])
    return np.column_stack([xyz,vpar-mean,vper])

def window(points,cl,interior=False):
    p=PARAMS[cl]
    cen=coordinates(np.array([p['ra']]),np.array([p['dec']]),np.array([p['parallax']]),np.array([p['pmra']]),np.array([p['pmdec']]),cl)[0,:3]
    r=p['r'] if interior else 0
    hemisphere=np.ones(len(points),bool) if cl=='hyades' else points[:,0]<-r
    return (hemisphere&(np.abs(points[:,2]-cen[2])<p['zmax']-r)&
            (abs(points[:,3])<50-(p['a'] if interior else 0))&
            (abs(points[:,4])<10-(p['b'] if interior else 0)))

def neighbors(points,query,cl,self_indices=None):
    p=PARAMS[cl];scale=np.array([p['r']]*3+[p['a'],p['b']])
    pos=points/scale;q=query/scale;tree=cKDTree(pos);answer=np.zeros(len(q),int)
    for start in range(0,len(q),1000):
        lists=tree.query_ball_point(q[start:start+1000],np.sqrt(2)+1e-12,workers=2)
        for k,ix in enumerate(lists):
            if not ix:continue
            delta=pos[ix]-q[start+k]
            good=(np.sum(delta[:,:3]**2,axis=1)<=1+1e-12)&(np.sum(delta[:,3:]**2,axis=1)<=1+1e-12)
            answer[start+k]=good.sum()
            if self_indices is not None and self_indices[start+k]>=0:answer[start+k]-=1
    return answer

def published(cl):
    with closing(sqlite3.connect(C/'cluster_catalogues.sqlite')) as con:
        records=con.execute("SELECT original FROM source_records WHERE source='risbud2025' AND cluster=?",(cl,)).fetchall()
    return pd.DataFrame([json.loads(r[0]) for r in records])

def main():
    checks=[]
    for cl in PARAMS:
        p=published(cl)
        pp=coordinates(*(p[k].astype(float).to_numpy() for k in ['RA_ICRS','DE_ICRS','plx','pmRA','pmDE']),cl)
        delta=pp[:,3:]-p[['Delvpll','vperobs']].astype(float).to_numpy()
        # Published CP coordinates are rounded to 0.1 degree, so small velocity
        # differences are expected; record rather than silently recalibrate.
        summary=dict(cluster=cl,published_members=len(p),velocity_residual_median_kms=np.median(delta,axis=0).tolist(),
                     velocity_residual_abs90_kms=np.quantile(abs(delta),.9,axis=0).tolist())
        path=W/'raw'/f'{cl}_gaia_cp_parent.csv'
        if not path.exists():checks.append(summary);continue
        start=time.monotonic();frame=pd.read_csv(path,dtype={'source_id':str})
        extra=W/'raw/hyades_gaia_cp_positive_x.csv'
        if cl=='hyades' and extra.exists():
            frame=pd.concat([frame,pd.read_csv(extra,dtype={'source_id':str})],ignore_index=True)
            assert frame.source_id.is_unique
        pts=coordinates(*(frame[k].to_numpy() for k in ['ra','dec','parallax','pmra','pmdec']),cl)
        eligible=window(pts,cl);target=window(pts,cl);index=np.flatnonzero(eligible)
        local=np.full(len(frame),-1,int);local[index]=np.arange(len(index))
        # For reproduction, query the published members rather than millions
        # of irrelevant field stars. The complete background remains in tree.
        lookup=pd.Series(np.arange(len(frame)),index=frame.source_id)
        found=p.GaiaDR3.map(lookup);valid=found.notna().to_numpy();fi=found[valid].astype(int).to_numpy()
        nn=neighbors(pts[index],pts[fi],cl,self_indices=local[fi])
        result=p[['GaiaDR3','5DNeighcount','Method']].copy()
        result['parent_present']=valid;result['automatic_neighbors']=np.nan;result['inside_automatic_window']=False
        result.loc[valid,'automatic_neighbors']=nn;result.loc[valid,'inside_automatic_window']=target[fi]
        result['automatic_selected']=result.inside_automatic_window&(result.automatic_neighbors>=6)
        # Save the more restrictive reading of the paper's boundary sentence;
        # compare both with its published members rather than assume exactness.
        result['eroded_boundary_selected']=False
        result.loc[valid,'eroded_boundary_selected']=window(pts[fi],cl,True)&(nn>=6)
        result.to_csv(W/f'{cl}_cp_reproduction.csv',index=False)
        summary.update(parent_rows=len(frame),parent_inside_window=int(eligible.sum()),
                       published_members_in_parent=int(valid.sum()),automatic_recovered_published=int(result.automatic_selected.sum()),
                       published_neighbors_median=float(pd.to_numeric(p['5DNeighcount']).median()),
                       reproduced_neighbors_median=float(np.median(nn)),
                       eroded_boundary_recovery=int(result.eroded_boundary_selected.sum()),
                       published_CP_members=int(p.Method.isin(['CP','CP, SCCP']).sum()),
                       automatic_recovered_CP_members=int((result.automatic_selected&p.Method.isin(['CP','CP, SCCP'])).sum()),
                       wall_seconds=time.monotonic()-start)
        np.savez_compressed(W/f'{cl}_cp_parent.npz',points=pts[index],source_id=frame.source_id.to_numpy(dtype=str)[index],
                            G=frame.phot_g_mean_mag.to_numpy()[index],parallax_error=frame.parallax_error.to_numpy()[index],
                            pmra_error=frame.pmra_error.to_numpy()[index],pmdec_error=frame.pmdec_error.to_numpy()[index])
        checks.append(summary);print(summary,flush=True)
    (W/'cp_reproduction.json').write_text(json.dumps(checks,indent=2)+'\n')
    # Simple neighborhood test checks two separate balls rather than a single
    # 5D ellipsoid; a point can sit inside one and outside the other.
    a,b,r=1.5,.9,20
    toy=np.array([[0,0,0,0,0],[r*.9,0,0,a*.9,0],[r*1.01,0,0,0,0],[0,0,0,a*1.01,0]])
    assert neighbors(toy,toy[:1],'pleiades',np.array([0]))[0]==1

if __name__=='__main__':main()
