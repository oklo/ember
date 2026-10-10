"""Observed ratios by aperture, preserving the full member union as the target."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
W=Path(__file__).resolve().parent;C=W.parent/'cluster-mass-functions-oct9-v1'
rows=[];coverage=[]
for cl,ra,dec,dist in [('hyades',66.75,15.87,47.03),('pleiades',56.75,24.12,136.2)]:
    a=np.load(C/f'{cl}_mass_sensitivity_draws.npz');m=a['mass_Msun']
    d=pd.read_csv(C/f'catalogues/{cl}_members_and_ember_masses.csv').set_index('object_id').loc[a['object_id']].copy()
    sky=SkyCoord(d.ra_2016_deg.to_numpy()*u.deg,d.dec_2016_deg.to_numpy()*u.deg)
    center=SkyCoord(ra*u.deg,dec*u.deg);ang=sky.separation(center).deg
    distance=d.distance_for_luminosity_pc.fillna(dist).to_numpy()
    radius=np.sqrt(distance**2+dist**2-2*distance*dist*np.cos(np.deg2rad(ang)))
    bd=(m>=.03)&(m<.08);star=(m>=.1)&(m<.2)
    for kind,limits,values in [('angular_deg',[2,3,5,10,20,180],ang),('radius_pc',[5,10,20,30,60,200,2000],radius)]:
        for limit in limits:
            mask=values<=limit;nb=bd[:,mask].sum(axis=1);ns=star[:,mask].sum(axis=1)
            ratio=np.divide(nb,ns,out=np.full_like(nb,np.nan,dtype=float),where=ns>0)
            q=np.nanquantile(ratio,[.16,.5,.84]);rows.append(dict(cluster=cl,aperture=kind,limit=limit,N_BD=np.median(nb),N_star=np.median(ns),ratio_lo=q[0],ratio=q[1],ratio_hi=q[2]))
    for kind,weight in [('BD',bd.mean(axis=0)),('stars',star.mean(axis=0))]:
        for source in sorted(set(';'.join(d.source_names.dropna()).split(';'))):
            present=d.source_names.fillna('').str.split(';').apply(lambda x:source in x).to_numpy()
            coverage.append(dict(cluster=cl,population=kind,source=source,expected_count=float(weight[present].sum())))
    d['angle_to_center_deg']=ang;d['distance_to_center_pc']=radius
    d['p_mass_BD']=bd.mean(axis=0);d['p_mass_reference_star']=star.mean(axis=0)
    d.to_csv(W/f'{cl}_spatial_selection.csv')
pd.DataFrame(rows).to_csv(W/'observed_aperture_ratios.csv',index=False)
pd.DataFrame(coverage).to_csv(W/'observed_source_contributions.csv',index=False)
(W/'spatial_audit_scope.json').write_text(json.dumps(dict(
    uncertainty='Only existing mass/age sensitivity draws; not counting, membership or completeness errors.',
    centers='Illustrative catalogue centers: Hyades (66.75,15.87 deg,47.03 pc); Pleiades (56.75,24.12 deg,136.2 pc).',
    caution='Distances without measurements may be membership-conditioned. Radius cuts are diagnostic only; all associated members remain in the target population.'),indent=2)+'\n')
print(pd.DataFrame(rows).to_string(index=False,float_format=lambda x:f'{x:.4g}'))
