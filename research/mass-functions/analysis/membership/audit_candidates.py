"""Apply the survey prerequisites to actual retained Gaia counterparts.

This diagnoses the candidate catalogue. It is not an independent completeness
sample: the candidates were found by heterogeneous, incomplete searches.
"""
from pathlib import Path
import json,sqlite3
import numpy as np
import pandas as pd
from cp_selection import coordinates,window,neighbors
W=Path(__file__).resolve().parent;C=W.parent/'cluster-mass-functions-oct9-v1'

def main():
    audit=pd.read_csv(W/'candidate_selection_audit.csv')
    with sqlite3.connect(C/'cluster_catalogues.sqlite') as con:
        obj=pd.read_sql('SELECT object_id,gaia_dr3_id FROM objects',con)
        sources=pd.read_sql("SELECT object_id,original FROM source_records WHERE source='risbud2025' AND selected=1",con)
    gaia=pd.read_csv(W/'candidate_gaia.csv',dtype={'source_id':str})
    audit=audit.merge(obj,on='object_id',validate='one_to_one').merge(gaia,left_on='gaia_dr3_id',right_on='source_id',how='left',validate='many_to_one')
    cpids=set(sources[sources.original.map(lambda x:json.loads(x)['Method'] in ['CP','CP, SCCP'])].object_id)
    audit['published_CP_selected']=audit.object_id.isin(cpids)
    output=[];counts=[];published=[]
    for cl in ['hyades','pleiades']:
        a=audit[audit.cluster==cl].copy().reset_index(drop=True)
        flags={};flags['All retained candidates']=np.ones(len(a),bool)
        flags['Gaia counterpart']=a.source_id.notna().to_numpy()
        astrometry=np.isfinite(a[['parallax','pmra','pmdec']].to_numpy()).all(1)&a.parallax.gt(2).to_numpy()
        photometry=np.isfinite(a[['phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag']].to_numpy()).all(1)
        quality=(a.ruwe.lt(1.4)&(a.astrometric_excess_noise.lt(1)|a.astrometric_excess_noise_sig.lt(2))).to_numpy()
        precise=(1000*a.parallax_error/a.parallax**2<10).to_numpy()
        flags['Astrometry and distance']=flags['Gaia counterpart']&astrometry
        flags['Three Gaia bands']=flags['Astrometry and distance']&photometry
        flags['Astrometric quality']=flags['Three Gaia bands']&quality
        flags['Distance precision']=flags['Astrometric quality']&precise
        pts=np.full((len(a),5),np.nan);ok=np.flatnonzero(astrometry)
        pts[ok]=coordinates(*(a[k].to_numpy()[ok] for k in ['ra','dec','parallax','pmra','pmdec']),cl)
        flags['Search volume']=flags['Distance precision']&window(pts,cl)
        data=np.load(W/f'{cl}_cp_parent.npz');index=pd.Series(np.arange(len(data['points'])),index=data['source_id'])
        ni=np.flatnonzero(flags['Search volume'])
        own=a.loc[ni,'source_id'].map(index).fillna(-1).to_numpy(dtype=int)
        nn=neighbors(data['points'],pts[ni],cl,own)
        flags['Automatic membership']=np.zeros(len(a),bool)
        flags['Automatic membership'][ni]=nn>=6
        a['automatic_neighbors']=np.nan;a.loc[ni,'automatic_neighbors']=nn
        for stage,(label,mask) in enumerate(flags.items()):
            a[f'selection_stage_{stage}']=mask
            for pop,col in [('BD','expected_BD'),('reference','expected_reference')]:
                counts.append(dict(cluster=cl,stage=stage,selection=label,population=pop,
                                   expected_count=float(a.loc[mask,col].sum()),source_count=int(mask.sum())))
        for pop,col in [('BD','expected_BD'),('reference','expected_reference')]:
            published.append(dict(cluster=cl,population=pop,published_CP_expected_count=float(a.loc[a.published_CP_selected,col].sum()),
                                  luminosity_usable_CP_objects=int(a.published_CP_selected.sum())))
        a['passes_BP_RP_requirement']=photometry;a['passes_quality_requirement']=quality
        output.append(a)
    pd.concat(output).to_csv(W/'candidate_recovery_details.csv',index=False)
    pd.DataFrame(counts).to_csv(W/'candidate_recovery_counts.csv',index=False)
    pd.DataFrame(published).to_csv(W/'published_CP_counts.csv',index=False)
    print(pd.DataFrame(counts).pivot(index=['cluster','stage','selection'],columns='population',values='expected_count').to_string(float_format=lambda x:f'{x:.4g}'))
    print(pd.DataFrame(published).to_string(index=False,float_format=lambda x:f'{x:.4g}'))

if __name__=='__main__':main()
