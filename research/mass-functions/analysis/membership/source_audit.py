"""Audit actual source eligibility without changing the retained candidate union."""
from pathlib import Path
from contextlib import closing
import json,sqlite3
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
W=Path(__file__).resolve().parent;C=W.parent/'cluster-mass-functions-oct9-v1';V1=W.parent/'cluster-dissolution-oct9-v1'

def main():
    with closing(sqlite3.connect(C/'cluster_catalogues.sqlite')) as con:
        s=pd.read_sql('SELECT object_id,source,selected,membership_probability,original FROM source_records',con)
    iso=pd.read_csv(C/'cluster_isochrones_provisional.csv');inventory={r['case']:r for r in json.loads((C/'track_inventory.json').read_text())}
    rows=[]
    for r in iso.itertuples():
        h=[json.loads(line) for line in (C/'tracks'/r.case/'evolution/history.jsonl').read_text().splitlines()]
        t=r.assumed_age_years-inventory[r.case]['initial_contraction_scale_years']
        teff=np.interp(t,[x['years'] for x in h],[x['Teff_K'] for x in h])
        rows.append(dict(cluster=r.cluster,cloud=r.cloud,mass_Msun=r.mass_Msun,logL_Lsun=r.logL_Lsun,Teff_K=teff,case=r.case))
    fulliso=pd.DataFrame(rows);fulliso.to_csv(W/'isochrones_with_temperature.csv',index=False)
    summary=[];objects=[]
    for cl in ['hyades','pleiades']:
        d=np.load(C/f'{cl}_mass_sensitivity_draws.npz');mass=d['mass_Msun'];ids=d['object_id']
        frame=pd.read_csv(V1/f'{cl}_spatial_selection.csv',low_memory=False).set_index('object_id').loc[ids]
        for kind,lo,hi in [('BD',.03,.08),('reference',.1,.2)]:frame['expected_'+kind]=((mass>=lo)&(mass<hi)).mean(0)
        ss=s[s.object_id.isin(ids)&s.selected.eq(1)]
        isboyle=frame.index.isin(ss[ss.source=='boyle2025'].object_id)
        hasother=frame.index.isin(ss[ss.source!='boyle2025'].object_id)
        frame['boyle_selected']=isboyle;frame['boyle_only']=isboyle&~hasother
        boyle={r.object_id:json.loads(r.original) for r in ss[ss.source=='boyle2025'].itertuples()}
        for name,key in [('boyle_Teff_K','Teff'),('boyle_CMD_consistent','FlagP'),('boyle_binary_flag','FlagB'),('boyle_probability','Prob')]:
            frame[name]=[boyle.get(i,{}).get(key,np.nan) for i in frame.index]
        for name,sel in [('all_associated',np.ones(len(frame),bool)),('boyle_present',isboyle),('boyle_only',isboyle&~hasother),
                         ('other_catalogues',hasother),('boyle_only_CMD_flag_false',isboyle&~hasother&(frame.boyle_CMD_consistent==0)),
                         ('outside_30deg',frame.angle_to_center_deg>=30),('outside_90deg',frame.angle_to_center_deg>=90)]:
            for kind in ['BD','reference']:
                summary.append(dict(cluster=cl,sample=name,kind=kind,expected_count=float(frame.loc[sel,'expected_'+kind].sum()),objects=int(np.sum(sel))))
        columns=['cluster','angle_to_center_deg','distance_to_center_pc','distance_for_luminosity_pc','logL_observed','logL_method','mass_median_Msun',
                 'expected_BD','expected_reference','boyle_selected','boyle_only','boyle_Teff_K','boyle_CMD_consistent','boyle_binary_flag','boyle_probability']
        objects.append(frame[columns].reset_index())
    pd.concat(objects).to_csv(W/'candidate_selection_audit.csv',index=False)
    pd.DataFrame(summary).to_csv(W/'catalogue_count_audit.csv',index=False)
    # Published group centers (Boyle et al., Table 1). Assign nearest centers
    # for a descriptive partition only: this is NOT a new membership classifier.
    centers={'Pleiades core':[-120.37,29.33,-33.12],'UPK 303':[-156.10,110.72,-43.89],
             'AB Dor':[-25.37,16.53,2.39],'HSC 1964':[-33.85,-106.80,-20.83],
             'UPK 545':[28.53,-320.95,-35.92],'GPC-1':[10.15,118.14,13.20],'GPC-2':[148.94,-101.76,-2.32]}
    b=pd.read_csv(C/'raw/pleiades_boyle2025_physical.csv',dtype={'GaiaID':str})
    xyz=b[['xp','yp','zp']].to_numpy();xyz[:,0]+=8122
    dist=np.linalg.norm(xyz[:,None,:]-np.array(list(centers.values()))[None,:,:],axis=2)
    b['nearest_center']=np.array(list(centers))[dist.argmin(1)];b['center_separation_pc']=dist.min(1)
    b['object_id']='pleiades:gaia:'+b.GaiaID
    allobj=pd.concat(objects).set_index('object_id')
    for k in ['expected_reference','expected_BD','boyle_only','distance_for_luminosity_pc']:
        b[k]=b.object_id.map(allobj[k])
    b.to_csv(W/'pleiades_group_partition.csv',index=False)
    group=b.groupby('nearest_center').agg(candidates=('GaiaID','size'),temperature_min_K=('Teff','min'),
            reference_stars=('expected_reference','sum'),BD_mass_assignments=('expected_BD','sum'),
            distance_median_pc=('distance_for_luminosity_pc','median'),CMD_fraction=('FlagP','mean'))
    group.to_csv(W/'pleiades_groups.csv')
    cuts=[]
    for cloud,g in fulliso[fulliso.cluster=='pleiades'].groupby('cloud'):
        m=np.geomspace(.03,.25,2000);te=PchipInterpolator(g.mass_Msun,g.Teff_K)(m)
        eligible=(te>3000)&(te<6500)
        cuts.append(dict(cloud=cloud,lowest_mass_temperature_cut_Msun=float(m[eligible].min()),
                         maximum_BD_temperature_K=float(te[m<.08].max()),
                         BD_temperature_eligibility_fraction=float(eligible[m<.08].mean())))
    (W/'temperature_cut.json').write_text(json.dumps(cuts,indent=2)+'\n')
    print(pd.DataFrame(summary).to_string(index=False,float_format=lambda x:f'{x:.4g}'))
    print(group.to_string(float_format=lambda x:f'{x:.4g}'));print(cuts)

if __name__=='__main__':main()
