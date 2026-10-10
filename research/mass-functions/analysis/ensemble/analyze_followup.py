"""Fixed-aperture structural checks and held-out BD predictions.

No IMF tuning, arbitrary population rescaling, or fitted completeness factor.
The model grid is a sensitivity ensemble, not a posterior distribution.
"""
from pathlib import Path
from contextlib import closing
import json, sqlite3
import numpy as np
import pandas as pd
import analysis_functions as a

W=Path(__file__).resolve().parent;V1=W.parent/'cluster-dissolution-oct9-v1'
C=W.parent/'cluster-mass-functions-oct9-v1'
PLAN=json.loads((W/'analysis_plan.json').read_text())

def observed():
    rows=[];profiles=[];drawrows=[];angular=[];coverage=[];mf=[]
    for cl in ['hyades','pleiades']:
        raw=np.load(C/f'{cl}_mass_sensitivity_draws.npz');mass=raw['mass_Msun'];ids=raw['object_id']
        d=pd.read_csv(V1/f'{cl}_spatial_selection.csv').set_index('object_id').loc[ids]
        if cl=='hyades' and (W/'hyades_2mass.csv').exists():
            phot=pd.read_csv(W/'hyades_2mass.csv',dtype={'source_id':str}).drop_duplicates('source_id').set_index('source_id')
            with closing(sqlite3.connect(C/'cluster_catalogues.sqlite')) as con:
                exact=pd.read_sql('SELECT object_id,gaia_dr3_id FROM objects',con).set_index('object_id').gaia_dr3_id.reindex(ids)
            for key,value in [('Jmag','j_m'),('Kmag','ks_m')]:
                extra=exact.map(phot[value]).to_numpy()
                d[key]=np.where(d[key].notna(),d[key],extra)
        r=d.distance_to_center_pc.to_numpy();ang=d.angle_to_center_deg.to_numpy()
        masks={'all_union':np.ones(len(d),bool),'r200_union':r<200,
            'central_matched':(r<10 if cl=='hyades' else ang<3)}
        if cl=='pleiades':
            with closing(sqlite3.connect(C/'cluster_catalogues.sqlite')) as con:
                rec=pd.read_sql("select object_id,selected,membership_probability from source_records where source='olivares2018_bouy2015'",con)
            rec=rec.groupby('object_id').agg({'selected':'max','membership_probability':'max'}).reindex(ids)
            masks['ir_matched']=(rec.selected.to_numpy()==1)&(ang<3)&(d.Kmag.to_numpy()<17)
            geometry=(rec.selected.to_numpy()==1)&(ang<3);missing=d.Kmag.isna().to_numpy()
        else:
            geometry=(r<18)&(d.ra_2016_deg.to_numpy()>=46)&(d.ra_2016_deg.to_numpy()<=85)&(d.dec_2016_deg.to_numpy()>=-4.5)&(d.dec_2016_deg.to_numpy()<=38.5)
            missing=d.Jmag.isna().to_numpy()
            masks['ir_matched']=geometry&(d.Jmag.to_numpy()<17.5)
        for pop,lo,hi in [('BD',.03,.08),('stars',.1,.2)]:
            member=(mass>=lo)&(mass<hi)
            coverage.append(dict(cluster=cl,population=pop,expected_in_geometry=float(member[:,geometry].sum(1).mean()),
                expected_missing_band=float(member[:,geometry&missing].sum(1).mean()),
                expected_selected=float(member[:,masks['ir_matched']].sum(1).mean())))
        # Source-specific overlap is diagnostic only. These lists are correlated
        # and must never be multiplied as independent membership evidence.
        for source in ['liu2025','risbud2025','kos2024']:
            with closing(sqlite3.connect(C/'cluster_catalogues.sqlite')) as con:
                selected=pd.read_sql('SELECT DISTINCT object_id FROM source_records WHERE cluster=? AND source=? AND selected=1',con,params=(cl,source)).object_id
            masks['source_'+source]=d.index.isin(selected)
        rng=np.random.default_rng(942+(cl=='pleiades'))
        for scenario,sel in masks.items():
            nb=((mass[:,sel]>=.03)&(mass[:,sel]<.08)).sum(1)
            ns=((mass[:,sel]>=.1)&(mass[:,sel]<.2)).sum(1)
            ratios=np.divide(nb,ns,out=np.full(len(nb),np.nan),where=ns>0)
            # Gamma-Poisson ratio draws add counting noise to each mass scenario.
            kk=rng.integers(len(nb),size=32768)
            counting=rng.gamma(nb[kk]+.5)/rng.gamma(ns[kk]+.5)
            qq=np.quantile(counting,[.16,.5,.84])
            rows.append(dict(cluster=cl,scenario=scenario,N_BD=np.median(nb),N_star=np.median(ns),
                ratio=np.nanmedian(ratios),ratio_lo=qq[0] if np.any(ns) else np.nan,ratio_hi=qq[2] if np.any(ns) else np.nan))
            for k,(b,s) in enumerate(zip(nb,ns)):
                drawrows.append(dict(cluster=cl,scenario=scenario,rep=k,N_BD=b,N_star=s,ratio=b/s if s else np.nan))
                counts=np.histogram(mass[k,sel],a.MASS_BINS)[0]
                for z,n in enumerate(counts):
                    mf.append(dict(cluster=cl,scenario=scenario,rep=k,bin=z,N=n,
                        mass_lo_Msun=a.MASS_BINS[z],mass_hi_Msun=a.MASS_BINS[z+1]))
            for pop,lo,hi in [('BD',.03,.08),('stars',.1,.2)]:
                for k,m in enumerate(mass):
                    hist=np.histogram(r[sel&(m>=lo)&(m<hi)],[0,5,10,20,50,200,np.inf])[0]
                    for b,n in enumerate(hist):profiles.append(dict(cluster=cl,scenario=scenario,rep=k,population=pop,bin=b,N=n))
                    abins=[0,5,10,20,40,80,180] if cl=='hyades' else [0,1,2,3,5,10,30,90,180]
                    hist=np.histogram(ang[sel&(m>=lo)&(m<hi)],abins)[0]
                    for b,n in enumerate(hist):angular.append(dict(cluster=cl,scenario=scenario,rep=k,population=pop,bin=b,angle_hi_deg=abins[b+1],N=n))
    pd.DataFrame(rows).to_csv(W/'observed_matched_counts.csv',index=False)
    pd.DataFrame(drawrows).to_csv(W/'observed_count_draws.csv',index=False)
    pd.DataFrame(profiles).to_csv(W/'observed_profile_draws.csv',index=False)
    pd.DataFrame(angular).to_csv(W/'observed_angular_profile_draws.csv',index=False)
    pd.DataFrame(coverage).to_csv(W/'matched_photometry_coverage.csv',index=False)
    pd.DataFrame(mf).to_csv(W/'observed_mass_function_draws.csv',index=False)

def models():
    cal=a.calibrations();ends=[];mocks=[]
    for base in [V1,W]:
        for path in sorted((base/'runs').glob('*')):
            if not (path/'initial.json').exists():continue
            meta=json.loads((path/'initial.json').read_text())
            if 'name' not in meta:continue
            snapshots=[p for d in [path,path/'continuation'] for p in d.glob('data.*') if p.name[5:].isdigit()]
            endpoint=None
            for f in sorted(snapshots,key=lambda p:(int(p.suffix[1:]),'continuation' in str(p)),reverse=True):
                try:
                    with f.open() as stream:head=stream.readline().split()
                    if len(head)<3 or abs(float(head[2])-meta['age_Myr'])>1e-6:continue
                    h,p=a.read_snapshot(f)
                    if len(p.mass)!=h.n:continue
                    endpoint=f;break
                except (ValueError,IndexError):continue
            if endpoint is None:continue
            row,p,rel,r=a.structure(h,p,meta)
            crit=PLAN['structural_screen'][meta['cluster']];ap=r<crit['aperture_pc']
            visible=~np.isin(p.star.type,[13,14,15]);luminous=p.star.type<10
            mm=p.mass[ap&visible];rr=r[ap&visible];order=np.argsort(rr)
            row['matched_aperture_pc']=crit['aperture_pc']
            row['aperture_mass_excluding_BH_NS_Msun']=float(mm.sum())
            row['aperture_luminous_mass_Msun']=float(p.mass[ap&luminous].sum())
            row['aperture_halfmass_pc']=float(np.interp(mm.sum()/2,np.cumsum(mm[order]),rr[order]))
            row['structural_screen_pass']=bool(crit['mass_excluding_BH_NS_Msun'][0]<=mm.sum()<=crit['mass_excluding_BH_NS_Msun'][1] and crit['half_mass_radius_pc'][0]<=row['aperture_halfmass_pc']<=crit['half_mass_radius_pc'][1])
            row['birth_mass_Msun']=meta['total_mass'];row['birth_ratio']=meta['birth_ratio'];row['workspace']=base.name
            row['expected_birth_ratio']=meta.get('expected_birth_ratio',.63)
            row['snapshot']=str(endpoint)
            mr,ledger=a.mock(p,rel,r,meta,cal);mocks+=mr;ends.append(row)
            ledger.to_csv(W/f'{meta["name"]}_descendants.csv',index=False)
            print(meta['name'],f'M={mm.sum():.4g}',f'rh={row["aperture_halfmass_pc"]:.4g}',row['structural_screen_pass'],flush=True)
    pd.DataFrame(ends).to_csv(W/'structural_endpoints.csv',index=False)
    pd.DataFrame(mocks).to_csv(W/'mock_survey_ratios.csv',index=False)
    pd.DataFrame(a.COUNT_DRAWS).to_csv(W/'model_count_draws.csv',index=False)
    pd.DataFrame(a.PROFILE_DRAWS).to_csv(W/'model_profile_draws.csv',index=False)
    pd.DataFrame(a.ANGULAR_PROFILE_DRAWS).to_csv(W/'model_angular_profile_draws.csv',index=False)
    pd.DataFrame(a.MASS_FUNCTION_DRAWS).to_csv(W/'model_mass_function_draws.csv',index=False)
    pd.DataFrame(a.MAGNITUDE_DRAWS).to_csv(W/'model_magnitude_draws.csv',index=False)

if __name__=='__main__':
    observed();models()
