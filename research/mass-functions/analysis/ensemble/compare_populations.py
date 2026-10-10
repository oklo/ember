"""Stellar checks and conditional BD recovery, with no birth-IMF fitting."""
from pathlib import Path
from contextlib import closing
import json,sqlite3
import numpy as np
import pandas as pd
W=Path(__file__).resolve().parent

def main():
    obs=pd.read_csv(W/'observed_count_draws.csv')
    mod=pd.read_csv(W/'model_count_draws.csv')
    op=pd.read_csv(W/'observed_angular_profile_draws.csv')
    mp=pd.read_csv(W/'model_angular_profile_draws.csv')
    ends=pd.read_csv(W/'structural_endpoints.csv').set_index('run')
    rng=np.random.default_rng(104291);rows=[];radial=[]
    for (run,cloud,ext,scenario),d in mod.groupby(['run','cloud','faint_extension','scenario']):
        cl=d.cluster.iloc[0]
        # The below-grid light variants gave no material difference in v1.
        # Preserve both outputs, but do not count them as independent models.
        target='r200_union' if scenario.startswith('wide') else scenario
        o=obs[(obs.cluster==cl)&(obs.scenario==target)]
        if len(o)==0:continue
        # Published aggregate DANCe rates only: this gives an absolute-count
        # sensitivity check; it does not establish equal rates in each mass bin.
        classifier=.90/(1-.043) if cl=='pleiades' and scenario=='ir_matched' else 1.
        im=rng.integers(len(d),size=32768);io=rng.integers(len(o),size=32768)
        db=d.N_BD.to_numpy()[im]*classifier;ds=d.N_star.to_numpy()[im]*classifier
        ob=o.N_BD.to_numpy()[io];os=o.N_star.to_numpy()[io]
        obdraw=rng.gamma(ob+.5);osdraw=rng.gamma(os+.5)
        # Poisson draws represent finite simulated population and observation
        # counts. This is a conditional sensitivity interval, not the full grid
        # posterior or a calibrated mass-dependent selection estimate.
        mbdraw=rng.gamma(db+.5);msdraw=rng.gamma(ds+.5)
        relative=(obdraw/osdraw)/(mbdraw/msdraw)
        q=np.quantile(relative,[.16,.5,.84])
        means=d[['N_BD','N_star']].median()*classifier
        row=dict(run=run,cluster=cl,cloud=cloud,faint_extension=ext,scenario=scenario,
            structural_screen_pass=bool(ends.loc[run,'structural_screen_pass']),
            random_pair_stress=bool(ends.loc[run,'binary_mode']=='random' and ends.loc[run,'binary_fraction']>0),
            expected_birth_ratio=ends.loc[run,'expected_birth_ratio'],
            classifier_count_factor=classifier,model_N_BD=means.N_BD,model_N_star=means.N_star,
            observed_N_BD=o.N_BD.median(),observed_N_star=o.N_star.median(),
            model_ratio=np.median(d.ratio),observed_ratio=np.median(o.ratio),
            required_BD_recovery=o.N_BD.median()/means.N_BD,
            required_star_recovery=o.N_star.median()/means.N_star,
            relative_BD_recovery_lo=q[0],relative_BD_recovery=q[1],relative_BD_recovery_hi=q[2])
        pm=mp[(mp.run==run)&(mp.cloud==cloud)&(mp.faint_extension==ext)&(mp.scenario==scenario)&(mp.population=='stars')]
        po=op[(op.cluster==cl)&(op.scenario==target)&(op.population=='stars')]
        # Compare disjoint radial bins. Do not multiply correlated cumulative
        # aperture counts as though each were an independent observation.
        mm=pm.groupby('bin').N.mean()*classifier;oo=po.groupby('bin').N.mean()
        vm=pm.groupby('bin').N.var()*classifier**2;vo=po.groupby('bin').N.var()
        valid=np.ones(len(oo),bool)
        m=mm.loc[valid].to_numpy();n=oo.loc[valid].to_numpy()
        v=(vm.loc[valid]+vo.loc[valid]).to_numpy()+m+n
        row['stellar_count_chi2_diagnostic']=float(np.sum((n-m)**2/np.maximum(v,1)))
        row['stellar_CDF_max_difference']=float(np.max(abs(np.cumsum(n)/sum(n)-np.cumsum(m)/sum(m))))
        for k in mm.index:
            radial.append(dict(run=run,cluster=cl,cloud=cloud,faint_extension=ext,scenario=scenario,
                bin=k,model_stars=mm[k],observed_stars=oo[k],
                coordinate='angular_separation_deg',
                required_star_recovery=oo[k]/mm[k] if mm[k]>0 else np.nan))
        rows.append(row)
    result=pd.DataFrame(rows)
    result.to_csv(W/'joint_checks.csv',index=False)
    pd.DataFrame(radial).to_csv(W/'stellar_radial_checks.csv',index=False)
    summary=[]
    for (cl,sc,birth),g in result[(result.structural_screen_pass)&~result.random_pair_stress&(result.faint_extension!='continued_slope')].groupby(['cluster','scenario','expected_birth_ratio']):
        summary.append(dict(cluster=cl,scenario=sc,expected_birth_ratio=birth,n_models=g.run.nunique(),
            model_ratio_min=g.model_ratio.min(),model_ratio_max=g.model_ratio.max(),
            required_star_recovery_min=g.required_star_recovery.min(),required_star_recovery_max=g.required_star_recovery.max(),
            relative_BD_recovery_min=g.relative_BD_recovery.min(),relative_BD_recovery_max=g.relative_BD_recovery.max()))
    pd.DataFrame(summary).to_csv(W/'structurally_screened_ranges.csv',index=False)
    print(pd.DataFrame(summary).to_string(index=False,float_format=lambda x:f'{x:.4g}'))
    with closing(sqlite3.connect(W/'cluster_followup.sqlite')) as con:
        for name in ['structural_endpoints','mock_survey_ratios','observed_matched_counts',
                     'observed_count_draws','model_count_draws','joint_checks','stellar_radial_checks',
                     'structurally_screened_ranges']:
            pd.read_csv(W/f'{name}.csv').to_sql(name,con,if_exists='replace',index=False)

if __name__=='__main__':main()
