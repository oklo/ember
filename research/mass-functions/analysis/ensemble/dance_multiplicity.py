"""A catalogue-level binary sensitivity check, not a birth-MF correction.

Use published P_EMB to split source light, not published evolutionary masses.
Only the selected central DANCe sources are affected. Unequal companions not
recognized by the photometric classifier remain unknown. No missing sources
are invented and no correction is extended outside the observed footprint.
"""
from pathlib import Path
from contextlib import closing
import json,sqlite3
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
W=Path(__file__).resolve().parent;C=W.parent/'cluster-mass-functions-oct9-v1';V1=W.parent/'cluster-dissolution-oct9-v1'

def run():
    d=pd.read_csv(V1/'pleiades_spatial_selection.csv').set_index('object_id')
    with closing(sqlite3.connect(C/'cluster_catalogues.sqlite')) as con:
        rec=pd.read_sql("SELECT object_id, membership_probability, original FROM source_records WHERE source='olivares2018_bouy2015' AND selected=1",con)
    rec['pemb']=[json.loads(s)['pemb'] for s in rec.original]
    rec=rec.drop_duplicates('object_id').set_index('object_id')
    d=d.join(rec[['pemb','membership_probability']],how='inner')
    d=d[(d.angle_to_center_deg<3)&(d.Kmag<17)&d.logL_observed.notna()]
    iso=pd.read_csv(C/'cluster_isochrones_provisional.csv');rows=[]
    for cloud in ['nc','f2']:
        q=iso[(iso.cluster=='pleiades')&(iso.cloud==cloud)].sort_values('mass_Msun')
        # Every branch is retained, including any faint deuterium fold.
        mm=q.mass_Msun.to_numpy();ll=q.logL_Lsun.to_numpy()
        changes=np.flatnonzero(np.sign(np.diff(ll))[1:]!=np.sign(np.diff(ll))[:-1])+1
        limits=np.r_[0,changes,len(mm)-1]
        inverse=[]
        for lo,hi in zip(limits[:-1],limits[1:]):
            y=ll[lo:hi+1];x=np.log10(mm[lo:hi+1])
            if y[-1]<y[0]:x=x[::-1];y=y[::-1]
            inverse.append(PchipInterpolator(y,x,extrapolate=False))
        rng=np.random.default_rng(609+(cloud=='f2'))
        def invert(y):
            roots=np.array([f(y) for f in inverse]);weights=np.array([abs(f.derivative()(y)) for f in inverse])
            weights=np.where(np.isfinite(roots),weights,0);norm=weights.sum(0)
            ix=np.argmax(np.cumsum(weights,axis=0)>rng.random(len(y))*norm,axis=0)
            return np.where(norm>0,10**roots[ix,np.arange(len(y))],np.nan)
        error=np.hypot(d.e_logL_stat.fillna(.1).clip(lower=.01),d.logL_conversion_systematic_dex.fillna(.04))
        for rep in range(512):
            lum=d.logL_observed.to_numpy()+rng.normal(size=len(d))*error.to_numpy()+rng.normal(0,.04)
            binary=rng.random(len(d))<d.pemb.to_numpy()
            single=invert(lum)
            for flux in [0.,.6,1.]:
                if flux==0:mass=single
                else:
                    l1=lum[binary]-np.log10(1+flux);l2=l1+np.log10(flux)
                    mass=np.r_[single[~binary],invert(l1),invert(l2)]
                nb=np.sum((mass>=.03)&(mass<.08));ns=np.sum((mass>=.1)&(mass<.2))
                rows.append(dict(cloud=cloud,rep=rep,companion_flux_ratio=flux,N_BD=nb,N_star=ns,ratio=nb/ns,
                    recognized_binary_fraction_BD=np.mean(d.pemb.to_numpy()[(single>=.03)&(single<.08)]),
                    recognized_binary_fraction_star=np.mean(d.pemb.to_numpy()[(single>=.1)&(single<.2)])))
    out=pd.DataFrame(rows);out.to_csv(W/'dance_binary_sensitivity_draws.csv',index=False)
    summary=[]
    for (cloud,flux),g in out.groupby(['cloud','companion_flux_ratio']):
        qs=g.ratio.quantile([.16,.5,.84]).to_numpy()
        summary.append(dict(cloud=cloud,companion_flux_ratio=flux,ratio_lo=qs[0],ratio=qs[1],ratio_hi=qs[2],
            N_BD=g.N_BD.median(),N_star=g.N_star.median(),
            recognized_binary_fraction_BD=g.recognized_binary_fraction_BD.median(),
            recognized_binary_fraction_star=g.recognized_binary_fraction_star.median()))
    pd.DataFrame(summary).to_csv(W/'dance_binary_sensitivity.csv',index=False)
    print(pd.DataFrame(summary).to_string(index=False,float_format=lambda x:f'{x:.4g}'))
    (W/'dance_binary_scope.json').write_text(json.dumps(dict(
        source='https://arxiv.org/abs/1804.08360',objects=len(d),
        scope='Fixed nominal age/initial-age correction; two cloud cases and luminosity uncertainties. P_EMB recognizes high-mass-ratio pairs only. Flux-ratio assumptions 0.6 and 1 bracket a sensitivity, not a fitted q distribution. Not a correction for missing stars, low-q pairs or dissolution.',
        field_warning='This gives component counts in a central survey footprint, not a recovered birth component mass function.'),indent=2)+'\n')

if __name__=='__main__':run()
