"""Fit reference stars only, freeze the choices, then inspect brown dwarfs.

Scores are covariance-weighted diagnostics, not calibrated likelihoods or
Bayes factors. No population scale is applied to an existing N-body cluster.
"""
import hashlib,json
import numpy as np
import pandas as pd
from common import W,ANGLES

G=.00430091727003628 # pc (km/s)^2 / Msun

def gravity_limit(bound,diffuse):
    r=diffuse['size_pc'];v=diffuse['speed_kms'];kin=1.5*v*v
    cluster=G*bound['birth_mass_Msun']*np.sqrt(2/np.pi)/r
    # Gaussian self-potential per unit mass = -G M/(2 sqrt(pi) sigma_x).
    mass=max(0,(.2*kin-cluster)*2*np.sqrt(np.pi)*r/G)
    return min(2000.,mass/diffuse['mass_per_reference_star']),cluster/kin

def main():
    allmeta=[json.loads(p.read_text()) for p in sorted((W/'templates').glob('*.json'))]
    bound=[r for r in allmeta if r['kind']=='bound'];diffuse=[r for r in allmeta if r['kind']=='diffuse']
    counts={r['tag']:np.load(W/'templates'/f'{r["tag"]}.npz')['counts'] for r in allmeta}
    observations={cl:np.load(W/f'observed_{cl}.npz')['counts'] for cl in ANGLES}
    rows=[]
    for b in bound:
        cl=b['cluster'];o=observations[cl][:,:,2:].reshape(len(observations[cl]),-1)
        bc=counts[b['tag']][:,:,2:].reshape(len(counts[b['tag']]),-1)
        target=o.mean(0);base=bc.mean(0)
        covariance=np.cov(o,rowvar=False)+np.cov(bc,rowvar=False)+np.diag(target+base+.5)
        inverse=np.linalg.inv(covariance);residual=target-base
        chi0=float(residual@inverse@residual)
        common=dict(bound_tag=b['tag'],run=b['run'],cluster=cl,expected_birth_ratio=b['expected_birth_ratio'],cloud=b['cloud'],
                    membership_recovery=b['membership_recovery'],deep=b['deep'],structural_screen_pass=b['structural_screen_pass'],
                    mass_Msun=b['aperture_mass_Msun'],radius_pc=b['aperture_halfmass_pc'],star_bins=len(target),
                    observed_stars=float(target.sum()),bound_stars=float(base.sum()),score_bound_only=chi0)
        rows.append(dict(common,diffuse_tag='',initial_diffuse_reference=0.,score_stars=chi0,
                         predicted_stars=float(base.sum()),initial_diffuse_mass_Msun=0.,neglected_gravity_fraction=0.,
                         diffuse_mass_inside_aperture_upper_Msun=0.))
        for d in diffuse:
            if any(b[k]!=d[k] for k in ['cluster','cloud','expected_birth_ratio','membership_recovery','deep']):continue
            vector=counts[d['tag']][:,:,2:].mean(0).ravel()
            cap,external=gravity_limit(b,d)
            if cap<=0:continue
            denominator=vector@inverse@vector
            amp=float(np.clip((vector@inverse@residual)/max(denominator,1e-30),0,cap))
            r=residual-amp*vector;score=float(r@inverse@r)
            m=amp*d['mass_per_reference_star'];gravity=external+G*m/(2*np.sqrt(np.pi)*d['size_pc']*1.5*d['speed_kms']**2)
            rows.append(dict(common,diffuse_tag=d['tag'],initial_diffuse_reference=amp,score_stars=score,
                             predicted_stars=float((base+amp*vector).sum()),initial_diffuse_mass_Msun=m,
                             neglected_gravity_fraction=gravity,diffuse_mass_inside_aperture_upper_Msun=amp*d['aperture_initial_mass_per_reference_upper']))
    fits=pd.DataFrame(rows);fits.to_csv(W/'stellar_fits.csv',index=False)
    # Choose separately under every declared birth ratio, cloud and selection
    # scenario. The brown-dwarf observations are absent from the objective.
    keys=['cluster','expected_birth_ratio','cloud','membership_recovery','deep']
    good=fits[fits.structural_screen_pass]
    best=good.loc[good.groupby(keys,sort=False).score_stars.idxmin()].copy()
    best.to_csv(W/'frozen_stellar_choices.csv',index=False)
    freeze=dict(training_mass_Msun=[.1,.2],selection_uses_BD_counts=False,
                fit_sha256=hashlib.sha256((W/'frozen_stellar_choices.csv').read_bytes()).hexdigest(),
                objective='Reference-star disjoint angular x mass cells only; covariance-weighted diagnostic. No fitted bound-cluster normalization.',
                caveat='Scores are not calibrated likelihoods. Published Gaia detection does not calibrate cluster membership or infrared completeness.')
    (W/'stellar_fit_receipt.json').write_text(json.dumps(freeze,indent=2)+'\n')
    # Prediction stage begins only after the stellar choices are written.
    pred=[];profiles=[];meta={r['tag']:r for r in allmeta}
    for row in best.itertuples():
        b=counts[row.bound_tag].mean(0);d=counts[row.diffuse_tag].mean(0) if row.diffuse_tag else np.zeros_like(b)
        predicted=b+row.initial_diffuse_reference*d;o=observations[row.cluster];obs=o.mean(0)
        nb=predicted[:,0].sum();ns=predicted[:,2:].sum();ob=obs[:,0].sum();os=obs[:,2:].sum()
        pred.append(dict(row._asdict(),predicted_BD=nb,observed_BD=ob,predicted_ratio=nb/ns,observed_ratio=ob/os,
                         residual_relative_BD_recovery=(ob/os)/(nb/ns)))
        for a in range(len(obs)):
            for m in range(obs.shape[1]):
                profiles.append(dict(cluster=row.cluster,birth=row.expected_birth_ratio,cloud=row.cloud,
                                     membership_recovery=row.membership_recovery,deep=row.deep,angle_bin=a,mass_bin=m,
                                     observed=obs[a,m],observed_mass_sd=o[:,a,m].std(),bound=b[a,m],diffuse=row.initial_diffuse_reference*d[a,m],predicted=predicted[a,m]))
    prediction=pd.DataFrame(pred);prediction.to_csv(W/'brown_dwarf_predictions.csv',index=False)
    pd.DataFrame(profiles).to_csv(W/'fitted_profiles.csv',index=False)
    print(prediction[['cluster','expected_birth_ratio','cloud','membership_recovery','deep','run','initial_diffuse_reference','score_stars','star_bins','predicted_stars','observed_stars','predicted_BD','observed_BD','predicted_ratio']].to_string(index=False,float_format=lambda x:f'{x:.4g}'))

if __name__=='__main__':main()
