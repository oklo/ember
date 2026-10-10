"""Plot, document and archive the completed stellar-count test."""
import hashlib,json,sqlite3,shutil
from contextlib import closing
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator,ScalarFormatter
from common import W,V2,ANGLES

ROOT=W.parent.parent
DOC=ROOT/'docs/research/independent'
NAME='CLUSTER_STELLAR_TEST_20261009'
COLORS={.63:'#276D8A',.55:'#BA7836'}

def observed_summary(cluster):
    counts=np.load(W/f'observed_{cluster}.npz')['counts']
    result={}
    for kind,selection in [('star',slice(2,None)),('BD',slice(0,1))]:
        v=counts[:,:,selection].sum(2);mean=v.mean(0)
        sigma=np.sqrt(mean+v.var(0))
        result[kind]=(mean,np.array([np.maximum(0,mean-sigma),mean+sigma]))
    return result

def plot():
    frame=pd.read_csv(W/'fitted_profiles.csv')
    plt.rcParams.update({'font.size':12,'axes.labelsize':12,'xtick.labelsize':10.5,
                         'ytick.labelsize':11,'axes.spines.top':False,'axes.spines.right':False,
                         'font.family':'DejaVu Sans','savefig.facecolor':'white'})
    fig,axs=plt.subplots(2,2,figsize=(12.5,8.6),sharex='col',gridspec_kw={'hspace':.18,'wspace':.21})
    for col,cluster in enumerate(['hyades','pleiades']):
        angles=ANGLES[cluster];x=np.arange(len(angles)-1);obs=observed_summary(cluster)
        for row,kind in enumerate(['star','BD']):
            ax=axs[row,col];part=frame[(frame.cluster==cluster)&((frame.mass_bin>=2) if kind=='star' else (frame.mass_bin==0))]
            for birth,g in part.groupby('birth'):
                curves=g.groupby(['cloud','membership_recovery','deep','angle_bin']).predicted.sum().unstack('angle_bin')
                representative=curves.loc[('f2',1.,False)].to_numpy()
                ax.fill_between(x,curves.min().to_numpy(),curves.max().to_numpy(),color=COLORS[birth],alpha=.19,linewidth=0)
                ax.plot(x,representative,color=COLORS[birth],lw=1.9,marker='o',ms=3.5)
            y,quantiles=obs[kind]
            errors=np.maximum(0,np.array([y-quantiles[0],quantiles[1]-y]))
            ax.errorbar(x,y,yerr=errors,fmt='o',ms=5,color='#222222',capsize=3,lw=1.05,zorder=8)
            ax.set_yscale('asinh',linear_width=1)
            ax.set_ylim(0,260 if row==0 else 175)
            ax.yaxis.set_major_locator(FixedLocator([0,1,3,10,30,100]))
            ax.yaxis.set_major_formatter(ScalarFormatter())
            ax.minorticks_off();ax.grid(axis='y',alpha=.17,lw=.6)
            ax.set_xlim(-.25,len(x)-.75)
            if col==0:ax.set_ylabel('Reference-star sources' if row==0 else 'Brown-dwarf sources')
            left_label=(col==0 and row==0)
            ax.text(.035 if left_label else .965,.955,('Fit: 0.10–0.20' if row==0 else 'Prediction: 0.03–0.08')+r' $M_\odot$',
                    transform=ax.transAxes,ha='left' if left_label else 'right',va='top',fontsize=11,
                    bbox={'facecolor':'white','edgecolor':'none','alpha':.85,'pad':1.5})
            if row==1:
                ax.set_xticks(x,[f'{a:g}–{min(b,180):g}' for a,b in zip(angles[:-1],angles[1:])])
                ax.set_xlabel('Angular separation from cluster center (degrees)')
        axs[0,col].text(.5,1.07,cluster.capitalize(),ha='center',transform=axs[0,col].transAxes,fontsize=16)
    handles=[Line2D([],[],color='#222222',marker='o',lw=0,label='Observed candidates'),
             Line2D([],[],color=COLORS[.63],lw=2,label='Model: birth R = 0.63'),
             Line2D([],[],color=COLORS[.55],lw=2,label='Model: birth R = 0.55')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.52,.975),ncol=3,frameon=False,fontsize=11.5)
    fig.subplots_adjust(left=.08,right=.985,top=.88,bottom=.15)
    fig.text(.08,.066,'R = N(0.03–0.08) / N(0.10–0.20) at birth. Models are chosen using stars alone. Shading: atmosphere and survey assumptions.',fontsize=10.5)
    fig.text(.08,.038,'Gaia source detection is modeled. Cluster-membership completeness remains uncertain; no birth mass function is inferred here.',fontsize=10.5)
    for suffix in ['png','pdf']:
        fig.savefig(DOC/f'{NAME}.{suffix}',dpi=220)
    plt.close(fig)

def report():
    p=pd.read_csv(W/'brown_dwarf_predictions.csv');check=pd.read_csv(W/'diagonal_objective_check.csv')
    rows=[]
    for (cl,birth),g in p.groupby(['cluster','expected_birth_ratio']):
        span=lambda key:f'{g[key].min():.4g}–{g[key].max():.4g}'
        rows.append(f'| {cl.capitalize()} | {birth:.2f} | {g.observed_stars.iloc[0]:.4g} | {span("predicted_stars")} | {g.observed_BD.iloc[0]:.4g} | {span("predicted_BD")} |')
    table='\n'.join(rows)
    cost=json.loads((V2/'validation_receipt.json').read_text())['cpu_hours_completed']
    wall=sum(r['wall_seconds'] for r in json.loads((W/'diffuse_execution.json').read_text()))
    text=f'''# Can the same models account for the stars and brown dwarfs?

**The restricted model family still fails the joint test.** Its best stellar
descriptions predict too many brown dwarfs, while also missing the observed
stellar spatial distribution. A smooth, diffuse associated population does not
repair both discrepancies. This does not establish different birth mass functions:
membership selection and the initial spatial structure remain insufficiently known.

![Stellar fits and brown-dwarf predictions]({NAME}.png)

Points show candidate-source counts inferred with Ember. Error bars combine
mass-inference variance with Poisson counting variance, showing an approximate
one-standard-deviation range clipped at zero; they
exclude unknown contamination and completeness. Curves show models chosen using
only the 0.10–0.20-solar-mass stars, with no adjustment to their brown-dwarf birth
ratio. Lines join angular bins for readability. The vertical scale is approximately
logarithmic above one source and linear close to zero. Shading spans the declared
atmosphere and survey scenarios, not a confidence interval. The 0.55 birth-ratio
alternative was simulated only for the Pleiades.

| Population | Birth R | Observed stars | Predicted stars | Observed brown dwarfs | Predicted brown dwarfs |
|---|---:|---:|---:|---:|---:|
{table}

R denotes N(0.03–0.08 solar masses)/N(0.10–0.20 solar masses), fixed comparison
intervals rather than a hydrogen-burning classification. Birth R counts components;
the observed and predicted counts above count sources after unresolved-binary
blending. Fractional counts are means over mass-inference draws. They can differ
slightly from medians in other summaries. The field-motivated birth ratios are
specified independently of these cluster counts; see the
[field analysis](FIELD_MASS_FUNCTION_20261009.md).

The Hyades model has about 110 reference stars within 10 degrees, compared with
87.54 candidates, but only 17.22 between 20 and 40 degrees, compared with 64.43.
Adding a positive diffuse population cannot fix every angular range at once.
The selected diffuse amplitude is zero in all eight survey/atmosphere scenarios.
The best Pleiades models add only 7.539–16.34 initial reference stars in the
diffuse population and continue to miss the extended stellar distribution.
Their covariance-weighted score improves by less than one; these scores are
diagnostics, not calibrated likelihoods or significance levels.

Removing cross-bin covariance from the fitting objective changes which Pleiades
realization is preferred and increases the diffuse amplitude. It still predicts
248.9–317.1 brown-dwarf sources against 195.7 observed, and does not recover the
angular pattern. The Hyades choice is unchanged. The conclusion therefore does
not depend on selecting one particular covariance weighting.

The distant Pleiades candidates illustrate why a single dissolving cluster is
too limited a starting description: essentially all candidates beyond 90 degrees
come from the Greater Pleiades Complex catalogue. These contribute 29.81 expected
reference stars. Of the 53.08 between 30 and 90 degrees, 52.05 appear in that
catalogue. Source lists overlap and these contributions must not be added across
catalogues. The published complex contains multiple coeval groups, plausibly
formed in one molecular cloud, rather than requiring every member to have escaped
one compact cluster. [Boyle, Bouma & Mann 2025](https://arxiv.org/abs/2511.07533).

## What is modeled

All 17 collisional simulations have reached their requested ages: 650 Myr for
the Hyades and 125 Myr for the Pleiades. The two random-pairing stress tests are
excluded from inference; 15 endpoints supply the bound-cluster alternatives.
Every simulated escaper remains available. Each existing cluster keeps its actual
population normalization. Broad final mass/size ranges screen candidates; these
are orientation ranges, not formal observational error bars. Dropping that screen
does not change the preferred models for the illustrated scenario.

The diffuse alternatives start at the cluster's Galactic birth location, with
one-axis Gaussian sizes of 10, 30 or 60 pc and one-axis speeds of 0.3, 0.7 or
1.5 km/s. Their trajectories are integrated through MWPotential2014 without
linearizing the Galactic force. The 27 combinations reuse existing mass and
binary draws; they are templates per initial reference star, not 27 additional
collisional clusters. This approximates an initially unbound association. It
does not simulate gas expulsion, primordial mass segregation or multiple distinct
birth centers.

Diffuse self-gravity and the bound cluster's force are omitted. Allowed amplitudes
must keep their estimated initial potential energy below 20% of kinetic energy.
Selected cases are below 2.200%, and their added mass in the central structural
aperture is below 0.009927 solar masses even using an initial-mass upper bound.
This check does not guarantee negligible deflection for every individual tracer.
A tighter orbit integration changes test positions by less than 0.000004831 pc.

Gaia detection now uses the published sky-dependent DR3 completeness function,
evaluated at each simulated source's location and magnitude. The original
order-7 M10 map, Zenodo record 8063930, is retained with its verified checksum.
This replaces a sharp Gaia magnitude cutoff, but supplies **source detection,
not recovery by a cluster-membership algorithm**.
[Cantat-Gaudin et al. 2023](https://arxiv.org/abs/2208.09335),
[GaiaUnlimited documentation](https://gaiaunlimited.readthedocs.io/en/latest/notebooks/dr3-empirical-completeness.html).

Gaia membership-recovery factors of 0.8 and 1 are sensitivity choices, not measured
completeness. Infrared selection uses the stated survey footprints and limiting
magnitudes, with a shallow/deep alternative. Its detailed color and proper-motion
selection is not reproduced. Within 3 degrees of the Pleiades, the comparison
uses DANCe candidates with probability above 0.84 and Ks < 17, and its published
aggregate recovery/contamination adjustment. Outside 3 degrees, all associated
and tail candidates are retained. The full candidate union remains in the data;
the matched central selection removes 10.87 expected reference stars and 5.258
expected brown dwarfs from this particular comparison.
[Olivares et al. 2018](https://arxiv.org/abs/1804.08360).

Ember supplies the existing cluster-composition isochrones and clear/cloud-clearing
alternatives. Mass inference retains nonmonotonic branches. The photometric
calibration, age, initial-state, reddening, and binary limitations remain those
documented in the [cluster mass-function analysis](CLUSTER_MASS_FUNCTIONS_20261009.md).
No new stellar tracks or changes to Ember physics were required.

## What this says about the apparent deficit

Survey reach remains an important part of the contrast, but the present models
do not explain the full data. Their brown-dwarf excess cannot yet be converted
into an intrinsic birth deficit because their stellar comparison already fails.
The next discriminating calculation is to measure recovery of simulated members
by the actual membership selections, and to represent the observed associated
groups rather than fitting a single smooth halo. More random seeds of the same
initial model would not resolve either limitation.

## Reproduction and data

The 17 collisional calculations used {cost:.4g} CPU-hours in total, including the
saved-state continuation accounting. They were reused here. The 27 diffuse orbit
integrations took {wall:.4g} seconds of measured wall time in total; this excludes
catalogue fitting and plotting and is not the cost of a new collisional ensemble.
No cluster jobs remain running.

Workspace: `out/cluster-stellar-fit-oct9-v1/`; reproduction order and dependencies
are in its `README.md`. `cluster_stellar_test.sqlite` contains the fit candidates,
frozen stellar choices, held-out predictions, angular/mass profiles, source
attribution and independent objective check. `plan.json`, map/source checksums,
orbit checks, and `validation_receipt.json` record the assumptions and checks.
The compact result tables are also saved beside this report.
'''
    (DOC/f'{NAME}.md').write_text(text)
    for source,suffix in [('brown_dwarf_predictions.csv','predictions'),('diagonal_objective_check.csv','objective_check'),('fitted_profiles.csv','profiles')]:
        shutil.copyfile(W/source,DOC/f'{NAME}_{suffix}.csv')

def archive():
    with closing(sqlite3.connect(W/'cluster_stellar_test.sqlite')) as con:
        for name in ['stellar_fits','frozen_stellar_choices','brown_dwarf_predictions','fitted_profiles','catalogue_attribution','diagonal_objective_check']:
            pd.read_csv(W/f'{name}.csv').to_sql(name,con,index=False,if_exists='replace')
        data=[]
        for path in sorted((W/'templates').glob('*.json')):data.append(json.loads(path.read_text()))
        pd.DataFrame(data).to_sql('template_metadata',con,index=False,if_exists='replace')
        metadata=[{'name':p.name,'json':p.read_text()} for p in W.glob('*.json')]
        pd.DataFrame(metadata).to_sql('receipts',con,index=False,if_exists='replace')
        con.execute('CREATE INDEX IF NOT EXISTS idx_fit_population ON stellar_fits(cluster,expected_birth_ratio,cloud,membership_recovery,deep)')
        con.execute('CREATE INDEX IF NOT EXISTS idx_profile_population ON fitted_profiles(cluster,birth,angle_bin,mass_bin)')
        assert con.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
        con.commit()

if __name__=='__main__':
    plot();report();archive()
    print(DOC/f'{NAME}.png')
