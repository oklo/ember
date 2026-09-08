# Low-mass equilibrium: published models and observed stars

Checked 2026-09-07 against the authors' original solar-metallicity
[BHAC15 tracks and internal structure](https://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/BHAC15_tracks%2Bstructure)
and [BCAH98 solar grid](https://perso.ens-lyon.fr/isabelle.baraffe/BCAH98_models.1).
Mass is exactly .100 Msun; published columns are linearly interpolated in
log age to 5 Gyr. Ember is the 4096-point static model documented in
[EQUILIBRIUM.md](EQUILIBRIUM.md), with no assigned evolutionary age.

| Quantity | Ember | BHAC15, 5 Gyr | Ember offset |
|---|---:|---:|---:|
| R/Rsun | .1286 | .1241 | +3.7% |
| L/Lsun | .0009756 | .0008715 | +12.0% |
| Teff (K) | 2844.6 | 2811 | +33.6 K (+1.2%) |
| Tc (K) | 4.566e6 | 4.581e6 | -.33% |
| rhoc (g/cm³) | 366.9 | 401.0 | -8.5% |

BHAC15 uses Rsun=6.96e10 cm and Lsun=3.839e33 erg/s; these are converted
to ember's nominal solar units before taking differences. The published
radius is rounded to .001 Rsun, so extra comparison digits are not accuracy.
The bracketing BHAC15 log ages are 9.689862 and 9.700961. Both have Teff=2811,
log L/Lsun=-3.061, R/Rsun=.124, log Tc=6.661 and log rhoc=2.6031.
The same grid changes very little from 1 to approximately 10 Gyr at this
mass: Teff=2810..2812 K, rounded radius=.124, log L/Lsun near -3.061.

BCAH98 at 5 Gyr gives Teff=2812 K, log g=5.251 and log L/Lsun=-3.07:
R≈.124 Rsun inferred from GM/g, and L≈8.51e-4 Lsun in its quoted solar units.
Its bracketing ages are 4.49349 and 5.03219 Gyr, with identical values of
these three columns. This is consistent with the newer grid at the
precision useful for the current experimental model.

The principal physical comparison limitations are the grey atmosphere and
EOS consistency defects in ember, plus unmatched composition. BHAC15 uses
interior Y=.28 and a revised solar mixture with Z=.0153, while ember uses
X=.7, Y=.28, Z=.02 with GS98 opacity and an explicit metals-as-helium EOS
approximation. BHAC15 couples its interior to a detailed atmosphere at
tau=100; ember currently uses a radiative Eddington grey boundary at tau=2/3.
See [Baraffe et al. (2015)](https://arxiv.org/abs/1503.04107) for their physics.

Ember also freezes He3=0. This is not an age-matched composition, but it
would be incorrect to impose full pp-chain equilibrium merely because a
.1 Msun star is on the main sequence: He3 destruction is slow in this
regime. See [Chabrier & Baraffe (1997)](https://arxiv.org/abs/astro-ph/9704118).
We have not quantified which approximation causes each offset. Close bulk
properties, including Tc, do not validate the local thermodynamic responses.

Downloaded source SHA-256 values:

- BHAC15: `b95474c5d4284373a2fed3f06d969a44bcd925ac0e5b226cc0235acb7e068d2a`
- BCAH98: `4dfd86699c60a5fb68bf863416f4ea2530d4ab26de9c8f50eca30d83dc542996`

The downloads are temporarily in `/tmp/ember-BHAC15-tracks.txt` and
`/tmp/ember-BCAH98-models1.txt`. These grids are comparison references,
not inputs to ember's stellar solve. No physics parameters were adjusted
to match them.

## Observational comparison

Additional 4096-point static solves use each star's adopted central mass,
with the same fixed composition and physics as the .100 Msun reference.
These are conditional predictions, not fits to age, composition, activity,
or the joint observational posterior. No physics parameters were calibrated
against these stars. Quoted observational uncertainties are the papers'
reported uncertainties; they are not independent in all columns.

| Star | Adopted M/Msun | Observed R/Rsun | Ember R/Rsun | Radius offset |
|---|---:|---:|---:|---:|
| EBLM J2114−39 B | .0993 ± .0033 | .1250 ± .0016 | .12787 | +2.3% |
| TRAPPIST-1 | .0898 ± .0023 | .1192 ± .0013 | .11769 | -1.3% |
| Proxima Centauri | .120 ± .003 | .146 ± .007 | .14900 | +2.1% |

| Star | Observed L/Lsun | Ember L/Lsun | Luminosity offset | Observed Teff (K) | Ember Teff (K) |
|---|---:|---:|---:|---:|---:|
| TRAPPIST-1 | .000553 ± .000019 | .0006379 | +15.4% | 2566 ± 26 | 2674 |
| Proxima Centauri | .00151 ± .00008 | .001681 | +11.3% | 2980 ± 80 | 3028 |

Sources and measurement dependencies:

- **EBLM J2114−39 B:** [Davis et al. (2024), published version](https://academic.oup.com/mnras/article/530/3/2565/7645099).
  Eclipse photometry and radial velocities are anchored to an external
  primary-star radius inferred using SED/parallax and stellar models. This
  is not a wholly model-independent absolute mass/radius determination.
  Use the final publication: the older preprint has different central
  values. The paper does not determine the secondary's Teff or bolometric
  luminosity; its spectroscopic temperature belongs to the primary.
- **TRAPPIST-1:** [Agol et al. (2021), Table 7](https://arxiv.org/html/2010.01074v2).
  Mass uses Mann et al.'s empirical K-band mass-luminosity relation; radius
  combines that mass with transit density. Luminosity is adopted from
  Ducrot et al. (2020); Teff follows from luminosity and radius, so those
  three quantities are correlated. This is not a dynamical absolute mass
  determination from the planetary system.
- **Proxima:** [Ribas et al. (2017)](https://arxiv.org/abs/1704.08449).
  Radius and Teff use the SED and interferometric angular-size constraint.
  The adopted mass is explicitly inferred by matching BHAC15 to L and
  Teff at 4.8 Gyr, so this row is not an independent test of the BHAC15
  mass-luminosity relation. Use the luminosity uncertainty .00008 from the
  abstract and summary table (the text has a decimal-place typo).

For an example with dynamical masses and interferometric radii,
[Kervella et al. (2016)](https://arxiv.org/abs/1607.04351) measure
GJ65 B at M=.1195 ± .0043 Msun, R=.159 ± .006 Rsun. They find its radius
12 ± 4% above BHAC15; component A is 14 ± 4% above. The authors propose
magnetic inhibition of convection in these rapid rotators. We have not
run an exact-mass ember model for GJ65 B. Its inflation demonstrates an
astrophysical source of disagreement, not a universal tolerance for
errors in nonmagnetic models.

Numerical checks for the three additional models:

- All converged at 4096 points with residuals below 2.4e-10. Increasing
  resolution from 2048 to 4096 changes R by <.01% and L by <.07%.
- Proxima uses tau_top=.002 because the shallower atmosphere encounters
  strict CMS19 support limits. At 2048 points, increasing tau_top from
  .002 to .004 changes R by .0012% and L by .0083%, much less than the
  observational offsets. The other two models retain tau_top=.001.
- Full solve outputs are `out/compare-eblm-j2114-fine.json`,
  `out/compare-trappist1-fine.json`, and `out/compare-proxima-fine.json`.
  [Compact numerical results and provenance](results/observational_comparison.json)
  record the adopted data, reproduction commands, and refinement checks.

The few-percent radius offsets and 10–15% luminosity excesses are much
larger than these numerical changes. They cannot yet be assigned to a
single physical approximation. Mass uncertainty, composition, age, and
covariance must enter any formal goodness-of-fit calculation. In
particular, do not divide the fixed-mass offsets by a single observed
error bar and report those ratios as a model rejection significance.
The approximately 24% local EOS Maxwell-identity defect in the .100 Msun
reference is a separate thermodynamic consistency diagnostic, not a 24%
uncertainty in radius or luminosity. Close radii do not remove the need
to repair that defect and replace the grey atmospheric boundary.
