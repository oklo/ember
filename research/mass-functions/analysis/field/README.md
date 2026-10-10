# Ember field mass function

Corrected nearby present-day field population, 0.01–0.25 solar masses.
Scientific discussion: `docs/research/independent/FIELD_MASS_FUNCTION_20261009.md`.
The lower cutoff and the smooth brown-dwarf power law are assumptions.

## Products

- `field_mass_function.png/pdf`: absolute density and normalized cluster comparison.
- `field_mass_function_checks.png/pdf`: component luminosities, apparent-single
  luminosities and temperatures, fitted separately.
- `field_mass_function.csv`: densities per cubic parsec per dex. `density_nominal`
  is the plotted clear-atmosphere line; the other columns summarize the equal-weight
  sensitivity ensemble. Stellar rows instead summarize the census inversion.
- `field_census.sqlite`: queryable source observations, empirical luminosities,
  individual stellar estimates, model histories and population fits.
- `selected_tracks.csv`: exact selected model inputs, history paths and hashes.
- `population_final_fits.json`: fit assumptions, intervals and predicted observables.
- `fit_*.npz`: joint slope/age weights and fitted normalizations.
- `component_luminosity_function.npz`: Monte Carlo counts, covariance and separate
  measured/imputed luminosity contributions.
- `validation.json`: population quadrature, metallicity-node and lower-cutoff checks.
- `*_receipt.json`: selection, uncertainty and database summaries.

The selected grid contains 69 distinct tracks and 6551 saved states. Of these,
44 reach the requested 10 Gyr, 24 reach the 500 K atmosphere boundary, and one
metal-rich clear 0.05-solar-mass run reaches 9.383 Gyr. This is not a grid of
69 complete 10 Gyr tracks. Old cold objects are censored, with a tested bound
for interpolation; the final incomplete warm interval uses adjacent completed
mass tracks. Both approximations are described in the report.

## Reproduction from saved inputs

Run from the repository root. These steps consume saved models; they do not
launch new evolution jobs. Environment variables apply only to these commands.

```sh
export PYTHONPATH=out/cluster-mass-functions-oct9-v1/deps
export MPLCONFIGDIR=/private/tmp/ember-cluster-mpl
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
field_python=out/science-env-sept18-v1/bin/python
field_work=out/field-mass-function-oct9-v1
$field_python $field_work/prepare_observations.py
$field_python $field_work/component_luminosity_function.py
$field_python $field_work/run_inference.py
$field_python $field_work/stellar_mass_function.py
$field_python $field_work/validate_population.py
$field_python $field_work/plot_results.py
$field_python $field_work/finalize_catalogue.py
```

`prepare_observations.py` uses immutable downloaded files in `raw/` and
`literature/`, plus the empirical bolometric calibrations and Sanghi catalogue
in the cluster workspace. `downloads.json` records URLs, retrieval times and
SHA-256 values. Source evolutionary mass columns are preserved for provenance
but never consumed in the inference.

Evolution inputs are in `tracks/*/configuration.txt`; launch/completion records
give commands, executable hashes and stop reasons. Source derives from frozen
commit `087351c` in `out/luhman16-independent-oct8-v1/source` and the cluster's
private initial-model continuation. The private post-deuterium executable
allows exactly zero initial deuterium; it does not alter the evolution equations.
Its build command and provenance are in `physics/post_deuterium_build.json`.
Shared `src/` and `scripts/` were not changed for this request.

Reused atmosphere data comprise SPHINX, Sonora Diamondback and Sonora Bobcat,
with a specified cool cloud-clearing join. Metallicities −0.3 and +0.3 use
tables in `out/field-metallicity-oct9-v1/physics/`. Their builder scripts and
source receipts are preserved there. Table construction reused existing
published atmospheres rather than computing new spectra.

## Interpretation and scope

The field fit uses a present-age density proportional to exp(beta*age/Gyr).
A negative beta means a younger local population. It need not mean a changing
Galactic birthrate. The luminosity likelihood retains the Monte Carlo covariance.
The separate single-object check uses published asymmetric errors. The young
L-dwarf constraint is transferred to the component fit as an explicit approximation.

The nominal model and uncertainty ensemble assume a 0.01-solar-mass lower
cutoff. A 0.014 control is saved separately and is not silently included in
the plotted curve below that cutoff. The equal-weight sensitivity cases have
no posterior model probabilities. The displayed band is not a complete
systematic-error interval. The low-mass power-law fit and stellar inversion
are separate; no physical discontinuity is inferred at their 0.08 join.

Atmosphere composition stays within declared fixed-table tolerances. Some
boundary controls use a shallower envelope and no fingering mixing; microscopic
transport and the 1% envelope thermal-reservoir limit remain enabled. Initial
contraction ages and cloud clearing remain model assumptions. These limitations
matter especially at young ages and near the hydrogen-burning boundary.

The all-sky census uses 20 pc for temperatures at least 600 K and 15 pc for
450–600 K; colder objects are excluded. Objects with missing component fluxes
are explicitly flagged. Known multiples are not counted as both systems and
children. Unknown binaries are bounded only by the declared equal-light-pair
scenario; this is not a full multiplicity selection model. Parent-catalogue
omissions and a metallicity–age relation cannot be reconstructed from these data.

## Example database queries

```sql
SELECT luminosity_method, COUNT(*)
FROM luminosity_observations
GROUP BY luminosity_method;

SELECT "case", mass_Msun, metallicity, cloud, end_age_Gyr,
       end_Teff_K, reached_requested_age
FROM ember_tracks
ORDER BY metallicity, cloud, mass_Msun;

SELECT "case", alpha_median, beta_median, chi2_distribution
FROM population_fits
WHERE observable = 'component_luminosity';
```

Raw machine-readable data retain precision. Human-readable discussion and
figures use at most four significant figures. No work was assigned to Opus;
his Luhman analysis and the paused WD work were left separate.
