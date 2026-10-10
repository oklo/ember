# Cluster structure and survey follow-up

This continues the `cluster-dissolution-oct9-v1` experiment. The same privately
built PeTar executable supplies unsoftened collisional dynamics with SSE/BSE and
MWPotential2014. Ember supplies the low-mass light. No shared stellar source or
table is changed; no distant escaper is deleted.

`analysis_plan.json` fixes the broad mass/size checks before the new endpoint
results. They are screening ranges, not statistical confidence intervals.
`ensemble_plan.json` specifies six new runs: three component-mass draws each
for Hyades and Pleiades, with the existing binary assumptions. A second,
separate three-run Pleiades test uses the lower independently obtained field
normalization, R = 0.55; `field_sensitivity_plan.json` records that decision.
It is not a best-fitting IMF or a Hyades test at the same normalization.

Do not relaunch an ensemble. Inspect `runs/*/pid`, `receipt.json`, and any
`continuation/` before starting anything. Hyades seed 43 was deliberately
stopped after a saved 200-Myr state and continued with its tree timestep
doubled. Its original receipt has exit -15; the continuation is the final
result. `step_validation_s43` compares this choice against the original step.
No change was made to the binary grouping radius in the adopted continuation.
The alternative that enlarged that radius was slower and was not adopted.

Reproduction, from the Ember root, after all runs have finished:

```sh
export PYTHONPATH=out/cluster-dissolution-oct9-v1/deps:out/cluster-mass-functions-oct9-v1/deps
export MPLCONFIGDIR=/private/tmp/ember-cluster-mpl
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
out/science-env-sept18-v1/bin/python out/cluster-dissolution-oct9-v2/analyze_followup.py
out/science-env-sept18-v1/bin/python out/cluster-dissolution-oct9-v2/compare_populations.py
out/science-env-sept18-v1/bin/python out/cluster-dissolution-oct9-v2/audit_magnitude_coverage.py
out/science-env-sept18-v1/bin/python out/cluster-dissolution-oct9-v2/dance_multiplicity.py
out/science-env-sept18-v1/bin/python out/cluster-dissolution-oct9-v2/plot_followup.py
out/science-env-sept18-v1/bin/python out/cluster-dissolution-oct9-v2/plot_mass_shapes.py
out/science-env-sept18-v1/bin/python out/cluster-dissolution-oct9-v2/validate_and_archive.py
out/science-env-sept18-v1/bin/python out/cluster-dissolution-oct9-v2/write_report.py
```

The Gaia–2MASS fetch is cached with ADQL queries and checksums. It uses exact
Gaia identifier strings from the original SQLite object table, not the numeric
identifier column in the older spatial CSV. It fills observed J/K magnitudes
for the Hyades aperture comparison; it does not change bolometric luminosities
or inferred masses. Missing band measurements remain explicitly counted.
The filled Hyades photometry also improves the empirical band-luminosity
calibration. The direct DANCe magnitude check includes candidates whose
luminosities remain unavailable; these objects are not assigned invented masses.

`cluster_followup.sqlite` contains structural endpoints, observed and simulated
counts, conditional recovery diagnostics, initial components, final descendants
and execution cost. CSV draws retain the mass, radial and angular distributions.
The all-associated-member union remains the target population. Restricted
apertures diagnose coverage; they do not redefine cluster membership.

The DANCe check uses selected catalogue objects, radius below 3 degrees and
Ks below 17. The published aggregate 90% classifier recovery and 4.3%
contamination are an absolute-count sensitivity only. They do not measure
mass-dependent recovery or correct unobserved tails. The multiplicity exercise
splits the observed light using published near-equal-binary probabilities and
two assumed flux ratios. It is not a complete binary census or a recovered IMF.

These are a few initial-condition realizations, not a fitted population
posterior. Gas expulsion, primordial mass segregation, exact union selection,
field contamination and the larger Pleiades birth complex remain unmodeled.
