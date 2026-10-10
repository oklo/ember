# Hyades and Pleiades mass functions

The requested range is 0.25 solar masses to each survey's detection limit.
The scientific objective is the birth mass distribution. Include identified
escaped members, tidal tails and associated populations. Keep probabilities
and association classifications rather than equating every candidate with a
confirmed member. Do not impose a bound-cluster cut on the main result.

Current progress: the catalogue union, source evidence, Gaia photometry and
provisional luminosities are in cluster_catalogues.sqlite. There are 1993
Hyades and 5703 Pleiades/associated candidates selected by at least one source,
across all masses. These are candidate counts, not confirmed membership or
the requested mass-function counts. No completeness-corrected mass function
or birth mass distribution is claimed yet.

Downloaded sources and hashes are in downloads.json and raw/. Literature is
in literature/. Run scripts using out/science-env-sept18-v1/bin/python with
this directory's deps/ on PYTHONPATH. Processing order:
normalize_catalogues.py, build_membership.py, fetch_gaia.py,
derive_luminosities.py, plot_catalogue_checks.py,
extend_substellar_luminosities.py, infer_observed_mass_functions.py,
export_cluster_results.py. The Gaia fetch reuses previous
results and only queries newly required IDs.

Database tables: objects, source_records (original rows, selection, source
probabilities and catalogue grades), objects_with_luminosity,
objects_for_mass_inference and ember_candidate_masses.
CSV views are in catalogues/. catalogue_coverage.png shows observed spatial
and luminosity distributions; it is a diagnostic, not the final mass function.

Important data details:

- The DANCe TAP join contains storage-unit columns mislabeled with physical
  units. The script corrects them and compares 23 shared objects with an
  independent ASU export. Use catalogues/pleiades_dance_resolved.csv.
  raw/pleiades_dance_join.csv is uncorrected; do not use it for inference.
- Rounded DANCe identifiers are not unique. There are 5444 membership rows,
  5441 unique position-based IDs and 5457 joined photometry rows. Ten IDs have
  multiple join rows. Original numeric Sarro identifiers now resolve all
  5444 rows one-to-one; counting the uncorrected join would overcount them.
- The Bouy ASU exports are truncated and serve only as unit checks. The TAP
  join contains the requested probability-selected sample.
- Source Gaia identifiers are preserved as strings. Published mass estimates
  are retained as original columns, not accepted as Ember masses.
- Exact Gaia IDs and propagated positions combine sources without adding
  their counts. Source selections and probabilities remain separate. A union
  of selections is not a calibrated combined membership probability.
- The Jadhav 2025 compilation preserves the original catalogue and its
  gold/silver/bronze assessment. These grades describe whole tail catalogues,
  not individual membership probabilities. Duplicate publications are not
  independent evidence. DR2 identifiers are not relabelled as DR3 identifiers.
- All Risbud 2025 candidates are retained. Its corrected-sample flag applies
  additional magnitude and symmetric tail-length cuts for comparisons; those
  cuts would exclude associated populations requested here.
- Deep central surveys and wide stellar/tail surveys have different selection
  functions. Membership-classifier recovery is not detection completeness.
  A missing brown-dwarf tail census cannot be corrected by a single sky-area
  factor or used as evidence for a deficit in the birth population.

Luminosities: use the published DANCe SED estimates where available, rescaled
to individual distance. An empirical Cifuentes 2020 Gaia bolometric correction
provides a stellar control; no published evolutionary masses enter it. Against
1176 Pleiades SED estimates its median difference is 0.01132 dex and robust
scatter 0.01844 dex. This does not validate faint-dwarf conversions, extinction,
unresolved binaries, or the full metallicity/gravity uncertainty budget.

Physics pilots: physics/ contains 100-bar atmosphere boundaries assembled from
published SPHINX I and Sonora Diamondback profiles, at [Fe/H]=0 and +0.15.
Interior abundances use GS98 Z/X=0.0231 and assumed Y=0.27, giving Z=0.01648
and 0.02307. Atmospheric helium and detailed metal mixtures remain those of
the published grids; this is not an exact element-by-element mixture match.
Clear and f_sed=2 cloudy cases remain separate. The 2000--2400 K join agrees
within 3.034% in temperature at 100 bar. Native-EOS adiabatic extensions are
explicitly identified and checked against withheld source depths.

Current outputs: observed_mass_functions.png/pdf; observed_mass_function_bins.csv;
observed_mass_function_receipt.json; cluster_isochrones_provisional.csv;
catalogues/{hyades,pleiades}_members_and_ember_masses.csv; mass_inference_checks.json.
The candidate exports retain sources without a luminosity or model-supported mass.
The figure is also copied to docs/research/independent/CLUSTER_MASS_FUNCTIONS_20261009.png.

57 cluster/control runs reached their targets; none remain active in this batch.
track_inventory.json supersedes pilot_track_summary.json. Cold tables now extend
to 500 K using existing Sonora Bobcat profiles; see physics/colder_boundary_receipt.json.
The cloudy case clears below the Diamondback range as an explicit scenario.
Published profile densities provide the cold geometric depth, with a native
Gibbs adiabat extending deep profiles to 100 bar where needed.

Primordial-D Pleiades tracks at 0.008, 0.010, 0.012 and 0.014 solar masses reach
200 Myr. Warm low-gravity coverage still blocks primordial-D 0.016/0.020 starts.
Post-D starts remain in that part of the provisional grid. Their join can give
multiple mass solutions; the inference preserves all branches, with conditional
inverse-slope weights for a flat log-mass prior. This is not a validated feature
of the cluster mass function.

Most tracks start after deuterium depletion. The mass inference explicitly
subtracts an estimated missing contraction time GM^2/(7RL), varying its scale
from 0 to 2. Earlier/compact 0.15-solar-mass models align after relative shifts
of 12.75 Myr (Pleiades) and 12.22 Myr (Hyades). The age-corrected 0.25 controls
agree in luminosity within 0.6565% and 0.1940% at their respective nominal ages.
See initial_age_control.json and initial_age_upper_mass_check.json. This does
not validate deuterium timing or determine the absolute formation age.

Sanghi 2023 observational SED luminosities and empirical K-band corrections
extend the faint luminosities. Seventeen Hyades sources have direct SED matches.
Kinematic distances used for selected faint candidates are conditional on
Hyades membership; measured distances take precedence. Source hashes and all
calibration choices are retained in substellar_luminosity_receipt.json.
Published evolutionary masses, radii and temperatures are not inference inputs.

The figure shows observed candidate counts, including tails, with sensitivity
to age, clouds, initial conditions and measured luminosity. It has no survey
completeness, contamination or binary correction. It cannot yet establish the
birth mass function. Next: resolve the faint initial-state join, finish the
deep-survey checks, then quantify selection and unresolved binaries. The
0.8-square-degree deep Pleiades survey is not a census of the whole association.

The clean executable is from commit 087351c. physics/prepare_seed_continuation.py
provides an optional private initialization-only entropy ramp; its disabled
control exactly matches the original recorded seed. No shared solver code or
paused WD development was changed.

Greg directs Opus on Luhman 16AB. The primary owns this cluster work, including
tables and tracks. The earlier Opus cluster request is withdrawn.
