# Hyades and Pleiades: provisional observed mass distributions

The first observed candidate mass distributions and individual Ember mass
estimates are available in [the figure](../figures/cluster-mass-functions.png).
They include identified escapees, tidal tails and associated populations.
The plotted decline toward low masses cannot yet be interpreted as the birth
mass distribution: faint-object coverage is much narrower than the stellar
and tail surveys, and contamination and unresolved binaries are uncorrected.

## Observations

The merged list contains 1993 Hyades and 5703 Pleiades/associated candidates
selected by at least one source, across all masses. These are neither confirmed
member counts nor counts below 0.25 solar masses. Exact Gaia identifiers and
epoch-corrected coordinates remove duplicate entries. Original probabilities,
selection flags, photometry and catalogue provenance remain available separately.

The sources include the [2026 Pleiades census](https://arxiv.org/abs/2607.17300),
[Liu 2025 memberships](https://arxiv.org/abs/2504.08179),
[Risbud 2025 tails](https://arxiv.org/abs/2501.17225),
[Jadhav 2025 catalogue assessments](https://arxiv.org/abs/2508.15056),
[the Greater Pleiades Complex](https://arxiv.org/abs/2511.07533),
DANCe, and deeper substellar candidate lists. Catalogue grades describe an
entire recovered tail; they are not individual membership probabilities.
The selection union likewise is not a combined membership posterior.

Gaia photometry and errors were retrieved for 7066 distinct selected sources.
Published SED luminosities are preferred. An empirical Gaia bolometric
correction, using [Cifuentes 2020](https://doi.org/10.1051/0004-6361/202038295)
luminosities rather than published masses, supplies a stellar control.
Against 1176 Pleiades SED estimates, its median difference is 0.01132 dex
and robust scatter 0.01844 dex. Extinction, binaries and substellar spectra
still require explicit treatment.

Luminosities are now available for 1480 selected Hyades and 4295 selected
Pleiades sources. Most sources without luminosities are brighter objects
outside the requested mass range. Seventeen Hyades sources have direct matches
to the observational SED compilation of
[Sanghi et al.](https://arxiv.org/abs/2309.03082); its photometry and parallaxes
also calibrate an empirical faint-star K-band correction. Its published
evolutionary masses, radii and temperatures are not inputs to the Ember masses.
For faint candidates without measured distances, a common Hyades space velocity
provides membership-conditional kinematic distances with a working 20% error.
Measured distances take precedence. All eight candidates published by
[Schneider et al. 2024](https://arxiv.org/abs/2408.10112) are included, preserving
their two different membership-probability estimates.

## Ember calculations

Fifty-seven cluster and control calculations have reached their requested ages;
this batch has finished. Successful tracks span 0.008 to 0.25 solar masses,
with separate clear, cloudy and initial-state controls. These are not 57
independent masses. The current inventory records each case.

The atmosphere inputs combine published SPHINX and Sonora profiles with the
native Gibbs EOS for deep adiabatic extensions. Their common warm range agrees
within 3.034% in temperature at 100 bar. Existing
[Sonora Bobcat profiles](https://doi.org/10.5281/zenodo.5063476) now extend the
tables to 500 K, using their published densities for cold geometric depths.
The cloudy scenario joins the cloud-free cold profiles over 900–1200 K;
this is a cloud-clearing assumption, not a new cloudy atmosphere calculation.
The blended tables cover 500–4000 K; separate cool tables extend to lower
gravity. Interior composition assumes
Y=0.27 and scaled GS98 metals: [Fe/H]=0 gives Z=0.01648, while +0.15 gives
Z=0.02307. The published atmosphere mixtures are not exact element-by-element
matches to that interior mixture.

Most tracks start after deuterium exhaustion. Their elapsed times cannot be
used directly as cluster ages: a polytropic contraction estimate gives missing
early-time scales of roughly 10–42 Myr. The inference subtracts this estimate
from the assumed cluster age and varies its scale from zero to twice the
estimate. The primordial-deuterium cases still need initial-entropy checks.

Earlier initial-state controls at 0.15 solar masses show that the two starts
follow essentially the same contracting sequence with different clock origins.
Shifts of 12.75 Myr (Pleiades composition) and 12.22 Myr (Hyades composition)
align the luminosities over 30–150 Myr, with RMS residuals in ln L of 0.0003429
and 0.00002587. Without this shift, luminosity differs by 5.386% at 125 Myr
and 0.004853% at 650 Myr. This measures relative starting-age sensitivity;
it does not independently establish the age since formation.

At 0.25 solar masses, applying the analytic starting-age correction brings the
earlier and compact starts into luminosity agreement within 0.6565% at 125 Myr
for the Pleiades composition, and 0.1940% at 650 Myr for the Hyades composition.
This supports the approximate clock correction in this mass range, without
establishing deuterium-burning ages at the faint end.

The cold extension lets primordial-deuterium Pleiades tracks at 0.008, 0.010,
0.012 and 0.014 solar masses reach 200 Myr. Warm low-gravity coverage still
blocks equivalent 0.016 and 0.020 starting models, so post-deuterium tracks
remain in that part of the provisional grid. The relation between luminosity
and mass folds near this join. All mass solutions are retained, weighted by
the inverse local slope for a conditional flat prior in log mass. The feature
could reflect both deuterium burning and inconsistent early conditions; it
is not yet a validated prediction.

## Interpretation and remaining work

The [normalized overlay](../figures/normalized-cluster-mass-functions.png) divides
each distribution by its own candidate count between 0.10 and 0.20 solar masses.
This compares shapes without using the unequal faint survey depths to set the
normalization or relying on the upper model boundary at 0.25 solar masses.
The same normalization is applied separately to every sensitivity realization;
its integral over the reference interval is checked to equal one. Median
reference counts are 361 for the Hyades and 590 for the Pleiades. This does not
establish completeness even in that interval. Associated members and tails
remain included. The normalized shading excludes counting uncertainty, as well
as the unmodeled selection and binary effects.

The figure uses 256 sensitivity scenarios: ages of 600–800 Myr for the Hyades
and 110–140 Myr for the Pleiades, clear/cloud-clearing alternatives, starting-age
uncertainty, observational luminosity errors and a common 0.04 dex model shift.
The shaded ranges describe these choices, not a calibrated membership posterior
or a complete error budget. Masses outside actual track coverage are not
clamped to the grid edges. Unresolved sources are treated as single objects.

The merged data do not define one detection limit across the sky. In particular,
the [deep Pleiades survey](https://arxiv.org/abs/1407.2849) covers only 0.8 square
degrees, whereas the associated population extends far beyond the core.
Including every identified tail member is essential but does not recover
undetected low-mass escapees. Next work must resolve the deuterium/start-state
join, finish the deep-survey selection audit, and quantify completeness,
contamination and binary effects before estimating a birth mass function.

## Reproducible files

Workspace: `out/cluster-mass-functions-oct9-v1/`.

- `cluster_catalogues.sqlite`: objects, original source evidence and luminosities.
- `catalogues/{hyades,pleiades}_members_and_ember_masses.csv`: complete selected
  lists with provisional mass estimates and explicit missing/out-of-range values.
- `observed_mass_functions.png/pdf`, `observed_mass_function_bins.csv`: figure
  and numerical distributions.
- `observed_mass_function_receipt.json`: inference choices and model provenance.
- `mass_inference_checks.json`: uniqueness, database integrity, inversion and
  multiple-solution checks; all pass.
- `catalogues/`: CSV exports, matching records and extraction checks.
- `downloads.json`, `raw/`, `literature/`: source records and hashes.
- `physics/`: atmosphere builders, tables and comparison receipts.
- `tracks/`: configurations, histories and completion/timing receipts.
- `track_inventory.json`: updated results, active cases and initial-age scales.
- `initial_age_control.json`: earlier-start comparisons and relative clock fits.
- `initial_age_upper_mass_check.json`: corrected 0.25-solar-mass comparison.
- `catalogue_coverage.png`: observed spatial and luminosity distributions.
- `README.md`: processing order and outstanding checks.

The standard executable is from commit 087351c. A private initialization-only
variant uses static entropy continuation to obtain earlier seeds; its source,
patch and disabled-control check are recorded under physics/. Physical evolution
is unchanged. Shared solver code and paused WD development were not changed.

See the [research index](../README.md) for the combined interpretation and
[reproduction guide](../REPRODUCING.md) for published and locally retained data.
