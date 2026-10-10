# Nearby field mass function with Ember

The first corrected field comparison favors substantially more objects at
0.03–0.08 solar masses per 0.10–0.20 solar-mass star than either cluster
catalogue currently contains. This is a present-day local population estimate,
not a reconstruction of the clusters' birth populations.

[Mass functions and cluster comparison](../figures/field-mass-function.png) ·
[Fits to the measured distributions](../figures/field-checks.png)

| Population | N(0.03–0.08) / N(0.10–0.20) |
|---|---:|
| Field, nominal clear-atmosphere fit | 0.67; statistical interval 0.60–0.74 |
| Field, across declared sensitivity cases | 0.63; central range 0.55–0.71 |
| Pleiades, observed associated candidates | 0.34 |
| Hyades, observed associated candidates | 0.13 |

The cluster selections include identified escaped members and tails. Their
missing members, contamination and unresolved binaries remain uncorrected.
The apparent deficits therefore cannot yet establish different efficiencies
of brown-dwarf formation. The mass interval is fixed for comparison; its upper
edge is not an assertion of a universal hydrogen-burning limit.

## Observations and corrections

The principal sample is the published
[Kirkpatrick et al. 2024 census](https://doi.org/10.3847/1538-4365/ad24e2),
downloaded from the [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/docs/20pcCensus.html).
Removing combined-system rows leaves 3582 component entries, of which 3496
pass a hard 20 pc distance cut before further selection. White dwarfs are
excluded. Published evolutionary masses are never inference inputs. This is
the archived 2024 release, not an assertion that every discovery through 2026
has been incorporated.

The calculation uses the published completeness volumes: 20 pc above 600 K
and 15 pc at 450–600 K, excluding colder objects, with the census's Galactic
plane corrections. Parallax uncertainty, a local-volume distance prior,
bolometric errors, missing photometry and counting uncertainty are propagated.
Missing-photometry corrections distinguish apparent singles from known
multiples. Known components are counted once; parent-system light is not
assigned to both children.

Luminosities preferentially use integrated SEDs and published bolometric
estimates, then empirical photometric corrections. Missing component
luminosities contribute about 20% of the fitted counts, reaching 40% in one
bin. These are explicitly flagged estimates conditional on catalogue
temperature, with at least 0.18 dex scatter and a shared calibration error;
they are not direct flux measurements. This is a material limitation of the
component result. The raw catalogue, source matches and method flags remain
available in the database.

As a separate check, Ember fits the completeness-corrected 25 pc
single-object luminosity function of
[Best et al. 2024](https://doi.org/10.3847/1538-4357/ad39ef).
That sample excludes known close binaries and wide companions. The two
catalogues overlap and are **never multiplied as independent likelihoods**.
Both luminosity fits use 14 ≤ Mbol < 20. The single-object fit gives a
consistent shape, with lower normalization as expected from its selection.

## Ember inference

Brown dwarfs are fitted as an evolving population: dN/dM ∝ M⁻ᵅ and the local
present-age distribution p(t) ∝ exp(βt/Gyr), over ages up to 10 Gyr. A single
age is not assigned to the field. The nominal clear-atmosphere fit gives
α ≈ 0.3 and β ≈ −0.37. Negative β means fewer old objects in the local volume;
it does not uniquely measure the Galactic birthrate. The smooth substellar
curve is an assumed power-law population fit, not a histogram of individually
measured brown-dwarf masses.

The young L-dwarf fraction from Best et al. constrains the age distribution.
Transferring that single-object age constraint to the component population is
an approximation. Varying its age threshold from 100 to 300 Myr moves the
clear-atmosphere median α from about 0.7 to about 0.0, the largest tested
shape sensitivity. The adopted lower population cutoff is 0.01 solar masses;
raising it to 0.014 shifts the nominal median α to 0.4. The cutoff itself is
not measured by this analysis, and the census is not complete to that mass
at every age.

Both clear and cloudy Ember tracks are used. Metallicity is represented by
[Fe/H] = −0.3, 0 and +0.3 with weights 0.235, 0.63 and 0.135, fixed Y = 0.27
and scaled metal abundances. Their mean is −0.03 dex, consistent with the
nearby M-dwarf result of [Dittmann et al. 2016](https://doi.org/10.3847/0004-637X/818/2/153).
The adopted spread is a scenario, not an object-by-object measurement; an
age–metallicity correlation is omitted. The cloudy prescription clears toward
the existing cloud-free boundary below 1200 K.

Above 0.08 solar masses, individual luminosities are mapped through the Ember
age/metallicity grid. The 0.10–0.20 reference density is about 0.0151 objects
per cubic parsec; an inner 15 pc check gives 0.0153. The stellar inversion and
brown-dwarf population fit are separate, so their small mismatch near 0.08
solar masses must not be interpreted as a physical discontinuity.

Residual unknown binaries are explored by treating 0–5% of apparent singles
as equal-light pairs. This brackets one plausible bias; it is not a complete
mass-ratio and separation selection model. The shaded field band combines
statistical draws with age, cloud, metallicity, initial-age, cold-boundary and
multiplicity scenarios at equal weight. It is **not** a calibrated posterior
covering every systematic uncertainty. Cluster shading has the more limited
meaning documented in the cluster report.

## Checks and remaining limits

The clear-atmosphere component luminosity fit has χ² = 2.57 over six bins;
the separate single-object fit has χ² = 0.85. Temperature-function fits are
less successful: χ² = 14.8 for clear atmospheres and 48.9 for cloudy ones over
ten bins. These diagnostics remain visible. A good luminosity fit alone does
not validate the cloud or temperature calibration.

Doubling the population integration's mass and age sampling changes predicted
luminosity counts by less than 0.1%. Interpolation reproduces calculated
metallicity nodes. Envelope-depth and initial-composition controls change
overlapping model luminosities by at most 1.9% in the comparisons performed.
These are useful checks, not a full evolutionary convergence study.

Some low-mass tracks stop at the real 500 K atmosphere edge, already fainter
than the fitted luminosity range. Their later states are censored; alternative
cold bounds test interpolation between masses. One metal-rich 0.05 solar-mass
track ends at 9.383 Gyr. Later populations interpolate neighboring completed
masses; a conservative upper bound on the affected density is 3.88 × 10⁻⁵
objects per cubic parsec at the 97.5th percentile, below the observational bin
errors. Neither case is presented as a computed continuation.

Several tracks use compact post-deuterium starts with an approximate omitted
contraction age. Controls vary that correction from zero to twice its estimate.
Private shallow-envelope runs retain the 1% envelope thermal-reservoir limit;
some disable fingering mixing and use a fixed-composition atmosphere within
declared composition tolerances. Thus this is a corrected, reproducible first
field result, with specific remaining population and model approximations.

The next improvement with the most direct value is better component
luminosities and a joint multiplicity/selection model, followed by a common
completeness treatment for the cluster comparison.

## Files

Workspace: `out/field-mass-function-oct9-v1/`.

- `README.md`: processing order and reproduction commands.
- `field_census.sqlite`: source catalogue, luminosities, Ember tracks and fits.
- `field_mass_function.csv`: nominal curve and sensitivity summaries.
- `selected_tracks.csv`: exact inputs and histories selected for inference.
- `field_mass_function_receipt.json`, `validation.json`: numerical results.
- `component_luminosity_receipt.json`: corrections and missing-flux contributions.

Underlying luminosity data also use [Sanghi et al. 2023](https://doi.org/10.3847/1538-4357/acff66)
and its [published catalogue](https://doi.org/10.5281/zenodo.10086810).

See the [research index](../README.md) for the combined interpretation and
[reproduction guide](../REPRODUCING.md) for published and locally retained data.
