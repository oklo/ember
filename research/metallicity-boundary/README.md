# Hydrogen burning near the stellar boundary

Ember compares clear and cloudy atmospheres from one third to three times
Z_ref = 0.01648, with Y = 0.27 and solar metal ratios. The cloudy models use
f_sed = 2 and clear between 900 and 1200 K. Tracks start on Hayashi contraction
with deuterium.

Bars span tested cooling and burning masses with metal ionization and
matching physical inputs:

| Initial Z | Clear (solar masses) | Cloudy (solar masses) |
|---|---|---|
| 0.005494 | 0.075–0.07525 | 0.07025–0.0705 |
| 0.009516 | 0.0735–0.074 | 0.068–0.069 |
| 0.01648 | 0.07125–0.0725 | 0.0665–0.06675 |
| 0.02855 | 0.069–0.070 | 0.0635–0.0645 |
| 0.04945 | 0.066–0.0665 | 0.060–0.061 |

Squares have 99–101% nuclear support at 1 Tyr, varying by less than two
percentage points over the final 50 Gyr. Open circles have declined below
0.1%; diamonds are unfinished trials. These intervals measure mass sampling,
not physical uncertainty or permanent main-sequence membership. Shading
connects the five sampled metallicities. The
neutral-metal EOS comparison is separate;
different EOS treatments must not be combined into a tighter interval.

The high-Z clear cooling endpoint extends the permitted atmospheric H/He
range slightly. Its composition response remains provisional. An independent
test found a normalization error in the local native gas-opacity calculation
used to assess H/He responses. Corrected individual atmospheres now pass
their checks; a family covering formation and evolution is being assembled.
Passing flux checks under the old opacity does not validate the corrected model.
The imported Bobcat, Diamondback, SPHINX and MARCS tables were not generated
by this local calculation path. New continuations do not tighten the plotted
mass intervals. Cloud physics, composition response, mixing and dense-matter
thermodynamics remain larger uncertainties than the mass spacing in several
intervals.

## Replot

Install the [research-record dependencies](../record/requirements.txt), then:

```sh
python3 research/metallicity-boundary/plot_transition_mass.py --include-comparison
python3 research/metallicity-boundary/plot.py
```

The scripts read the retained JSON/CSV samples and write local figures.
Generated figures are not versioned. `transition_mass.json` and `provenance.json`
retain source identities and assumptions. `numerical_checks.json` and
`conductive_band_checks.json` retain compact numerical comparisons; the
experimental high-metallicity and cloudy continuations have separate data
under their respective subdirectories. They do not tighten the intervals above.
