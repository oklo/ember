# Conductive layers at the hydrogen-burning boundary

Test whether a nonconvective layer separates an inner convective region from
the envelope near the burning boundary, and whether this association persists
as initial metallicity changes. Establish the layer's existence separately
from any claim that it eventually causes burning to fail.

## Present evidence

At initial Z = 0.08241 and Y = 0.27, the clear 0.06125-solar-mass control
retains a separating layer through 359.8 Gyr with the repaired, extended
metal EOS. Its main layer spans approximately m/M = 0.3577–0.3830; a second,
one-cell interval lies at 0.3984–0.4035. The main layer has five nonconvective
mesh intervals. Conduction supplies 95.73–97.26% of its radiative-plus-conductive
transport. Composition gradients also stabilize the layer.

A repeat from 173.0 to 239.9 Gyr using the extended EOS retains the same
convective boundaries. The nearby 0.0614-solar-mass control remains fully
convective at 1 Tyr. Reducing conductivity by 10% in the 0.06125 control
also prevents separation through 1 Tyr. That perturbation is a sensitivity
test, not an established conductivity uncertainty.

These calculations still use the atmosphere response assessed with an
incorrect native opacity normalization. Corrected pairs at 913, 1184 and
1657 K now pass their column checks, but do not yet constitute a corrected
family covering formation and subsequent evolution. The present evidence
establishes persistence under the EOS change,
not formation under every corrected input, robustness across metallicities,
or eventual burning failure. The older cloudy 0.05835-solar-mass example
omitted metal electrons from its EOS and cannot serve as a matched metallicity
comparison.

Fresh runs with 536, 1072 and 2144 points all form a separating layer near
31 Gyr. Its later width is sensitive to the timestep at formation: at 60 Gyr,
the 1072-point run with a 1 Gyr step ceiling is 11.32% fainter than the
536-point run, whereas the two meshes agree within 0.9170% in luminosity
with a 200 Myr ceiling. The 2144-point, 1 Gyr calculation recovers the broader
layer. At 200 Myr steps, the 2144-point model retains 34 nonconvective
intervals and differs from the 1072-point model by 1.953% in luminosity.
Turning off Richardson extrapolation does not remove the layer; with short
steps it changes luminosity by 1.382% and one band boundary by one cell.
These controls retain the old atmosphere response.

At first separation the short-step model's layer is thermally stable even
without a composition gradient. By 60 Gyr its composition gradient is needed
to prevent convection. Slow mixing in that later state therefore needs its
own physical sensitivity test; persistence under instantaneous convection
alone is insufficient. The [compact numerical record](conductive_band_checks.json)
contains endpoints, band coordinates, comparisons and input hashes.
The full local diagnostics are in
`out/highz-native-opacity-oct12-v1/band_night_oct12/assessment.json` and
`face_stability_assessment.json` in the same directory.

A mixing-only Langer sensitivity with efficiency 0.1 retains separation
through 150 Gyr, but ends at 1381 K rather than the control's 1713 K, with
65.91% less luminosity. Both remain nuclear supported at that age. Each has
two distinct nonconvective intervals separated by a small convective region.
An otherwise identical efficiency-0.01 calculation also retains two layers
and reaches 1559 K, between those controls, at 150 Gyr.
The experiment uses total thermal conductivity and composition-dependent
buoyancy; material transport carries enthalpy, but turbulent heat transport
is omitted. It tests sensitivity rather than validating a unique mixing
rate in degenerate matter. The result is recorded in
`slow_mixing_review_oct13.json` alongside the numerical diagnostics.

Atmosphere inputs also matter. Adding a third metallicity source changes
the 60 Gyr luminosity by 1.458% and moves first separation from 30.82 to
26.24 Gyr through the changed composition interpolation. Separately,
multiplying boundary pressure by 1.095 throughout a fresh calculation,
using the measured 1184 K offset, still produces a band but yields a maximum
nuclear fraction of 0.9081 followed by cooling. This pressure-only experiment
is not a corrected atmosphere family: the measured offset varies along the
track. Neither test determines the corrected burning boundary. That boundary
must be bracketed again before comparing the fate of borderline masses.

The first matched three-times-reference controls use the same executable,
EOS, interior opacity method and expanded atmosphere family as the five-times
controls. At 0.064 solar masses, a separating band appears at 8.978 Gyr,
spanning m/M = 0.3040–0.3728 with 14 mesh intervals. The object subsequently
cools; by 16.17 Gyr the conductive region includes the center and its surface
reaches the 1100 K table edge. At 0.0665 solar masses the model remains fully
convective and nuclear supported through 150 Gyr. At 0.06588 solar masses,
a band forms at 31.91 Gyr and 1222 K, spanning m/M = 0.3627–0.4087 with
nine intervals and about 99.75% conductive transport. Nuclear luminosity
reaches 110.2% of surface luminosity near 99.51 Gyr. At 150 Gyr the star is
at 1787 K, nearly back to energy balance, and retains two separating bands.
An EOS reconstruction finds eight of the nine first band intervals thermally
stable without composition gradients; by the nuclear-fraction maximum, all
nine require composition stabilization. The reconstruction retains the birth
metal ratios because the profile output records total Z.
It has consumed only 0.3337% of its initial hydrogen, while central X has
fallen by 12.39% and central Z has risen from 0.04945 to 0.07573. Composition
redistribution therefore matters substantially. These are results at a second
metallicity with matched interior physics; the old H/He atmosphere response
remains, and the corrected long-term burning boundary is not yet located.

Additional native atmosphere checks qualify at 1100 and 1300 K. Compared
with the expanded three-source atmosphere family, the 5Z matching pressure
is higher by 7.703% at 1100 K, 7.764% at 1184 K, and 11.32% at 1300 K.
The 3Z correction at 1300 K is 10.52%. Their temperature dependence rules
out treating a uniform pressure factor as a corrected atmosphere grid.
A sparse interpolation test with independent temperature, gravity and
composition checks is in progress before stellar incorporation.

## Matched sequences

1. Start with Z/Z_ref = 1/3, 1, 3 and 5, where Z_ref = 0.01648. Hold
   Y = 0.27 and the metal abundance ratios fixed. Treat clear and cloudy
   atmospheres separately, with one specified cloud prescription. Extend
   to lower and higher Z only where the actual evolving compositions have
   supported EOS, radiative opacity and atmosphere inputs. Later vary helium
   independently for enriched birth mixtures.
2. Use one executable and physical treatment across the comparison. Reuse
   completed tracks only when their inputs match. In particular, use corrected
   H/He atmosphere responses before interpreting new precise mass boundaries.
   Check coverage along the evolving surface and enriched core, not only at
   the birth composition.
3. For each composition, bracket cooling and sustained burning, then add
   intermediate masses adaptively. Compare masses relative to that
   composition's own transition. Do not reuse one absolute mass as the
   borderline case at every Z. Include at least one clearly cooling and one
   clearly burning control beyond the closely spaced transition models.
4. Continue through 1 Tyr even after nuclear luminosity reaches or exceeds
   the surface luminosity. A falling nuclear fraction below 0.1% permits an
   early stop. Record early stops explicitly; nuclear balance at one age does
   not establish permanent main-sequence membership.

## Record the layer and its consequences

Save convective partitions and the intervening nonconvective intervals,
their mass and radial widths, number of mesh intervals, onset age and duration.
Retain local conduction fractions, temperature, density, composition and
degeneracy. Distinguish the required radiative-plus-conductive gradient from
the adiabatic gradient and the composition-dependent stability threshold.
High conductivity alone does not establish a separating band.

Measure hydrogen and metal transfer across both boundaries, the mass and
fuel inventory of the inner convective region, and nuclear luminosity versus
surface luminosity. Electron conduction carries heat across the layer;
microscopic diffusion can still transport ions. The layer is not an
impermeable chemical boundary or a thermal insulator.

## Focused robustness checks

- For the decisive layer-bearing mass at each Z, repeat with approximately
  twice the spatial resolution and shorter timesteps. Resolve the narrowest
  layer; do not treat a persistent one-cell interval as spatially converged.
  Start before separation so remapping cannot create the proposed feature.
- Repeat selected cases with modest conductivity perturbations, and with
  justified finite mixing or additional slow mixing where the local stability
  criteria permit it. Diagnose formation and later composition stabilization
  separately. Do not assume a missing mixing process has an arbitrary rate.
- Compare corrected atmosphere responses and their coverage/interpolation
  against the old-response controls. A continuation after a boundary change
  is not a substitute for a matched run that develops its own composition
  gradients from an earlier state.

Report, for each Z and atmosphere choice, the mass interval showing a
resolved, persistent band alongside the independently measured burning
boundary. Show numerical sampling bounds separately from physical
sensitivity. A band that disappears or moves away from the burning boundary
is a result to report, not a failed test to hide.

The first priority is the corrected-atmosphere comparison for the current
five-times-reference controls. Broadening the grid before that comparison
would multiply runs with a known input limitation.
