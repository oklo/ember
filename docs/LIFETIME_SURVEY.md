# Toward converged lifetimes as a function of mass and composition

Checkpoint: 2026-09-10. This is the scientific target and remaining work,
not a claim that Ember has produced a complete lifetime or a metallicity survey.

The current shared-mixture calculation follows a fixed baryonic 0.1 Msun
star at Z=.02 from its specified static main-sequence initial model. At
3.40 trillion years it remains fully convective, with XH=.2034,
R=.1459 Rsun, L=.002016 Lsun and Teff=3202 K. The active objective
is an evolved helium remnant through Teff=100 K; see
[COLD_REMNANT.md](COLD_REMNANT.md). Earlier numerical controls below refer
to the preserved input families and do not establish later endpoint convergence.
At the two-trillion-year endpoint, an independent repeat reproduces the full
history and final profile byte for byte; see
[the repeat check](results/evolution_metal_2tyr_gas_repeat.json).
Doubling the mesh at two trillion years changes luminosity by +.06163% and
radius by -.01234%; fourfold tighter time tolerances change L by -.003610%
and R by -.000594%. See
[the two-trillion-year comparison](results/evolution_metal_m010_2tyr_gas_convergence.json).
For the completed 1.3-trillion-year control,
doubling the mesh from 512 to 1024 points changes the final luminosity by
+.01917%; fourfold tighter time tolerances change it by +.0006203%. These
quantify sensitivity of that segment, not the uncertainty of a hydrogen-burning lifetime.
See [the comparisons](results/evolution_metal_1300gyr_convergence.json).

[Native restarts](RESTART.md) now preserve the internal state and timestep
controller. Short uninterrupted and restarted trajectories agree exactly,
including the final profiles; hashes bind the executable and inputs. The
checkpoint-enabled executable also repeats the entire 2.5-trillion-year
track byte for byte, and its saved state has continued through 2.75 to 2.85 trillion years.
These are reproducibility
checks, separate from physical validation or endpoint-age convergence.

The evolution driver still fixes the initial mass, XH=.7 and Z=.02.
Its age excludes formation and pre-main-sequence contraction. The detailed
GS98 EOS and non-grey atmosphere families also fix Z=.02. The helium-rich
extension lowers hydrogen at fixed Z; it is not a high-metallicity grid.
Interior opacity families at Z=.01/.02/.03 cover only a modest interval.

The present work establishes explicit non-equilibrium H/He3 burning,
thermodynamically consistent EOS interpolation, molecular non-grey boundaries,
electron conduction, and conservative burning/mixing with Ledoux transport.
Condensate depletion is being tested as a zero-grain-opacity limiting case.
The atmosphere grids, source data and final checks are described in
[FORWARD_EVOLUTION.md](FORWARD_EVOLUTION.md). Mixture differences between
the source EOS and opacities, and the depletion/gas-chemistry hybrid, remain
physical approximations despite tight numerical checks.

## What constitutes the lifetime result

1. Finish a reference track through radiative-core formation and hydrogen
   exhaustion. Extend temperature, density and composition support along the
   actual path, including low XH. Resolve the moving convective boundary and
   test mesh adaptation if a fixed mesh cannot do so. Check the significance
   of microscopic diffusion, additional nuclear branches and thermal neutrino
   losses in the newly encountered states before choosing approximations.
2. Define the endpoints before comparing models: central hydrogen exhaustion
   at a stated threshold and the end of sustained hydrogen burning should be
   recorded separately. A remaining envelope or shell can retain fuel after
   core exhaustion. Vary the exhaustion threshold to quantify its influence.
   Near the hydrogen-burning minimum mass, demonstrate a stable main sequence
   rather than classifying a cooling object as a short-lived star.
3. Converge those endpoint ages themselves using spatial and time refinements,
   atmosphere/opacity/EOS resolution controls, repeated calculations, restart
   checks, and independent initial models. Archive the executable, exact
   inputs and source artifacts. Separate numerical sensitivity from physical
   uncertainty; compare with independent stellar codes where their validated
   microphysics overlaps. Match mass and age conventions as well as abundance
   labels. Ember currently uses conserved baryonic mass in Newtonian gravity;
   it records nuclear rest-mass release but does not update gravitating mass
   from that release. The 2-trillion-year integral is .001810 of baryonic mass.
   This approximation must be assessed before a precise cross-code lifetime
   comparison; it is not part of the mesh/timestep error estimate.
4. Generalize the driver and generate mixture-specific source families for a
   pilot mass/composition survey. No unsupported extrapolation or solar-opacity
   scaling substitutes for a source calculation. Check the high-Z EOS and
   opacity against independent calculations and appropriate physical limits;
   source availability alone does not establish physical accuracy.

## Metallicity does not uniquely specify composition

The target is initially a lifetime surface in M, X, Y and the metal abundance
pattern. A curve labelled lifetime(M,Z) selects a helium-enrichment relation
Y(Z) and relative metal abundances, with X=1-Y-Z. Alternative enrichment
histories and metal patterns therefore belong in the physical sensitivity
study, especially at high Z; a solar pattern is only one controlled choice.

Adams & Laughlin's [1997 review, section II.B.1](https://sites.astro.caltech.edu/ay1/RevModPhys.69.337.pdf)
illustrates a maximum near Z=.04 under approximate opacity/homology scalings
and Y=Yp+2Z. Greater opacity competes with diminished hydrogen fuel and
composition-driven changes in luminosity. That estimate motivates a numerical
test. Neither a maximum nor its position is imposed on Ember. The eventual
survey must bracket any measured maximum and show how it shifts with mass,
helium enrichment, metal pattern, and microphysical uncertainty.
