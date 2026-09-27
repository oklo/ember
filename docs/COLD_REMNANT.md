# The 0.1-solar-mass calculation through remnant cooling

## Current priority — September 27, 2026

Continue the same Hayashi-started star through hydrogen exhaustion and onto
white-dwarf cooling, testing whether a helium-3 pulse occurs with a consistent
atmosphere history. The current model has reached **4.006 Tyr** and **4401 K**.
Convection occupies **2.214%** of the mass in an outer envelope; hydrogen
burning peaks at **m/M = 0.9704** and supplies **46.24%** of the surface
luminosity. Its power is declining and no helium-3 runaway has occurred.
The first radiative region formed at **3.558 Tyr**. Read the
[current handoff](../HANDOFF.md) before starting or changing jobs.
[Current state and opacity response](results/dense_envelope_density_extension_sept27_v1.json).

The common lifetime program includes initial deuterium, pp and explicit CN
burning, plasma neutrino losses, composition-dependent EOS and opacity,
wavelength-dependent atmospheres, hot H/He/metal diffusion, and conservative
heat transport during composition changes. Each connected convective region
is mixed instantaneously, with a check on the gradient needed to carry its
species flux. Finite mixing is available in the same program but is not
selected in this trajectory. Unsupported cool radiative transport remains a
limitation; no missing microscopic law is silently replaced.
[Finite-mixing checks](results/finite_mixing_driver_sept27_v1.json).

All accepted full and half intervals pass isotope and energy audits. The
current cooling calculation also uses checked Richardson extrapolation.
Across a **20 Myr** comparison, its selected settings reduce wall time from
**26.56 to 19.25 seconds**, with differences below **0.1%** against a
smaller-step reference. A mixing-step limit compares new composition against
the previous composition averaged over the current instantaneous mixed regions,
so prescribed homogenization does not force arbitrarily small steps. Burning
and transport changes still count; finite mixing retains its pointwise limit.
These controls do not establish convergence through a flash.
[Time accuracy](results/richardson_cooling_sept27_v1.json),
[mixing limit](results/instantaneous_mixing_step_cap_sept27_v1.json).

The atmosphere match remains at optical depth **100** throughout. Solved
mixed hydrogen–helium columns cover **4400–4550 K** at **log g = 6.1–6.3**,
with hydrogen fractions **0.98–0.9955**. A continuous composition overlap
carries the star through the former **0.5% helium** limit. Only cells with
complete source support are interpolated; cooler columns are being calculated.
The EOS supports the accepted state. A tested interior-opacity extension permits
further cooling; completed atmosphere coverage currently ends at **4400 K**. Gas-only and helium-isotope approximations remain explicit. The separate flash histories retain their atmosphere
adjustment and cannot establish atmosphere-independent ignition.
[Cold source checks](results/cold_atmosphere_4400_sept27_v1.json).

In the dense hydrogen envelope, the opacity approximation uses the measured
hydrogen slope of the retained source while remaining within its temperature
and density coverage. The selected domain extends through **log R = 2.5**
and **log T = 6.3**. At the latest density limit, convection carries **99.95%**
of the thermal luminosity. Matched **1 Myr** half/double-opacity controls
produce negligible structural changes. This supports the local approximation without establishing opacity
accuracy in radiative layers.
[Opacity response](results/dense_envelope_density_extension_sept27_v1.json).

The [working paper](reports/2026-09-27/ember_status_and_future.pdf) is **30 pages**
and currently plots the trajectory through **4560 K**. It retains the
physical qualifications, model comparisons and lifetime timeline. The
[completed review](research/fable/ADVERSARIAL_REVIEW_20260926.md) records
remaining code, reproduction and physics limitations.

## Objective and next milestones

Use one modular evolution program, with the same physical state and conservation
accounting, from Hayashi contraction through main-sequence burning and WD cooling.
The immediate scientific test is whether the helium-3 instability occurs without
the abrupt atmosphere adjustment used in the separate late-evolution experiments.
The program should also support neighboring masses, including the boundary
between stars and brown dwarfs.

1. Continue through the growth of the radiative interior with H/He and
   metal diffusion. Prepare atmosphere and interior table coverage in advance;
   preserve explicit source limits and verify overlap values and derivatives.
2. Continue through hydrogen exhaustion with a consistent joining depth between
   interior and atmosphere. Treat any rapid burning episode with suitable finite
   mixing and time resolution, conserving species and energy. Follow both the
   nuclear power and the remaining helium-3 reservoir before declaring final cooling.
3. Establish the WD cooling solution and add the physical processes appropriate
   to its composition and density: diffusion, conductive transport, dense-matter
   thermodynamics, crystallization and envelope chemistry. Grain opacity and its
   thermodynamic consequences must be coupled where relevant; the selected stellar
   atmosphere calculations currently contain gas opacity only.
4. At very low luminosity, compare isolation with evolving environments, including
   accretion, unbound tidal encounters and conditional dark-matter heating. These
   are scenario calculations, not universally imposed heating floors.
5. Extend the lifetime account beyond 100 K under explicit nucleon-decay models,
   retaining stable baryons as an alternative. An assumed particle lifetime must
   remain separate from established stellar physics. Use reduced descriptions
   where uncertainty in the physical inputs outweighs the benefit of resolving
   additional stellar detail.

The nuclear network presently includes pp burning and explicit CN conversion,
not a complete CNO network. Microscopic metal transport uses a common GS98 metal
response with a fully ionized collision approximation in its declared hot domain.
Convection homogenizes each connected region in the selected continuous run;
finite mixing and other envelope processes require appropriate selection and
checks before their effects can be claimed. Unsupported physical regimes stop
the calculation rather than silently acquiring a different atmosphere or heat law.

## Supporting documents

- [Atmosphere physics and coverage](ATMOSPHERE.md) and
  [white-dwarf atmosphere plan](WHITE_DWARF_ATMOSPHERES.md).
- [Dense equation of state](DENSE_EOS.md), [grains](GRAINS.md), and
  [neutrino losses](NEUTRINOS.md).
- [Very cold remnant and environmental physics](ULTRACOLD_PHYSICS.md) and
  [appropriate levels of approximation](APPROXIMATION_REMIT.md).
- [Computational cost](COMPUTATIONAL_COST.md),
  [data reproduction](DATA_REPRODUCTION.md), and
  [common driver configuration](LIFETIME_DRIVER.md).

The current job and review instructions belong in [HANDOFF.md](../HANDOFF.md).
Dated numerical evidence is retained in `docs/results`; the working paper
presents scientific conclusions and their limitations.
