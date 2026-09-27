# The 0.1-solar-mass calculation through remnant cooling

## Current priority — September 27, 2026

The furthest completed continuous Hayashi-started sequence reaches **2.102 Tyr**.
The successor with hot microscopic transport is running; consult the
[current handoff](../HANDOFF.md) and live history for its age. It remains fully
convective in the latest check. The common program now also supports recovered
metal-dependent atmospheres and the existing hydrogen/metal opacity extensions.
All **76 tests pass** in working and isolated checkouts, including an actual
radiative-core evolution, conservation checks and an exact restart.

Local composition solve precision follows the actual mixing partition: tighter
for the single fully mixed region, practical precision for stratified diffusion.
Physical inventory/energy checks are unchanged. A fresh-start negative control
with a uniformly loose tolerance accepted no evolved model; it is not progress
along the track. The replacement input set passes its initial PMS checks.

Remaining preparation is the continuous H-rich atmosphere interpolation and
smooth WD boundary. Fable is recovering old source inputs, computing tau=100
high-gravity hydrogen atmospheres, and preparing the shorter paper. The separate
flash histories inherit a forced atmosphere change; atmosphere-independent
ignition is still unestablished.

[Atmosphere integration](results/atmosphere_lifetime_atlas_sept27_v1.json),
[real diffusion/restart check](results/screened_radiative_core_sept27_v1.json),
[opacity support](results/opacity_lifetime_extension_sept27_v1.json),
[completed review](research/fable/ADVERSARIAL_REVIEW_20260926.md).

## Objective and next milestones

Use one modular evolution program, with the same physical state and conservation
accounting, from Hayashi contraction through main-sequence burning and WD cooling.
The immediate scientific test is whether the helium-3 instability occurs without
the abrupt atmosphere adjustment used in the separate late-evolution experiments.
The program should also support neighboring masses, including the boundary
between stars and brown dwarfs.

1. Finish a continuous sequence through the first radiative core with H/He and
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
