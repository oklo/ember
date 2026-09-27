# The 0.1-solar-mass calculation through remnant cooling

## Current priority — September 27, 2026

The continuous Hayashi-started sequence forms its first radiative region at
**3.558 Tyr**, **3250 K**, with surface hydrogen **X=.1722**. The checked
continuation passes **3.585 Tyr**, with H/He and metal
diffusion active across radiative boundaries. A radiative shell first surrounds
the convective centre; the centre then becomes radiative. About **54.23%** of the
mass is now radiative. All accepted intervals pass the
unchanged isotope and energy audits. The fixed-metal-atmosphere comparison
is complete at **3.078 Tyr**. Both use the common lifetime program.
[First radiative region](results/continuous_radiative_core_sept27_v1.json),
[internal structure](results/radiative_shell_geometry_sept27_v1.json).
The updated paper and Figure 1 use the continuous history through **3.585 Tyr**.
Read the [current handoff](../HANDOFF.md) for the live process and next work.

Failed trial audits now trigger a smaller timestep, with a bounded number of
retries. Every full and half interval still has to pass the same checks before
acceptance. Controller recovery, physical diffusion and exact restart pass
their tests in the isolated checkout. Local composition tolerances remain
1e-15 for one instantaneously mixed region and a checked 1e-13 setting otherwise;
the global inventory budget remains 1e-14. A matched 11.76-Myr comparison reduces
accepted intervals from 19 to 10, with a 0.00004235% luminosity difference.
[Composition-solve comparison](results/radiative_core_solve_accuracy_sept27_v1.json).

A matched **100 Myr** comparison supports species time-error **1e-6** in this
phase, with **19** accepted intervals instead of **31** and a **0.0002978%**
difference in total helium-3. Separating local Newton correction accuracy from
integrated conservation also removes a solver stall: a checked **100 Myr**
continuation passes in **17** steps with no rejections. Integrated species and
reconstructed continuity still use **1e-14**; local corrections use **1e-13**.
[Time accuracy](results/radiative_shell_time_accuracy_sept27_v1.json),
[diffusion convergence](results/species_correction_accuracy_sept27_v1.json).

A matched **1 Gyr** comparison supports a local Newton correction criterion ten
times the outer coupling accuracy while retaining the same conservation bounds.
It removes five stalled local solves and uses **27.12%** less CPU in this comparison;
the total helium-3 inventory differs by **0.0002092%**. Twelve existing burning,
transport, finite-mixing and restart tests pass. The same star continues with
this setting beyond **3.573 Tyr**. This is evidence for the present gradual
evolution, not for time convergence through a flash.
[Local-solve accuracy](results/local_correction_comparison_sept27_v2.json).

The continuous trajectory passes all isotope and energy checks and is beyond
**3.585 Tyr**. The atmosphere match remains at optical depth 100. The 4600 K
gravity row and 4625 K, log g = 6.1 corner now have checked source solutions.
All 5932 checked early evaluations remain exact, and all 67 saved late
structures are supported. Matching pressure changes by at most 1.101% on
those late structures; three stellar intervals and the final physical
checkpoint match exactly. The actual continuation also preserves its physical
starting state. Gas-only and trace-helium qualifications remain, and the
high-gravity source domain has not yet been extended. Separate flash histories
still inherit a forced atmosphere change.
[Atmosphere replacement checks](results/highg_atmosphere_integration_sept27_v2.json),
[hydrogen-rich boundary](results/hydrogen_envelope_integration_sept27_v1.json).

Exact reuse of repeated EOS evaluations preserves six complete stellar
intervals and the final physical checkpoint. Six EOS, transport and restart
tests pass. The matched timing pair uses about **5%** less CPU, with scheduling
variation included. The continuous calculation now uses this implementation.
[EOS comparison](results/exact_eos_reuse_sept27_v1.json).

The [September 27 paper](reports/2026-09-27/ember_status_and_future.pdf) is
**30 pages**, reduced from 63 with physical qualifications, comparisons and the
lifetime timeline retained.

[Hydrogen-rich atmosphere recovery](results/atmosphere_hydrogen_recovery_sept27_v1.json),
[real diffusion/restart check](results/screened_radiative_core_sept27_v1.json),
[opacity support](results/opacity_lifetime_extension_sept27_v1.json),
[completed review](research/fable/ADVERSARIAL_REVIEW_20260926.md).

A 50-Gyr comparison supports a species time-error tolerance of 1e-7 in the
smooth fully convective phase. It required 25 intervals instead of 73;
structure accuracy and conservation checks are unchanged. The same star
continues with this setting from 1.846 Tyr.
[Timestep comparison](results/lifetime_time_accuracy_sept27_v1.json).

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
