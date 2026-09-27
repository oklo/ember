# The 0.1-solar-mass calculation through remnant cooling

## Current priority — September 27, 2026

Continue the same Hayashi-started star through hydrogen exhaustion and onto
white-dwarf cooling, testing whether the helium-3 pulse occurs with a consistent
atmosphere history. The published trajectory reaches **3.615 Tyr**. Its centre
is radiative, beneath a convective envelope; **72.32%** of the mass is radiative
at the plotted endpoint. The first radiative region formed as a shell at
**3.558 Tyr**. The live calculation has continued farther; read the
[current handoff](../HANDOFF.md) before starting or changing jobs.
[Trajectory and inputs](reports/2026-09-27/pms_figure_inputs.json),
[first radiative region](results/continuous_radiative_core_sept27_v1.json),
[structure at the plotted endpoint](results/radiative_interior_sept27_v2.json).

The common lifetime program includes initial deuterium, pp and explicit CN
burning, plasma neutrino losses, composition-dependent EOS and opacity,
wavelength-dependent atmospheres, hot H/He/metal diffusion, and conservative
heat transport during composition changes. Each connected convective region
is mixed instantaneously, with a check on the gradient needed to carry its
species flux. Finite mixing is available in the engine but is not selected
in this trajectory. It is now selectable in the same lifetime program, with
[conservation and restart checks](results/finite_mixing_driver_sept27_v1.json).
Unsupported cool radiative transport remains a limitation;
there is no silent replacement for the missing microscopic law.

All accepted full and half intervals must pass isotope and energy audits.
For the present radiative-interior phase, a matched **1 Gyr** comparison
supports absolute species time accuracy **1e-4**. Compared with 1e-5 it uses
**18.03%** less CPU, changes total helium-3 by **0.000184%** and surface
temperature by less than **0.001 K**. No cell is excluded from the time-error
estimate. A separate **1 Gyr** comparison
supports local Newton corrections of **1e-12**, while retaining **1e-13** outer
composition coupling and **1e-14** integrated species balance. It uses
**27.12%** less CPU and changes total helium-3 by **0.0002092%**. These controls
do not establish time convergence through a flash.
[Time accuracy](results/envelope_time_accuracy_sept27_v1.json),
[local solve](results/local_correction_comparison_sept27_v2.json).

The atmosphere match remains at optical depth **100** throughout. The selected
high-gravity grid covers **4600–6000 K**, **log g = 5.9–6.2**, using checked
source columns and declared interpolation. The extension preserves all 5999
saved atmosphere queries and three stellar intervals exactly, including the
final physical checkpoint. Its actual continuation preserves the starting
physical state. All 67 saved late comparison structures remain supported.
Those separate histories inherit an atmosphere adjustment and cannot establish
atmosphere-independent ignition. Gas-only and trace-helium assumptions remain
explicit, and additional coverage is being prepared in advance.
[Atmosphere checks](results/highg_atmosphere_extension_sept27_v2.json),
[EOS coverage](results/eos_lifetime_coverage_sept27_v1.json).

Exact reuse of repeated EOS evaluations preserves six stellar intervals and
the final physical checkpoint, using about **5%** less CPU in one timing pair.
Six EOS, transport and restart tests pass. No interpolation or source-domain
checks change.
[EOS comparison](results/exact_eos_reuse_sept27_v1.json).

The [working paper](reports/2026-09-27/ember_status_and_future.pdf) is **30 pages**,
with the physical qualifications, model comparisons and lifetime timeline retained.
The [completed review](research/fable/ADVERSARIAL_REVIEW_20260926.md) records
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
