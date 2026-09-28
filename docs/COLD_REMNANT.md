# The 0.1-solar-mass calculation through remnant cooling

## Current calculation — September 28, 2026

The common lifetime program follows Hayashi contraction, main-sequence burning
and helium-white-dwarf cooling. The reference reaches **4235 K** without a late
helium-3 runaway. A fresh calculation reproduces its luminosity to **0.08179%**
and total He3 to **0.8907%**, using **5.361 CPU-hours** with tables prepared.
Both reference calculations stop at the same dense-envelope EOS limit.
[Fresh calculation](results/fresh_hayashi_cooling_sept28_v1.json).

A trial with a revised cold-envelope EOS present from Hayashi initialization
reaches **4.009 Tyr**, **3728 K**, and **0.0001943 Lsun**. Shell burning supplies
**5.200%** of luminosity and peaks at **m/M = 0.9693**; convection occupies an
outer **2.214%** of the mass. No late He3 runaway occurs. The latest continuation
adds **38.46 Myr** in **48.10 seconds**, including failed trials at its endpoint.
This is provisional: the cold-EOS join distorts heat capacity and the adiabatic
gradient, and the outer mesh also differs from a resolved envelope integration.
A direct-source outer boundary is being tested. The final solver failure occurs
near an EOS support gap; a direct EOS exception is not established by that log.
[Trial structure and limitations](results/cooling_trial_3728_sept28_v1.json).

The trial selects an optional quantum nuclear-screening correction at **3965 K**;
it is therefore not yet an unchanged-physics calculation from Hayashi contraction.
A matched **50 Myr** comparison changes Teff by **0.02668 K** and luminosity by
**0.003276%**, with unchanged convective mass. The correction retains explicit
validity limits; pycnonuclear burning remains unsupported.
[Screening and transport checks](results/quantum_screening_sept28_v1.json).

The code includes initial deuterium, pp and explicit CN burning, plasma neutrino
losses, composition-dependent EOS and opacity, wavelength-dependent atmospheres,
hot H/He/metal diffusion, and conservative heat transport during composition
changes. Connected convective regions mix instantaneously in this trajectory.
Finite mixing is available in the same program. Full-step/two-half-step checks,
isotope and energy audits remain active; checked Richardson extrapolation falls
back when its physical checks fail. A matched **20 Myr** comparison reduces wall
time from **26.56 to 19.25 seconds** with differences below **0.1%** against a
smaller-step reference. This does not establish convergence through a flash.
[Time accuracy](results/richardson_cooling_sept27_v1.json),
[mixing step limit](results/instantaneous_mixing_step_cap_sept27_v1.json),
[finite mixing](results/finite_mixing_driver_sept27_v1.json).

The atmosphere joins at optical depth **100**. Mixed H/He coverage reaches
**3200 K**, with **log g = 6.3–6.7** and **X = 0.98–0.9955** over the cold cells.
Independent **3300 K** columns differ from interpolation by at most **0.3912%**
in pressure and **0.009569%** in temperature; a deeper **3200 K** column changes
matching pressure by **0.02925%**. Old supported values remain unchanged.
The cooling trial uses the extended table. Gas-only and helium-isotope
approximations remain explicit.
[Atmosphere coverage](results/cold_atmosphere_3200_sept28_v1.json),
[source optimization](results/tlusty_russel_sept28_v1.json).

The low-metal EOS uses a quadratic free energy joined smoothly to the original
composition interpolation. Direct source calculations extend its supported
domain without accepting failed source states. The experimental cold-envelope
replacement still needs physical validation; energy consistency alone does not
establish the accuracy of its heat response.
[Low-metal EOS comparison](results/low_metal_eos_sept27_v1.json).
Dense radiative opacity uses the retained source's hydrogen slope within its
actual temperature/density coverage. Over **200 Myr** of cooling to **4035 K**,
halving or doubling it changes luminosity by at most **0.06752%**, with unchanged
convective mass. That supports this approximation locally in the convective
envelope, not in a future radiative envelope.
[Opacity comparison](results/dense_cooling_opacity_sept28_v1.json).

Separate flash histories retain an atmosphere adjustment before ignition and
cannot establish atmosphere-independent ignition. The [working paper](reports/2026-09-27/ember_status_and_future.pdf)
remains **30 pages** and currently plots the continuous trajectory through
**4560 K**; its figures lag these calculations. Current jobs and continuation
instructions are in [HANDOFF.md](../HANDOFF.md).

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
