# The 0.1-solar-mass calculation through remnant cooling

## Current calculation — September 28, 2026

The furthest diagnostic with a source-integrated outer envelope reaches
**3069 K**, **4.010 Tyr**, and **8.045e-5 Lsun**. Nuclear burning supplies
**0.4510%** of luminosity and is declining; no late helium-3 runaway occurs.
The latest segment advances **300 Myr** in **13.15 seconds**, with **16**
accepted intervals and no rejections. It reached its requested age.

Revised molecular collision-induced absorption is selected throughout the
new cooling atmosphere table. In a matched **150 Myr** stellar comparison,
it changes effective temperature by **3.219 K** and luminosity by **0.3513%**,
with unchanged convective extent. The supporting columns span **2500–4550 K**;
some edge cells remain unsupported. The outer-envelope map has useful coverage
to about **2550–2600 K** along the expected cooling corridor.

The outer **0.15%** of the mass is integrated hydrostatically, with its mass,
composition and approximate thermal reservoir retained in the common evolution
program. The estimated thermal error during steady cooling is **0.02–0.03%**
of luminosity. The boundary was inserted near **3757 K**; that readjustment is
not established at this accuracy. A fresh Hayashi run with the final boundary
throughout is still required. [Boundary checks](results/cooling_boundary_sept28_v1.json).

A conditional radiative-opacity continuation removes the dense-core table edge
where conduction dominates. Varying that radiation opacity by factors of
**0.1** and **10** changes luminosity by at most **0.008869%** over **150 Myr**.
The source slope and all first derivatives are retained in a smooth overlap;
the allowed change in total heat conductivity is bounded by **0.1%** under the
declared opacity uncertainty. [Transport checks](results/conductive_interior_opacity_sept28_v1.json).

The common lifetime program follows Hayashi contraction, main-sequence burning
and helium-white-dwarf cooling. The reference reaches **4235 K** without a late
helium-3 runaway. A fresh calculation reproduces its luminosity to **0.08179%**
and total He3 to **0.8907%**, using **5.361 CPU-hours** with tables prepared.
Both reference calculations stop at the same dense-envelope EOS limit.
[Fresh calculation](results/fresh_hayashi_cooling_sept28_v1.json).

The trial selects an optional quantum nuclear-screening correction at **3965 K**;
it is therefore not yet an unchanged-physics calculation from Hayashi contraction.
A matched **50 Myr** comparison changes Teff by **0.02668 K** and luminosity by
**0.003276%**, with unchanged convective mass. The correction retains explicit
validity limits; pycnonuclear burning remains unsupported.
[Screening and transport checks](results/quantum_screening_sept28_v1.json).

The selected electron-collision extension increases
the degeneracy range from **64 to 128**. Independent source checks give a
maximum transport-response difference of **0.06157%**. A matched **5 Myr**
comparison changes the largest tracked global quantity by **0.0004176%**.
This extends the same collision model; it does not add correlated or relativistic
scattering. [Collision coverage](results/electron_pair_eta128_sept28_v1.json).

The selected interior EOS extension fills missing source calculations
across the composition planes and improves agreement with direct FreeEOS
queries: the largest checked heat-capacity error falls from **0.3800%** to
**0.001378%**. A matched **5 Myr** stellar comparison changes luminosity by
**0.004556%**, with unchanged convective mass. The outer cold-envelope gap
still requires separate treatment.
[EOS coverage and source comparison](results/cooling_eos_fill_sept28_v1.json).

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

The atmosphere joins at optical depth **100**. The mixed H/He table uses
**X = 0.98–0.9955** and covers the cooling gravity range. Seven independently
calculated points differ from interpolation by at most **0.0186%** in
matching temperature and **0.08572%** in pressure. Gas-only opacity and
the helium-isotope approximation remain explicit.
[Atmosphere checks](results/cooling_boundary_sept28_v1.json),
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

1. Apply one physically supported envelope treatment from Hayashi contraction,
   including its mass, composition and heat reservoir from initialization.
   Complete the molecular-opacity update along the actual trajectory and test
   the cold/interior EOS overlap where the mesh uses it.
2. Repeat the complete trajectory with that treatment, H/He/metal diffusion
   and prepared tables. Monitor helium-3 burning with suitable mixing and time
   resolution; the absence of a flash in the present sequence is evidence, not
   proof that every physically allowed history avoids ignition.
3. Continue the WD cooling solution and add the physical processes appropriate
   to its composition and density: diffusion, conductive transport, dense-matter
   thermodynamics, crystallization and envelope chemistry. Quantum-ion heat
   capacity is the next core-EOS correction; its current estimate is **0.1754%**
   of the integrated heat capacity. Grain opacity and its
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
