# Ember

Ember is a one-dimensional stellar evolution code written in C++23. It calculates the structure, composition and energy transport of very low-mass stars.

[Read the working paper](docs/reports/2026-09-27/ember_status_and_future.pdf) ([LaTeX, figures and plotting data](docs/reports/2026-09-27/README.md)).

## Research status — September 28, 2026

The objective is one continuous calculation from Hayashi contraction through
hydrogen burning, helium-white-dwarf cooling and extreme cold, with conditional
nucleon-decay scenarios. Varying the initial mass and composition should also
resolve the stellar/brown-dwarf boundary. The full calculation remains unfinished.

The common `ember-evolve --lifetime` program follows an initially
**0.1-solar-mass** star from Hayashi contraction onto helium-white-dwarf
cooling, without a late helium-3 runaway. A fresh calculation reaches
**4235 K** at **4.007 Tyr** and reproduces the reference cooling luminosity
within **0.08179%** and total He3 within **0.8907%**, using **5.361 CPU-hours**
with the physics tables prepared.
[Fresh-track comparison](docs/results/fresh_hayashi_cooling_sept28_v1.json).

A development calculation with a source-integrated outer envelope reaches
**3069 K** at **4.010 Tyr**. Nuclear burning supplies **0.4510%** of luminosity;
the core remains at **2.750 MK**. Its latest **300 Myr** took **13.15 seconds**
on four CPU threads, with no rejected intervals. The envelope treatment was
introduced near **3757 K**, so this continuation is provisional: the final
boundary must still be applied in a fresh calculation from Hayashi contraction.
The leading quantum-ion heat-capacity correction is estimated at **0.1754%**
when integrated over this structure; it is not yet selected in the calculation.
[Cooling boundary checks](docs/results/cooling_boundary_sept28_v1.json).

The selected physics includes deuterium and pp burning, explicit C12/C13/N14
conversion with Solar Fusion III rates, plasma neutrino losses,
composition-dependent EOS and opacities, wavelength-dependent atmospheres,
hot H/He/metal diffusion, and conservative heat transport during composition
changes. Connected convective regions mix instantaneously. The network is
not full CNO; grains, finite convective mixing and cool radiative microscopic
transport are not selected.
[Physics and remaining work](docs/COLD_REMNANT.md),
[configuration](docs/LIFETIME_DRIVER.md).

Atmospheres use wavelength-dependent gas absorption and match at optical
depth **100**. The cooling calculation integrates the outer **0.15%** of
the mass beneath that atmosphere. New mixed H/He source columns reach
**2500 K**, with some edge cells still unsupported. Revised molecular
collision-induced absorption changes the stellar effective temperature by
**3.219 K** and luminosity by **0.3513%** in a matched **150 Myr** comparison.
These tables and the integrated-envelope implementation remain local research
inputs; public runtime and test-data availability is incomplete.
[Atmosphere and stellar checks](docs/results/cooling_boundary_sept28_v1.json).
Removing repeated chemistry work reduced source CPU time by **1.601–3.455**
in two controls, with state differences below **1.628e-7**.
[Source optimization](docs/results/tlusty_russel_sept28_v1.json).

The solver uses zone threads, nearby EOS and collision responses, a structure
predictor and local burning feedback. Every accepted interval retains full-step/
two-half-step error checks and isotope and energy audits. Checked Richardson
extrapolation falls back when its physical checks fail. A matched **20 Myr**
cooling comparison reduced wall time from **26.56 to 19.25 seconds**, with
state differences below **0.1%** against a smaller-step reference. That two-thread
measurement does not establish complete-track timing or flash convergence.
[Numerical checks](docs/results/richardson_cooling_sept27_v1.json).

Separate calculations develop a helium-3 pulse after an atmosphere adjustment;
they do not establish ignition with a continuous atmosphere history. The
[working paper](docs/reports/2026-09-27/ember_status_and_future.pdf) is **30 pages**,
includes the model comparisons and lifetime timeline, and currently plots the
continuous track through **4560 K**. Public runtime and test-data availability
remains incomplete; see the reproduction instructions below.

## Building

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

Requires a C++23 compiler. The release build targets Apple M4 when supported,
uses CPU threads for independent zone work, and solves its small matrix blocks
with a pivoted kernel. The current stellar solver does not use Metal or the GPU.

New bulk physics tables and raw archives are generated locally and excluded from
Git. The repository includes their generators, source patches, configurations and
small provenance records; see [the reproduction guide](docs/DATA_REPRODUCTION.md).
A fresh clone must generate or separately restore the new input families before
running the forward stellar model or its table-dependent tests. The documented
full-suite checks use the development machine's installed inputs.

`-ffast-math` is deliberately *not* used: it licenses the compiler to assume no
NaN or infinity, and a stellar model legitimately probes states where a table
returns one. Those should surface, not be optimised away.

## Licence

MIT. The opacity, equation-of-state and conductivity data are redistributed under their
own terms; see `data/*/README.md`.
