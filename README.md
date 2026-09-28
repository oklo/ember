# Ember

Ember is a one-dimensional stellar evolution code written in C++23. It calculates the structure, composition and energy transport of very low-mass stars.

[Read the working paper](docs/reports/2026-09-27/ember_status_and_future.pdf) ([LaTeX, figures and plotting data](docs/reports/2026-09-27/README.md)).

## Research status — September 27, 2026

The objective is one continuous calculation from a Hayashi starting model,
through hydrogen burning and helium-white-dwarf cooling, to extreme cold and
conditional disappearance through nucleon decay. The same program should
resolve the stellar/brown-dwarf boundary as initial mass and composition vary.
The full calculation is not yet established.

The common `ember-evolve --lifetime` program carries the Hayashi-started star
through **4.006 Tyr**, at **4235 K** and surface hydrogen **X = 0.9922**.
A radiative shell first forms at **3.558 Tyr**; the centre subsequently becomes
radiative. The convective envelope now contains **2.214%** of the mass.
Shell burning supplies **25.50%** of the luminosity and is declining;
no helium-3 runaway has occurred in this continuous calculation. The atmosphere
matching depth remains at optical depth 100. A continuous composition overlap now includes mixed hydrogen–helium atmospheres.
Solved mixed-atmosphere cells are available down to **4000 K** at the current
gravity; the selected continuation uses the table through **4200 K**. The latest continuation
advanced **20.48 Myr** in **34.06 seconds**, including failed trials at a dense-envelope
EOS limit. Cooler atmosphere columns and a consistent dense-hydrogen EOS treatment are being assessed. The paper currently shows the track through **4560 K** and
accounts for **31.30 CPU-hours** in that plotted sequence.
[4000 K atmosphere checks](docs/results/cold_atmosphere_4000_sept28_v1.json),
[Current model and envelope checks](docs/results/cooling_envelope_4235_sept27_v1.json),
[plotted trajectory](docs/reports/2026-09-27/pms_figure_inputs.json).

The selected physics includes initial deuterium burning, pp reactions and
explicit C12/C13/N14 conversion with Solar Fusion III rates, plasma neutrino
losses, composition-dependent thermodynamics and opacities, wavelength-dependent
atmospheres, hot microscopic H/He/metal diffusion, and conservative heat
transport during composition changes. Convection mixes each connected region
instantaneously. This is not a full CNO network; grains, finite convective
mixing and cool radiative microscopic transport are not selected.
[Physics and remaining work](docs/COLD_REMNANT.md),
[driver configuration](docs/LIFETIME_DRIVER.md).

Every accepted interval retains full-step/two-half-step accuracy checks and
isotope and energy audits. The eight-thread solver reuses nearby EOS and
collision responses, predicts structure and includes local burning feedback
in its structure matrix. Richardson extrapolation is selected in the current
cooling phase, with a fallback to the half-step solution when its physical
checks fail. Over the same **20 Myr**, the selected settings reduced wall time
from **26.56 to 19.25 seconds**, with all measured state differences below
**0.1%** against a smaller-step reference. The comparison uses two threads;
it does not measure a complete track or validate flash onset.
[Cooling comparison](docs/results/richardson_cooling_sept27_v1.json),
[solver settings and earlier measurements](docs/LIFETIME_DRIVER.md).

The accepted structure is within EOS support, but cooler trial states encounter
excluded source values. Direct comparisons show substantial differences between
FreeEOS and independent dense-hydrogen models in this region; more source points
alone may not resolve the limitation. Completed mixed hydrogen–helium
atmosphere cells cover **4200–4550 K**, **log g = 6.1–6.3**, and hydrogen
fractions **0.98–0.9955**; the **4200–4400 K** interval also reaches **log g = 6.5**.
An independent atmosphere calculation inside the table measures a maximum
**0.9565%** pressure interpolation difference. Interpolation uses only supported
cells; the boundary remains gas-only with an explicit helium-isotope approximation.
At low metal abundance, an optional quadratic free-energy interpolation avoids
unnecessary dependence on the most metal-rich source planes. Narrower source
derivative stencils recover additional cool-envelope support without relaxing
source acceptance or changing previously supported table values.
[Cold atmosphere checks](docs/results/cold_atmosphere_4400_sept27_v1.json),
[low-metal EOS comparison](docs/results/low_metal_eos_sept27_v1.json).

Separate late-evolution calculations develop a helium-3 shell pulse after an
atmosphere adjustment and then turn toward cooling. They establish a possible
instability under that boundary treatment, not atmosphere-independent ignition.
The continuous calculation tests whether the pulse survives a consistent history.
[Working paper](docs/reports/2026-09-27/ember_status_and_future.pdf),
[Fortran comparison](docs/F77_LBA97_COMPARISON.md).

The working paper is **30 pages**, reduced from 63 while retaining the
physical qualifications and lifetime timeline. The large local datasets are
required for the full test suite; public runtime/test-data availability remains
incomplete.

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
