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
through **3.615 Tyr**, at **3497 K** and surface hydrogen **X = 0.2777**.
A radiative shell first forms at **3.558 Tyr**; the centre subsequently becomes
radiative. At the plotted endpoint the radiative interior contains **72.32%**
of the mass beneath a convective envelope. The atmosphere matching depth
remains at optical depth 100 throughout. The retained sequence has used
**17.37 CPU-hours**, including rejected timestep trials.
[Current trajectory and inputs](docs/reports/2026-09-27/pms_figure_inputs.json),
[internal structure](docs/results/radiative_interior_sept27_v2.json).

The selected physics includes initial deuterium burning, pp reactions and
explicit C12/C13/N14 conversion with Solar Fusion III rates, plasma neutrino
losses, composition-dependent thermodynamics and opacities, wavelength-dependent
atmospheres, hot microscopic H/He/metal diffusion, and conservative heat
transport during composition changes. Convection mixes each connected region
instantaneously. This is not a full CNO network; grains, finite convective
mixing and cool radiative microscopic transport are not selected.
[Physics and remaining work](docs/COLD_REMNANT.md),
[driver configuration](docs/LIFETIME_DRIVER.md).

Every accepted interval passes full-step/two-half-step accuracy checks and
isotope and energy audits. Solves reuse the latest composition and can predict
the next structure. Optional local Taylor reuse of collision responses retains
the table limits and conservative exchange. Together, these changes reduced
a matched **1 Gyr** calculation from **297.8 to 81.76 CPU seconds** (**3.642×**),
with very small structural and fuel differences. Full/two-half timestep checks
remain selected. The physics and accuracy checks were unchanged; this timing
comparison does not establish the speed of a complete track.
[Solver comparison](docs/results/solver_starting_guesses_sept27_v1.json),
[collision reuse](docs/results/collision_reuse_sept27_v1.json),
[time accuracy](docs/results/envelope_time_accuracy_sept27_v1.json).

The atmosphere and EOS inputs cover all 67 saved settling comparison
structures. The selected high-gravity atmosphere table covers **4600–6000 K**
and **log g = 5.9–6.2**, using checked source solutions and explicit
interpolation. Its extension preserves all 5999 saved atmosphere queries,
three stellar intervals and the actual continuation's starting state exactly.
The boundary remains gas-only, with declared trace-helium and composition limits.
[Atmosphere checks](docs/results/highg_atmosphere_extension_sept27_v2.json),
[EOS coverage](docs/results/eos_lifetime_coverage_sept27_v1.json).

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


The lifetime solver reuses nearby EOS and nuclear-screening responses and
parallelizes independent zone and face calculations. A matched 5-Gyr
shell-burning segment on the M4 Max took 31.53 wall seconds with eight threads
versus 54.95 with the preceding code at four threads. All full/two-half
timestep checks were retained, with matching stellar results.
[Settings, measurements, and limits](docs/LIFETIME_DRIVER.md).
