# Ember

Ember is a one-dimensional stellar evolution code written in C++23. It follows
the structure, composition and energy transport of very low-mass stars over
trillions of years.

[Read the current working paper](docs/reports/2026-09-17/ember_status_and_future.pdf)
([LaTeX and figure files](docs/reports/2026-09-17/README.md)).

## Research status

The calculation with microscopic diffusion and metal settling reaches
**4.002 trillion years**, with effective temperature **4612 K** and central
hydrogen **X = 3.528e-8**. Hydrogen burning continues outside the depleted core,
and supplies **97.28%** of the surface luminosity. Both luminosity and effective
temperature have passed local maxima.

A continuation using a white-dwarf atmosphere develops a helium-3 burning pulse.
Finite-rate mixing and consistent mass volumes retain ignition on three spatial
grids, but the pulse's strength and timing depend on resolution. The finest
diagnostic calculation reaches **6803 years** from the model before shell
convection and stops when a trial surface composition exceeds the available
atmosphere coverage. The flash maximum, shutdown and subsequent cooling remain
unresolved.

Current work addresses hydrogen–helium atmosphere coverage, mixing across
composition gradients, nuclear burning and energy conservation during the pulse.
The paper compares Ember with the published LBA97 result and MESA calculations,
and distinguishes the physical assumptions of each track.

This publication updates the paper and its figures. The checked-in source is
an earlier development state; source and numerical input updates for the new
calculations are still being prepared.

## Aim

Follow an initially **0.1-solar-mass** star through hydrogen exhaustion and helium
white-dwarf cooling below **100 K**, then explore its later evolution under
explicit environmental and nucleon-decay assumptions. The complete cooling
trajectory has not yet been calculated.

## Background

The comparison with Laughlin, Bodenheimer and Adams (1997) uses their stated
results and approximate curves read from the published figures. The separate
[Fortran reconstruction](https://github.com/oklo/Henyey) seeks to reproduce
that calculation's assumptions and methods.

## Building

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

Requires a C++23 compiler. On Apple silicon the build tunes for the host core
and links Accelerate for LAPACK; it is not otherwise platform-specific.

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
