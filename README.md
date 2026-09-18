# Ember

Ember is a one-dimensional stellar evolution code written in C++23. It follows
the structure, composition and energy transport of very low-mass stars over
trillions of years.

[Read the current working paper](docs/reports/2026-09-18/ember_status_and_future.pdf)
([LaTeX and figure files](docs/reports/2026-09-18/README.md)).

## Research status

The calculation with microscopic diffusion and metal settling reaches
**4.002 trillion years**, with effective temperature **4612 K** and central
hydrogen **X = 3.528e-8**. Hydrogen burning continues outside the depleted core,
and supplies **97.28%** of the surface luminosity. Both luminosity and effective
temperature have passed local maxima.

A continuation using a white-dwarf atmosphere develops a helium-3 burning pulse.
Finite-rate mixing and consistent mass volumes retain ignition on three spatial
grids, but the pulse's strength and timing depend on resolution. The finest
diagnostic calculation reaches **10.27 thousand years** from the model before shell
convection. Its surface is helium enriched, with hydrogen mass fraction
**X = 0.8159** and effective temperature **4765 K**. Nuclear power is
**6.682e34 erg/s**; the shell remains active. The flash maximum, shutdown
and subsequent cooling remain unresolved.

Thirty solved hydrogen–helium atmosphere columns cover **4500–5500 K**,
**X = 0.78–0.9955** and **log g = 5.55–5.80**. Independent calculations check
interpolation, the helium-isotope approximation and lower-boundary sensitivity.
The pulse continuation uses the calculated composition response on a published
hydrogen boundary; its atmosphere dependence and spatial convergence remain
uncertainties. A separate gradual-mixing sensitivity reaches **3004 years**
from the pulse starting model; its merger timing and energy remain uncertain.

A fixed-duration timestep comparison supports retained steps of at most
**25 years** during the more gradual pulse evolution, with **0.5%** nuclear-energy
refinement and the existing structure, composition and conservation checks.
A short interior mixing readjustment also passes a two-step/four-step
comparison: deposited nuclear energy agrees to **0.006375%**. The paper states
the comparisons and their limits. Current work follows the active
burning shell toward white-dwarf cooling and tests the remaining physical
approximations. Comparisons with the published LBA97 result and MESA distinguish
the physical assumptions of each track.

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
