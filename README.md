# Ember

Ember is a one-dimensional stellar evolution code written in C++23. It follows
the structure, composition and energy transport of very low-mass stars over
trillions of years.

[Read the current working paper](docs/reports/2026-09-19/ember_status_and_future.pdf)
([LaTeX and figure files](docs/reports/2026-09-19/README.md)).

## Research status

The calculation with microscopic diffusion and metal settling reaches
**4.002 trillion years**, with effective temperature **4612 K** and central
hydrogen **X = 3.528e-8**. Hydrogen burning continues outside the depleted core,
and supplies **97.28%** of the surface luminosity. Both luminosity and effective
temperature have passed local maxima.

A continuation using a white-dwarf atmosphere develops a helium-3 burning pulse.
Finite-rate mixing and consistent mass volumes retain ignition on three spatial
grids, but the pulse's strength and timing depend on resolution. The 1695-point
diagnostic calculation reaches **13.45 thousand years** from the model before shell
convection. Its surface is helium enriched, with hydrogen mass fraction
**X = 0.8159** and effective temperature **4848 K**. Nuclear power is
**2.539e34 erg/s**; the shell remains active. The flash maximum, shutdown
and subsequent cooling remain unresolved.

A local refinement to **1923** mass points distributes half the burning power
over **18** cells instead of **one** in a matched **688.5-year** comparison.
Deposited nuclear energy is **55.17%** lower and endpoint nuclear power is
**59.09%** lower, while surface temperatures differ by only **0.6349 K**.
An initial structural readjustment means this measures both resolution and
restart effects; it does not establish spatial convergence.

A further refinement to **2010** mass points changes released nuclear energy
by **1.468%** and endpoint nuclear power by **1.980%** over a matched
**978.3-year** interval. This supports the late burning segment; the earlier
onset and peak remain uncertain.

The **2010-point** sequence reaches **57.47 thousand years**, **5059 K** and
nuclear power **8.376e32 erg/s**. Power declines by **75.31%** during the
additional **3.516e4 years**. All **196** intervals pass independent timestep,
composition and energy checks. Shell burning continues.

Thirty solved hydrogen–helium atmosphere columns cover **4500–5500 K**,
**X = 0.78–0.9955** and **log g = 5.55–5.80**. Independent calculations check
interpolation, the helium-isotope approximation and lower-boundary sensitivity.
The local mixed-atmosphere grid now contains **45** solved columns covering
**4700–5500 K** and **log g = 5.10–5.80**. The absolute hydrogen reference
below **log g = 5.5** remains an explicitly inferred continuation. A matched
**10,000-year** calculation testing the extension below **log g = 5.25**
reaches **log g = 5.231**. The two references change surface temperature
by **1.476 K** and luminosity by **0.1171%**, with released nuclear energy
differing by **7.675e-6%** and the same convective regions. A separate matched
**1000-year** comparison at **log g = 5.335–5.331** changes surface temperature
by **11.27 K** and luminosity by **0.9186%**; released nuclear energy differs
by **0.02255%**, with the same convective regions. This tests the region
occupied by the model without providing an uncertainty bound over the full grid.
Independent **4800 K** sources
check the helium correction to **0.09612%** in temperature and **0.02421%** in
pressure at that test point.
The pulse continuation uses the calculated composition response on a published
hydrogen boundary; its atmosphere dependence and spatial convergence remain
uncertainties. A separate gradual-mixing sensitivity reaches **3004 years**
from the pulse starting model; its merger timing and energy remain uncertain.

A tested timestep method reduces the stellar solves from **120** to **81**
over **2000 years**, including one rejected larger step and its smaller retry.
Released nuclear energy differs by **0.1353%** and endpoint nuclear power by
**0.08057%**. Physical convective regions agree at all three comparison ages.
The full-step ceiling is **200 years**; the existing local error checks still
require shorter steps where needed. This supports the declining pulse phase.
A short interior mixing readjustment also passes a two-step/four-step
comparison: deposited nuclear energy agrees to **0.006375%**. The paper states
the comparisons and their limits. Current work follows the active
burning shell toward white-dwarf cooling and tests the remaining physical
approximations. Comparisons with the published LBA97 result and MESA distinguish
the physical assumptions of each track.

This publication updates the paper and README; the figures are unchanged. The checked-in source is
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
