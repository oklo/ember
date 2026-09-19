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

The **2115-point** late sequence reaches **143.4 thousand years**, **5986 K**
and nuclear power **2.508e32 erg/s**, supplying **93.31%** of photon luminosity.
A matched **2000-year** mass-grid comparison changes deposited nuclear energy
by **0.5220%** and endpoint power by **0.5687%**. Half the burning power
occupies **39** cells instead of **12**. This supports the later shell burning;
the full pulse's onset and peak remain uncertain.

An **800/400/200-year** comparison supports larger time steps with unchanged
accuracy and conservation checks. The **17,200-year** segment takes
**31.21 CPU minutes**, including timestep comparisons and rejected trials.
A newly radiative envelope layer is included by extending microscopic species
transport to **1.500 MK**, supported by source ionization, actual exchange and
finite-step checks. This local approximation does not describe neutral-fluid
transport in a cold remnant.

The separate onset calculation reaches **2374 years**, **4606 K** and nuclear
power **2.610e34 erg/s**. Power reaches a local maximum of **3.719e35 erg/s**
near **1642 years**, then declines and rises again. The most recent **5.756-year**
convective adjustment passes an **8/16-step** comparison, with nuclear energy
differing by **0.2976%** and the same final convection. A complete, spatially
converged pulse and white-dwarf cooling sequence remain unresolved.

A matched **12,000-year** comparison of two atmosphere references reaches
**5701 K**. Effective temperatures differ by **0.2296 K** and surface luminosities
by **0.01569%**; deposited nuclear energies differ by **9.917e-6%**, with the
same final convective regions. This tests the visited interval without
bounding the shared atmosphere uncertainty.

Thirty solved hydrogen–helium atmosphere columns cover **4500–5500 K**,
**X = 0.78–0.9955** and **log g = 5.55–5.80**. Independent calculations check
interpolation, the helium-isotope approximation and lower-boundary sensitivity.
The local mixed-atmosphere grid now contains **54** solved columns covering
**4700–5500 K** and **log g = 4.95–5.80**. The absolute hydrogen reference
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
A matched **10,000-year** comparison below **log g = 5.10** reaches
**log g = 5.055**. The two reference prescriptions differ by **0.2991 K** and
**0.02206%** in surface luminosity, with the same convective regions. This is
sensitivity to the two prescriptions, not a bound on their shared uncertainty.
Nine additional **6000 K** mixture columns and a source at the published
hydrogen grid's gravity floor pass independent source and depth checks.
The local mixed-atmosphere grid contains **45** solved columns over
**5000–6000 K** and **log g = 4.65–5.25**, with independent depth checks
and unchanged matching states in every overlapping cell.
A matched **8000-year** comparison of two temperature extensions of the hydrogen
reference reaches **5553 K**. Their effective temperatures differ by **4.301 K**,
surface luminosities by **0.3070%**, and deposited nuclear energy by **5.576e-5%**,
with the same convective regions. This measures sensitivity over the visited interval.
The pulse continuation uses the calculated composition response on a published
hydrogen boundary; its atmosphere dependence and spatial convergence remain
uncertainties. A separate gradual-mixing sensitivity reaches **3004 years**
from the pulse starting model; its merger timing and energy remain uncertain.

A tested timestep method reduces the stellar solves from **120** to **81**
over **2000 years**, including one rejected larger step and its smaller retry.
Released nuclear energy differs by **0.1353%** and endpoint nuclear power by
**0.08057%**. Physical convective regions agree at all three comparison ages.
A **400-year** comparison using one, two and four steps also passes the
unchanged checks: the coarsest and finest paths differ by **0.3543%** in
released nuclear energy and **0.005697%** in endpoint nuclear power.
The full-step ceiling is **400 years**, with retained halves at most **200 years**;
the existing local error checks still require shorter steps where needed.
This supports the declining pulse phase.
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
