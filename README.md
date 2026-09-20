# Ember

Ember is a one-dimensional stellar evolution code written in C++23. It follows
the structure, composition and energy transport of very low-mass stars over
trillions of years.

[Read the working paper](docs/reports/2026-09-20/ember_status_and_future.pdf)
([LaTeX and figure files](docs/reports/2026-09-20/README.md)).

## Research status

The **0.1-solar-mass** calculation with microscopic diffusion and metal settling
reaches **4.002 trillion years**, **4612 K**, and central hydrogen **X = 3.528e-8**.
Hydrogen burning continues outside the depleted core. A continuation with a
white-dwarf atmosphere develops a helium-3 shell flash.

Two calculations examine the flash and subsequent evolution. Times below
are measured from their common model preceding shell convection.

| Calculation | Time | Effective temperature | Nuclear power |
| --- | ---: | ---: | ---: |
| Finite convective transport from flash onset, 1983 mass points | 38.79 thousand yr | 4887 K | 1.693e33 erg/s |
| Later shell burning, alternative history, 2115 mass points | 4.624 million yr | 6021 K | 1.428e31 erg/s |

**These are alternative histories, not consecutive pieces of one track.**
The later calculation includes different mixing treatments during its flash;
its subsequent accuracy checks do not validate that entire earlier evolution.
Both calculations use the supplied hydrogen atmosphere as their reference.
Whether the flash occurs under a different atmosphere prescription remains
untested. The flash appears about **0.7546 Gyr** after the atmosphere change;
that delay does not establish independence from the envelope and fuel history.

The onset calculation passes a local nuclear-power maximum of **3.719e35 erg/s**,
followed by an **81.64%** decline and renewed burning. Nuclear power declines
through the final **2433 years**, which release **1.387e44 erg** through **16**
independently checked intervals. The largest single-cell power contribution
is **3.195%**. A matched **100-year** comparison on **1927/1983 mass points**
changes nuclear heat by **8.261%** and final power by **15.24%**. Inserting
mass points produces an initial mechanical readjustment with an integrated
energy residual of **1.844e43 erg** in its first year. A controlled correction
of the initial radii leaves one-year nuclear heat and final power unchanged
within **1e-12** relative; the initial energy change remains explicit in the
budget. The full pulse's spatial accuracy, integrated energy and final
shutdown remain unresolved.

In the later calculation, radius and surface luminosity have declined
**70.19%** and **91.83%** from their saved maxima. Temperature falls from
**6399 K** by **377.6 K** over **3.536 million years**, while nuclear burning
supplies **32.66%** of the current surface luminosity. This supports a turn
toward cooling within that calculation.

The final **369,800 years** use **11.32 CPU minutes** and **3.971 minutes of
active elapsed time**, including timestep comparisons. Trial intervals are
at most **51,200 years**; the accepted solution retains the two shorter steps,
with unchanged accuracy and conservation criteria. The main track through
**4.002 Tyr** accounts for **15.13 elapsed hours** and **36.31 CPU-hours**,
excluding source construction and pauses between jobs. A complete full-physics
cooling runtime is unmeasured.

The selected late atmosphere grid contains **60** solved mixture columns over
**6000–6500 K** and **log g = 4.20–5.55**. It passes source, independent depth,
native and matched stellar checks. A **419,300-year** comparison of two inferred
hydrogen references changes temperature by **0.09968 K** and nuclear heat by
**0.0003911%**, with identical final convection. A shared **380,700-year** prefix
is reused. This tests local sensitivity after
the flash, not whether the atmosphere permits ignition in the first place.
The absolute hydrogen reference below the published gravity range remains inferred.

The paper compares Ember with MESA and with Laughlin, Bodenheimer and Adams
(1997), including the requested opacity-map comparison. The separate
[Fortran reconstruction](https://github.com/oklo/Henyey) seeks to reproduce
LBA97's assumptions and methods. A new [late cooling figure](docs/reports/2026-09-20/late_cooling.pdf) shows the checked temperature turn and luminosity decline, with its small data table and plotting script included.
The [standalone HR diagram of both flash histories](docs/figures/2026-09-20/ember_two_sequences_hr.png) shows the alternatives through 38.79 thousand years and 4.624 million years.
Standalone flash diagrams remain separate from the draft.

The public source does not yet contain all physics and numerical updates used
for these continuations. The paper and README describe the calculations ahead
of that source snapshot.

## Aim

Follow an initially **0.1-solar-mass** star through hydrogen exhaustion and helium
white-dwarf cooling to **100 K**, then examine later evolution under explicit
environmental and nucleon-decay assumptions. The complete cooling trajectory
has not yet been calculated.

A matched **1000-year** atmosphere comparison near **log g = 5.42** changes the final effective temperature by **5.918 K** and integrated nuclear heat by **0.00005299%**, with identical final convective regions. Both hydrogen references remain inferred; ignition with a different absolute atmosphere is untested.

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
