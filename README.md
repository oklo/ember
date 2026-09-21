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
| Finite convective transport from flash onset, 1983 mass points | 128.3 thousand yr | 5594 K | 2.393e32 erg/s |
| Later shell burning, alternative history, 2115 mass points | 17.91 million yr | 5183 K | 4.608e30 erg/s |

**These are alternative histories, not consecutive pieces of one track.**
The later calculation includes different mixing treatments during its flash;
its subsequent accuracy checks do not validate that entire earlier evolution.
Both calculations use the supplied hydrogen atmosphere as their reference.
Whether the flash occurs under a different atmosphere prescription remains
untested. The flash appears about **0.7546 Gyr** after the atmosphere change;
that delay does not establish independence from the envelope and fuel history.

The onset calculation passes a local nuclear-power maximum of **3.719e35 erg/s**,
followed by an **81.64%** decline and renewed burning. Nuclear power continues to decline through independently checked intervals. The largest single-cell power contribution
is **2.767%**. A matched **100-year** comparison on **1927/1983 mass points**
changes nuclear heat by **8.261%** and final power by **15.24%**. Inserting
mass points produces an initial mechanical readjustment with an integrated
energy residual of **1.844e43 erg** in its first year. A controlled correction
of the initial radii leaves one-year nuclear heat and final power unchanged
within **1e-12** relative; the initial energy change remains explicit in the
budget. The full pulse's spatial accuracy, integrated energy and final
shutdown remain unresolved.

In the later calculation, radius and surface luminosity have declined
**79.42%** and **97.86%** from their saved maxima. Temperature falls from
**6399 K** by **1216 K** over **16.82 million years**, while nuclear burning
supplies **40.25%** of the current surface luminosity. This supports a turn
toward cooling within that calculation.

Completed native calls over the final **3.740 million years** account for **159.3 CPU minutes**;
active elapsed time is **57.01 minutes**, including timestep comparisons. Trial intervals are
at most **51,200 years**; the accepted solution retains the two shorter steps,
with unchanged accuracy and conservation criteria. The main track through
**4.002 Tyr** accounts for **15.13 elapsed hours** and **36.31 CPU-hours**,
excluding source construction and pauses between jobs. A complete full-physics
cooling runtime is unmeasured.

The cooler atmosphere table contains **27** source columns over **5000–6000 K**
and **log g = 5.40–5.80**. Independent **5250 K** and **5750 K** columns test
interpolation between temperature and gravity nodes. At **5250 K**, the
composition response agrees to **0.000007587%** in temperature and **0.1582%**
in pressure.

A matched **2.1-million-year** cooling comparison over **5500–6000 K** changes
endpoint temperature by **2.292 K** and nuclear heat by **0.002680%** between
two hydrogen reference prescriptions. The detailed-onset comparison over
**10,000 years** changes the endpoint by **4.187 K** and nuclear heat by
**0.009269%**. Both comparisons retain identical final convective regions.
These local tests do not establish the absolute hydrogen boundary condition.

A further **10,000-year** comparison entering the **5500–6000 K** atmosphere
interval changes the detailed history's endpoint by **9.515 K**, surface
luminosity by **0.6680%**, and nuclear heat by **0.001667%**, with identical
final convective regions. This measures local sensitivity to the inferred
hydrogen reference; it does not validate its absolute normalization.


A further matched **2.000-million-year** comparison in the **5000–6000 K**
grid changes the endpoint by **6.020 K**, surface luminosity by **0.1435%**
and integrated nuclear heat by **0.009124%**. The surface convective boundary
differs by one mass interface. This tests local atmosphere sensitivity.


The paper compares Ember with MESA and with Laughlin, Bodenheimer and Adams
(1997), including the requested opacity-map comparison. The separate
[Fortran reconstruction](https://github.com/oklo/Henyey) seeks to reproduce
LBA97's assumptions and methods. A new [late cooling figure](docs/reports/2026-09-20/late_cooling.pdf) shows the checked temperature turn and luminosity decline, with its small data table and plotting script included.
The [standalone HR diagram of both flash histories](docs/figures/2026-09-20/ember_two_sequences_hr.png) shows the alternatives through 128.3 thousand years and 17.91 million years.
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
