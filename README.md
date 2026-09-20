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
are measured from the model preceding shell convection.

| Calculation | Time | Effective temperature | Nuclear power |
| --- | ---: | ---: | ---: |
| Flash resolved from its onset, 1899 mass points | 14.77 thousand yr | 4765 K | 4.177e34 erg/s |
| Later shell burning, separate history, 2115 mass points | 2.254 million yr | 6301 K | 3.143e31 erg/s |

**These calculations do not yet form one verified continuous history.**
The onset calculation passes a local nuclear-power maximum of **3.719e35 erg/s**,
followed by an **81.64%** decline and renewed burning. Over its final
**707.6 years**, power rises by **16.54%**, releasing **8.515e44 erg**. A matched **100-year** comparison
on **1763/1899 mass points** changes nuclear heat by **4.026%** and endpoint
power by **5.589%**. The full pulse's spatial accuracy, largest peak and final
shutdown remain unresolved.

In the later calculation, radius and surface luminosity have declined
**61.36%** and **83.53%** from their saved maxima. Temperature falls from
**6399 K** by **97.80 K** over **1.166 million years**, while nuclear burning supplies
**35.66%** of the current surface luminosity. This supports a turn toward
cooling in that calculation; its connection through the flash remains unverified.

The final **178,100 years** of the later calculation use **18.40 CPU minutes**,
and **8.567 minutes of active elapsed time**, including timestep comparisons.
Trial intervals in that completed segment are at most **25,600 years**;
the accepted solution retains the two shorter steps, with unchanged accuracy
and conservation criteria. The main track through **4.002 Tyr** accounts for
**15.13 elapsed hours** and **36.31 CPU-hours**, excluding source construction
and pauses between jobs. A complete full-physics cooling runtime is unmeasured.

The selected late atmosphere grid contains **48** solved mixture columns over
**6000–6500 K** and **log g = 4.20–5.25**. It passes source, independent depth,
native and matched stellar checks. A **102,400-year** comparison of two inferred
hydrogen references changes temperature by **0.07474 K**, surface luminosity
by **0.001740%**, and nuclear heat by **2.291e-5%**, with identical final convection.
Both histories decline in temperature. This tests local sensitivity; the absolute
hydrogen reference below the published gravity range remains inferred.

The paper compares Ember with MESA and with Laughlin, Bodenheimer and Adams
(1997), including the requested opacity-map comparison. The separate
[Fortran reconstruction](https://github.com/oklo/Henyey) seeks to reproduce
LBA97's assumptions and methods. A new [late cooling figure](docs/reports/2026-09-20/late_cooling.pdf) shows the checked temperature turn and luminosity decline, with its small data table and plotting script included.
The standalone flash diagrams remain separate from the draft.

The public source does not yet contain all physics and numerical updates used
for these continuations. The paper and README describe the calculations ahead
of that source snapshot.

## Aim

Follow an initially **0.1-solar-mass** star through hydrogen exhaustion and helium
white-dwarf cooling to **100 K**, then examine later evolution under explicit
environmental and nucleon-decay assumptions. The complete cooling trajectory
has not yet been calculated.

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
