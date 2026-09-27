# Ember

Ember is a one-dimensional stellar evolution code written in C++23. It calculates the structure, composition and energy transport of very low-mass stars.

[Read the working paper](docs/reports/2026-09-26/ember_status_and_future.pdf) ([LaTeX, figures and plotting data](docs/reports/2026-09-26/README.md)).

## Research status — September 26, 2026

The objective is one continuous calculation from a Hayashi starting model,
through hydrogen burning and helium-white-dwarf cooling, to extreme cold and
conditional disappearance through nucleon decay. The same program should
resolve the stellar/brown-dwarf boundary as initial mass and composition vary.
The full calculation is not yet established.

The common `ember-evolve --lifetime` program follows the same initially
0.1-solar-mass star to **2.102 Tyr**, with **5931 accepted intervals** in
**66.16 CPU minutes / 35.41 elapsed minutes**, including loading but excluding
source-table preparation. It remains fully convective at **3020 K**, with
hydrogen **X = 0.4170**. All accepted isotope and energy checks pass.
An overly restrictive atmosphere metal allowance caused **28 rejected trials**
and a final rounding-sensitive audit failure after repeated timestep reductions.
[Retained calculation](docs/results/continuous_lifetime_sept26_v2.json).

The driver now selects the existing hot microscopic H/He/metal transport with
conservative mixing heat in the cool convective envelope. A fresh Hayashi run
uses this prescription and a broader, measured atmosphere metal allowance.
Both changes are explicit configuration choices. Existing broad main-sequence
atmosphere tables join smoothly to the contraction boundary. Source masks and
conservation acceptance criteria remain unchanged. **All 74 tests pass** in an
isolated build with the required datasets.
[Common driver](docs/LIFETIME_DRIVER.md),
[transport checks](docs/results/envelope_transport_integration_sept26_v1.json),
[atmosphere comparisons](docs/results/archived_ms_atmosphere_check_sept26_v1.json).

Continuous boundary selection for substantial settling and WD cooling remains
unfinished. The microscopic metal law assumes complete ionization; unsupported
cool radiative layers are rejected. Initial deuterium is handled in a checked
whole-star mixing approximation. Nuclear burning includes the pp chain and
explicit C12/C13/N14 conversion, with Solar Fusion III rates in the common driver;
this is not a full CNO network.

Separate late-evolution calculations develop a helium-3 shell pulse after an
atmosphere adjustment and then turn toward cooling. They establish a possible
instability under that boundary treatment, not atmosphere-independent ignition.
The continuous calculation tests whether the pulse survives a consistent history.
[Working paper](docs/reports/2026-09-26/ember_status_and_future.pdf),
[Fortran comparison](docs/F77_LBA97_COMPARISON.md).

The coherent source update is prepared locally. Public runtime/test-data
availability and the shorter paper are still being completed.

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
