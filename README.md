# Ember

Ember is a one-dimensional stellar evolution code written in C++23. It calculates the structure, composition and energy transport of very low-mass stars.

[Read the working paper](docs/reports/2026-09-27/ember_status_and_future.pdf) ([LaTeX, figures and plotting data](docs/reports/2026-09-27/README.md)).

## Research status — September 27, 2026

The objective is one continuous calculation from a Hayashi starting model,
through hydrogen burning and helium-white-dwarf cooling, to extreme cold and
conditional disappearance through nucleon decay. The same program should
resolve the stellar/brown-dwarf boundary as initial mass and composition vary.
The full calculation is not yet established.

The common `ember-evolve --lifetime` program has passed **3.02 Tyr** from a
Hayashi start with hot transport selected and no rejected intervals. A completed
comparison with a stricter fixed-metal allowance reaches **2.102 Tyr**, with **5931 accepted intervals** in
**66.16 CPU minutes / 35.41 elapsed minutes**, including loading but excluding
source-table preparation. It remains fully convective at **3020 K**, with
hydrogen **X = 0.4170**. All accepted isotope and energy checks pass.
An overly restrictive atmosphere metal allowance caused **28 rejected trials**
and a final rounding-sensitive audit failure after repeated timestep reductions.
[Retained calculation](docs/results/continuous_lifetime_sept26_v2.json).

The common driver selects hot microscopic H/He/metal transport, conservative
mixing heat, metal-dependent atmospheres and the existing composition-dependent
interior opacities. A saved radiative-core model evolves with nonzero settling,
passes the physical conservation checks, and restarts exactly. Initial
fully mixed deuterium burning uses tighter composition convergence than the
stratified diffusion solve; the physical audits stay unchanged. **All 76 tests
pass** in the working and isolated builds with the required datasets.
[Common driver](docs/LIFETIME_DRIVER.md),
[diffusion and restart test](docs/results/screened_radiative_core_sept27_v1.json),
[atmosphere coverage](docs/results/atmosphere_lifetime_atlas_sept27_v1.json),
[opacity coverage](docs/results/opacity_lifetime_extension_sept27_v1.json).

A fresh Hayashi calculation now selects the restored hydrogen- and
metal-dependent atmospheres from the start. All 5932 saved early surface
queries are unchanged; the late-source tests cover the settling comparison
through 3.749 Tyr. Trace-helium and WD atmosphere coverage remains unfinished.
[Hydrogen-rich atmosphere recovery](docs/results/atmosphere_hydrogen_recovery_sept27_v1.json). New high-gravity
hydrogen columns use the same joining depth as the main-sequence atmospheres.
The microscopic metal law assumes complete ionization; unsupported cool
radiative layers are rejected. Initial deuterium uses a checked whole-star
mixing approximation. Nuclear burning includes the pp chain and explicit
C12/C13/N14 conversion with Solar Fusion III rates, rather than a full CNO network.

Separate late-evolution calculations develop a helium-3 shell pulse after an
atmosphere adjustment and then turn toward cooling. They establish a possible
instability under that boundary treatment, not atmosphere-independent ignition.
The continuous calculation tests whether the pulse survives a consistent history.
[Working paper](docs/reports/2026-09-27/ember_status_and_future.pdf),
[Fortran comparison](docs/F77_LBA97_COMPARISON.md).

The working paper is now 30 pages, reduced from 63 while retaining the
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
