# Ember

Ember is a one-dimensional stellar evolution code written in C++23. It calculates the structure, composition and energy transport of very low-mass stars.

[Read the working paper](docs/reports/2026-09-27/ember_status_and_future.pdf) ([LaTeX, figures and plotting data](docs/reports/2026-09-27/README.md)).

## Research status — September 27, 2026

The objective is one continuous calculation from a Hayashi starting model,
through hydrogen burning and helium-white-dwarf cooling, to extreme cold and
conditional disappearance through nucleon decay. The same program should
resolve the stellar/brown-dwarf boundary as initial mass and composition vary.
The full calculation is not yet established.

The common `ember-evolve --lifetime` program has reached **3.213 Tyr** from a
Hayashi start with composition-dependent atmospheres: **6099 accepted intervals,
none rejected**, at **3154 K** and hydrogen **X = 0.2385**. It remains fully
convective. The retained sequence has used **5.182 CPU hours**, including
contraction, step-doubling comparisons and startup, with separate controls excluded.
[Current calculation](docs/results/continuous_lifetime_sept27_v4.json).
A fixed-metal-atmosphere comparison reaches **3.078 Tyr**, also without rejected
intervals. [Comparison](docs/results/continuous_lifetime_sept26_v3.json).

The common driver selects hot microscopic H/He/metal transport, conservative
mixing heat, metal-dependent atmospheres and the existing composition-dependent
interior opacities. A saved radiative-core model evolves with nonzero settling,
passes the physical conservation checks, and restarts exactly. Initial
fully mixed deuterium burning uses tighter composition convergence than the
stratified diffusion solve; the physical audits stay unchanged. **All 77 tests
pass** in the isolated build with the required datasets.
[Common driver](docs/LIFETIME_DRIVER.md),
[diffusion and restart test](docs/results/screened_radiative_core_sept27_v1.json),
[atmosphere coverage](docs/results/atmosphere_lifetime_atlas_sept27_v1.json),
[opacity coverage](docs/results/opacity_lifetime_extension_sept27_v1.json).

The continuous calculation retains optical depth 100 through its composition-dependent
atmosphere family. Adding later table coverage leaves the accepted stellar state
unchanged and reproduces an independent ten-interval restart exactly. All 5932 tested early
surface values remain unchanged. The combined inputs cover all 67 saved
settling comparison states at the original optical depth of 100. Cool,
high-gravity rows include bounded estimates, which still need source checks;
continuous WD cooling is not yet established.
[Current calculation](docs/results/continuous_lifetime_sept27_v4.json),
[atmosphere coverage and restart evidence](docs/results/hydrogen_envelope_integration_sept27_v1.json).
Independent atmosphere columns now pass at 4800 K/log g = 6.2 and
4700 K/log g = 6.0, 6.1 and 6.2 with unchanged acceptance criteria. A density-inversion
correction also removes rejection caused solely by a remote unstable EOS
endpoint; the finite-mixing treatment of the cool envelope remains under test.
[Atmosphere checks](docs/results/highg_energy_balance_sept27_v1.json),
[density-inversion checks](docs/results/local_density_inversion_sept27_v1.json),
[adjacent atmosphere columns](docs/results/highg_adjacent_columns_sept27_v1.json).
The selected EOS supports the nodes and cell boundaries of all 67 saved settling
structures, including the composition derivatives used for diffusion and heat.
[EOS coverage](docs/results/eos_lifetime_coverage_sept27_v1.json).
A 50-Gyr timestep comparison used 25 intervals instead of 73, with a
0.0004605% luminosity difference and 0.001502% helium-3-reservoir difference.
The continuous star now uses that tested time-error setting; structure accuracy
and all conservation checks are unchanged.
[Timestep comparison](docs/results/lifetime_time_accuracy_sept27_v1.json).
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
