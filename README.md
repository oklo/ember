# Ember

Ember is a one-dimensional stellar evolution code written in C++23. It calculates the structure, composition and energy transport of very low-mass stars.

[Read the working paper](docs/reports/2026-09-26/ember_status_and_future.pdf) ([LaTeX, figures and plotting data](docs/reports/2026-09-26/README.md)).

## Research status — September 26, 2026

The target is one continuous calculation from a Hayashi starting model through hydrogen burning and helium-white-dwarf cooling to **100 K**. Later environmental heating and nucleon decay are conditional extensions. The same code should determine which lower-mass objects sustain hydrogen burning and which cool as brown dwarfs. The complete calculation is not yet available.

The common evolution program now follows the same star from a Hayashi start
to **4 Gyr**, with initial D/pp/CN burning, plasma losses and checked
whole-star convective mixing. Composition heat is included after initial D
exhaustion. At 1 Gyr it has **2765 K**, **0.1264 solar radii** and **99.17%**
nuclear support; by 4 Gyr it sustains hydrogen burning. The retained sequence
has **1034 accepted steps**, no rejected steps, and **9.785 CPU minutes**
including material loading, excluding atmosphere preparation and controls.
The [PMS figure](docs/reports/2026-09-26/pms_main_sequence.pdf) shows contraction
through main-sequence arrival. [Transport and conservation assessment](docs/results/convective_transport_integration_sept26_v1.json).

On the quiet main sequence, a matched timestep-ceiling increase reduces
1–4-Gyr evolution from **301** to **15** accepted steps and uses **4.534 times
less CPU**. The surface-temperature difference is **6.936e-5 K**, with all
local accuracy and independent conservation checks unchanged.
[Matched comparison](docs/results/lifetime_step_ceiling_sept26_v1.json).

The well-mixed approximation is checked against convective times, estimated
microscopic separation and stellar sensitivity to uncertain heat transport.
It rejects a radiative species boundary rather than silently suppressing its
flux. Radiative microscopic transport, extended atmosphere composition
coverage and the full continuous lifetime remain unfinished. New gas-opacity
sources at **X = 0.69** and **X = 0.65** are in preparation; they are not yet
accepted atmosphere boundaries. Deuterium exhaustion no longer forces an
unsupported cool diffusion calculation.

Separate long-term calculations include microscopic diffusion and metal settling. The pre-flash calculation reaches **4.002 Tyr**, **4612 K**, and central hydrogen **X = 3.528e-8**. Hydrogen burning continues outside the depleted core. A continuation with a white-dwarf atmosphere develops a helium-3 shell pulse.

| Flash and cooling calculation | Time from common pre-flash model | Effective temperature | Nuclear power |
| --- | ---: | ---: | ---: |
| Finite mixing from onset, 1983 mass points | 5.744 Myr | 5903 K | 1.126e31 erg/s |
| Alternative mixing history, 2115 mass points | 42.48 Myr | 4771 K | 3.166e30 erg/s |

These are alternative histories, not consecutive parts of one trajectory. Both turn toward cooling, but nuclear burning still supplies **55.85%** of the alternative model's light. The pulse's full strength, spatial accuracy and independence from the atmosphere adjustment remain unresolved. The continuous PMS-started calculation will test whether the pulse survives a consistent composition and atmosphere history.

The paper retains the Ember–MESA comparisons, the LBA97 opacity-map comparison, and the lifetime timeline. The [late cooling figure](docs/reports/2026-09-26/late_cooling.pdf) and [standalone HR comparison](docs/figures/2026-09-21/ember_two_sequences_hr.png) show the unchanged late-evolution endpoints. The separate [Fortran reconstruction](https://github.com/oklo/Henyey) seeks to reproduce LBA97's assumptions and methods.

The public source does not yet contain all physics and numerical updates used for these calculations. The paper and README describe results ahead of that source snapshot.

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
