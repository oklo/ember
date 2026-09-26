# Ember

Ember is a one-dimensional stellar evolution code written in C++23. It calculates the structure, composition and energy transport of very low-mass stars.

[Read the working paper](docs/reports/2026-09-26/ember_status_and_future.pdf) ([LaTeX, figures and plotting data](docs/reports/2026-09-26/README.md)).

## Research status — September 26, 2026

The target is one continuous calculation from a Hayashi starting model through hydrogen burning and helium-white-dwarf cooling to **100 K**. Later environmental heating and nucleon decay are conditional extensions. The same code should determine which lower-mass objects sustain hydrogen burning and which cool as brown dwarfs. The complete calculation is not yet available.

The contraction calculation reaches **1 Gyr**, **2765 K**, **0.1264 solar radii** and **8.414e-4 solar luminosities**. Nuclear burning supplies **99.23%** of its light. The sequence has **3881 accepted intervals**, with **23.31 CPU minutes** of recorded model calculation, excluding atmosphere preparation and separate controls. It includes homogeneous convective mixing and initial-deuterium/pp burning. The [new PMS figure](docs/reports/2026-09-26/pms_main_sequence.pdf) shows its trajectory and the transition to sustained hydrogen burning.

The common evolution engine now combines initial deuterium and explicit hydrogen/helium/carbon/nitrogen inventories, finite convective mixing, and material heat transport. A smooth metal-composition EOS extension supports all **512** main-sequence reference zones while preserving the checked lower-metallicity results. A bounded **1-million-year** integration with finite mixing and material heat converges on the actual 0.1-solar-mass structure; one full step and two half steps differ by **3.856e-8** in logarithmic structure variables. This control omits microscopic drift and does not extend the production trajectory. Some EOS coverage, cool transport and consistent atmospheres remain unfinished.

The common driver carries a new Hayashi model through **6.707 Myr**, to **2994 K** and **0.5146 solar radii**, exhausting initial deuterium. Its **433** accepted intervals pass isotope and energy checks, and exact restart replay passes. The retained calculation uses **4.241 CPU minutes**, including three material-loading stages. Homogeneous early mixing is bounded by its travel time and estimated carried heat. The next step stops because the screened transport prescription does not cover the cool envelope; the full continuous production track remains incomplete. [Driver result](docs/results/lifetime_driver_integration_sept26_v1.json), [early approximation](docs/results/pms_initial_transport_bound_sept26_v1.json).

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
