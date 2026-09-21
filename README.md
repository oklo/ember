# Ember

Ember is a one-dimensional stellar evolution code written in C++23. It follows
the structure, composition and energy transport of very low-mass stars over
trillions of years.

[Read the working paper](docs/reports/2026-09-21/ember_status_and_future.pdf)
([LaTeX and figure files](docs/reports/2026-09-21/README.md)).

## Research status

The **0.1-solar-mass** calculation with microscopic diffusion and metal settling
reaches **4.002 trillion years**, **4612 K**, and central hydrogen **X = 3.528e-8**.
Hydrogen burning continues outside the depleted core. A continuation with a
white-dwarf atmosphere develops a helium-3 shell flash.

Two calculations examine the flash and subsequent evolution. Times below
are measured from their common model preceding shell convection.

| Calculation | Time | Effective temperature | Nuclear power |
| --- | ---: | ---: | ---: |
| Finite convective transport from flash onset, 1983 mass points | 4.619 million yr | 6025 K | 1.371e31 erg/s |
| Later shell burning, alternative history, 2115 mass points | 42.48 million yr | 4771 K | 3.166e30 erg/s |

**These are alternative histories, not consecutive pieces of one track.**
The later calculation includes different mixing treatments during its flash;
its subsequent accuracy checks do not validate that entire earlier evolution.
Both calculations use the supplied hydrogen atmosphere as their reference.
Whether the flash occurs under a different atmosphere prescription remains
untested. The flash appears about **0.7546 Gyr** after the atmosphere change;
that delay does not establish independence from the envelope and fuel history.

The onset calculation passes a local nuclear-power maximum of **3.719e35 erg/s**,
followed by an **81.64%** decline and renewed burning. Nuclear power continues to decline through independently checked intervals. The largest single-cell power contribution
is **1.501%**. A matched **100-year** comparison on **1927/1983 mass points**
changes nuclear heat by **8.261%** and final power by **15.24%**. Inserting
mass points produces an initial mechanical readjustment with an integrated
energy residual of **1.844e43 erg** in its first year. A controlled correction
of the initial radii leaves one-year nuclear heat and final power unchanged
within **1e-12** relative; the initial energy change remains explicit in the
budget. The full pulse's spatial accuracy, integrated energy and final
shutdown remain unresolved.

In the later calculation, radius and surface luminosity have declined
**82.91%** and **98.94%** from their saved maxima. Temperature falls from
**6399 K** by **1628 K** over **41.39 million years**, while nuclear burning
supplies **55.85%** of the current surface luminosity. This supports a turn
toward cooling within that calculation.

The **24.48–29.48 million year** interval uses **74.91 CPU minutes** and **40.23 wall minutes**, including timestep comparisons and rejected trials. Nuclear power falls **7.869%**, while the whole-star helium-3 inventory grows **0.2052%**. At **29.48 million years**, helium-3 production is **1.263e12 g/s** and destruction is **4.106e11 g/s**. Helium-3 fusion supplies **23.12%** of nuclear power, down from **26.97%** five million years before. The discrete source sum reproduces the inventory change to better than **1e-10** relative. This establishes the balance of production and destruction; it does not establish stability against a later flash. The full cooling-track runtime remains unmeasured. Both histories now turn toward lower effective temperature. The finite-mixing calculation reaches **6422 K** near **1.139 million years** before declining.

A provisional pre-main-sequence calculation now follows **100,000 years** of contraction with its initial deuterium and helium-3 inventories intact. It uses **74.79 CPU seconds**, excluding atmosphere preparation, and ends at **2876 K** with radius **1.389 solar radii**. An independent atmosphere inside the starting cell differs from interpolation by **−1.304%** in matching temperature and **+2.156%** in pressure. A local stellar sensitivity test changes the initial effective temperature by **10.89 K** and radius by **0.5213%**. The tested trace-isotope effect is much smaller. Local grid refinement remains in progress; this is not yet the accepted physical reference track.
 Trial intervals are
at most **819,200 years**; the accepted solution retains the refined substeps,
with unchanged accuracy and conservation criteria. The main track through
**4.002 Tyr** accounts for **15.13 elapsed hours** and **36.31 CPU-hours**,
excluding source construction and pauses between jobs. A complete full-physics
cooling runtime is unmeasured.

The selected cooling atmosphere window has **18** columns over **4750–5500 K** and **log g = 5.80–6.00** at three hydrogen fractions. Twelve occupied rows are preserved exactly; six new columns pass independent source and depth checks. The combined table passes interpolation, EOS and full stellar-response comparisons, followed by **6 million years** of reviewed cooling. Higher-gravity columns through **log g = 6.20** pass source and response checks; the missing **4500 K** coverage currently limits further cooling. A matched **million-year** reference comparison gives **0.2632 K** and **0.0002222%** in nuclear heat, with the same final convection. Independent **5250 K** and **5750 K** columns test interpolation within the broader source collection; at **5250 K**, the composition response agrees to **0.000007587%** in temperature and **0.1582%** in pressure.

A matched **2.1-million-year** cooling comparison over **5500–6000 K** changes
endpoint temperature by **2.292 K** and nuclear heat by **0.002680%** between
two hydrogen reference prescriptions. The detailed-onset comparison over
**10,000 years** changes the endpoint by **4.187 K** and nuclear heat by
**0.009269%**. Both comparisons retain identical final convective regions.
These local tests do not establish the absolute hydrogen boundary condition.

A matched **10,000-year** numerical comparison reduces native CPU time from
**46.28 to 27.65 minutes**, a **40.25%** saving, with unchanged accuracy
criteria. Temperature differs by **0.07287 K**, nuclear heat by **0.03629%**,
and the convective regions agree. A separate **640-year** timestep comparison
passes. These measurements support the selected detailed continuation;
they do not establish the runtime or accuracy of a complete cooling track.

A matched **7500-year** comparison entering the **6000–6500 K** atmosphere interval gives differences of **1.008 K**, **0.06316%** in surface luminosity and **0.00002583%** in nuclear heat, with identical final convection. Each accepted interval passes the original full-step/two-half-step checks.

A further **25,000-year** comparison below **log g = 4.50** gives **0.4703 K** difference in endpoint temperature and **0.006257%** in nuclear heat. The final convective regions agree. The selected extension uses the measured gravity response of solved near-hydrogen atmospheres, preserving the occupied reference values. This checks local sensitivity, not the absolute atmosphere.

Microscopic species transport extends to **1.000 MK** under a locally checked ionization approximation. A **1600/800/400-year** comparison gives **0.2686%** difference in nuclear heat with unchanged accuracy and conservation criteria. The heat-transport blend remains at **2–3 MK**.

A matched **25,000-year** comparison of the corrected hot reference tables gives **1.366 K**, **0.06345%** in surface luminosity and **0.01279%** in nuclear heat, with identical final convection. Adaptive full intervals may reach **2560 years**; every accepted interval retains the original checks.

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
LBA97's assumptions and methods. A new [late cooling figure](docs/reports/2026-09-21/late_cooling.pdf) shows the checked temperature turn and luminosity decline, with its small data table and plotting script included.
The [standalone HR diagram of both flash histories](docs/figures/2026-09-21/ember_two_sequences_hr.png) shows the alternatives through 4.619 million years and 42.48 million years.
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
