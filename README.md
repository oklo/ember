# ember

A stellar evolution code for the lowest-mass stars — the ones that burn
hydrogen for trillions of years, turn blue instead of red when their fuel runs
out, and end as helium white dwarfs that outlive everything else in the galaxy.

`ember` is a C++23 development of the FORTRAN stellar-evolution line
reconstructed in [oklo/Henyey](https://github.com/oklo/Henyey), based on
Laughlin, Bodenheimer & Adams (1997). It is an active research code: numerical
checks, comparisons with the earlier model, and physical limitations are
reported alongside each result.

## Status — September 10, 2026

The **0.1 solar-mass model has reached 3.685 trillion years**, with a stable
central region and an outer convective envelope. Hydrogen burning continues.
**Hydrogen exhaustion and white-dwarf cooling have not been reached.** The
latest forward trial stopped at the atmosphere table's 3400 K upper limit.
A fresh trajectory with the checked 3600 K extension is now running toward
3.800 trillion years. That target is not yet a completed result.

| Latest accepted state | Value |
|---|---:|
| Initial composition | XH = 0.7, He3 = 0, Z = 0.02; GS98 metals |
| Mass mesh | 512 points; fixed baryonic mass of 0.1 Msun |
| Central hydrogen mass fraction | 0.03811 |
| Radius | 0.1405 Rsun |
| Luminosity | 0.002375 Lsun |
| Effective temperature | 3400 K |
| Central temperature | 9.290 million K |
| Convective mass fraction | 0.3251 |

The goal remains an evolved helium remnant at **100 K effective temperature**,
using material properties, heat transport and atmospheres valid along its
trajectory. The [current development plan](docs/COLD_REMNANT.md) describes
remaining physics; the [handoff](HANDOFF.md) records calculations and inputs.

The [September 10 working paper](docs/reports/2026-09-10/ember_status_and_future.pdf)
([source and figures](docs/reports/2026-09-10/README.md)) currently presents the
completed comparison through **3.560 trillion years**. It compares Ember
primarily with Laughlin, Bodenheimer and Adams (1997), using stated values and
approximate curves read from their figures. The reconstructed F77 results are
a separate, faint supplementary comparison.

The new hot-core opacity family extends to **zero hydrogen**. Twelve independent
composition checks differ by at most **0.05257%** on tested hot-profile states;
this is an opacity interpolation check, not a stellar lifetime error bound.
All previous hot data and the cool opacity inputs remain unchanged. The selected
EOS has 64 composition tables. The atmosphere extension contains 132 source
models; four low-gravity 3600 K points lie below the EOS density range, and
runtime guards reject those unsupported states. Nearby states at the evolving
star's gravity are supported. See the [opacity](docs/results/tops_exhaustion_hot_profile_all_v2.json)
and [atmosphere](docs/results/nongrey_t3600_extension_v1.json) checks.

Both 512- and 1024-point forward trials reached 3400 K, at 3.685 and 3.683
trillion years. Earlier fresh accuracy controls completed 3.600 trillion years:
a fourfold tighter time control changes central hydrogen by 0.1849%; doubling
the mesh changes it by 1.683%. These checks support continued calculation while
showing that endpoint ages still need refinement. See the
[accuracy comparison](docs/results/evolution_accuracy_3600gyr_v1.json) and
[latest stopped trials](docs/results/evolution_atmosphere_limit_3685gyr_v1.json).

**Production timing:** fresh runs through 3.600 trillion years took **39.78
minutes at 512 points** with the tighter time control and **59.26 minutes at
1024 points** with the original time control, excluding computer suspension.
Atmosphere source generation is separate, reusable work. These are partial
trajectories; the complete cooling-track cost remains unmeasured. Further reuse
of complete nuclear calculations now reduces CPU time by **42.21% in a short
benchmark**, with identical full outputs; a longer comparison is running. The
quoted 3.600-trillion-year timings precede that change. See
[the timing notes](docs/COMPUTATIONAL_COST.md).

New [grain-opacity calculations](docs/GRAINS.md) provide separate absorption
and scattering over wavelength, checked against an independent program and
small-particle analytic limits. A first fixed-atmosphere control combines them
with the computed alumina abundance. Coupling to the atmosphere, material
coverage, grain phase and size, settling and condensate heat transport remain
under development; the evolving star still uses gas-only atmospheres.

The [LBA97 comparison notes](docs/F77_LBA97_COMPARISON.md) distinguish the
published model from the current F77 reconstruction. The latter reaches its
central transition at a substantially different hydrogen abundance. The
[timing notes](docs/COMPUTATIONAL_COST.md) explain the computational costs.

**Verification:** 32 CTest suites pass with installed local inputs. Earlier
same-physics convergence checks at two trillion years found luminosity changes
of +0.06163% when doubling the mesh and −0.00361% with fourfold tighter time
tolerances. Those checks do not establish convergence at the later structural
transition or of a complete lifetime. [Native restarts](docs/RESTART.md) enforce
identical executable bytes, tables and physical selections.

Outstanding work includes zero-hydrogen and condensate-consistent atmospheres,
microscopic diffusion and settling, shell resolution, near-exhaustion nuclear
branches, full thermal-neutrino losses, and dense helium liquid/solid
thermodynamics and cooling boundaries. The [cold-electron and dense-EOS
audits](docs/DENSE_EOS.md) are preparation for those stages. No extrapolated
cooling law is counted as an evolved remnant.

Earlier controls and numerical comparisons are documented in
[FORWARD_EVOLUTION.md](docs/FORWARD_EVOLUTION.md),
[EXTENDED_EVOLUTION.md](docs/EXTENDED_EVOLUTION.md), and
[NONGREY.md](docs/NONGREY.md). The wider mass–composition survey is described in
[LIFETIME_SURVEY.md](docs/LIFETIME_SURVEY.md).

## Design

**Analytic derivatives throughout zone assembly.** Physics modules provide
their own partials, and the structure equations propagate them by the chain
rule. Assembly uses one EOS response per endpoint; finite differences remain
an independent check. An EOS without the required response derivatives must
report that explicitly.

**Logarithmic variables.** The solver works in `(ln r, ln rho, ln T, L)`.
Density spans eighteen decades between a giant's photosphere and a white
dwarf's core; linear variables cannot be conditioned across that.

**Physics behind interfaces.** `Eos`, `Opacity`, `Nuclear`, `Atmosphere` are
abstract. Swapping a table is a constructor argument, and two implementations
can be run against each other on the same track — which is how you find out
whether a result is physics or an artifact of one table.

**Tests that are physics statements.** Each equation-of-state test names a
limit where the answer is known independently of this code: the ideal gas, the
radiation-dominated limit, the Chandrasekhar constant, the Maxwell relation
behind `grad_ad`. A regression then reads as a broken law, not a changed
number.

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

## Physics sources

| Ingredient | Source |
|---|---|
| Composition | GS98 elemental inventory for forward evolution; legacy AAG21 carrier abundances retained |
| Low-T opacity | AESOPUS 2.1 gas (Marigo et al. 2024), nine X planes at Z=.01/.02/.03, log R to 6; Ferguson (2005) also available with grains |
| High-T opacity | LANL TOPS ATOMIC with refined hydrogen-poor composition planes, Z=.01/.02/.03; source-domain masks and smooth blends |
| Equation of state | FreeEOS 3.0 EOS1 represented by one C2 Helmholtz potential; 64-plane GS98 H/He3 family through XH=0 within its original thermal domain; isotope corrections and masked invalid states |
| Conduction | Ioffe pure-ion tables; Cassisi et al. (2007, 2021), Blouin et al. (2020); approximate mixture resistivities |
| Convection | Böhm–Vitense MLT; Schwarzschild or full-EOS Ledoux buoyancy; semiconvective and thermohaline composition transport |
| Nuclear rates | Reduced pp network, SFII S-factor quadrature, finite-degeneracy SVH screening, atomic mass-defect heating; legacy prescription retained |
| Atmosphere | TLUSTY208/SYNSPEC54 non-grey gas grid at tau=100; 96 nodes selected, 108-node candidate under refinement; AMES-COND/grey default control retained |
| Thermal neutrinos | Optional HRW plasma-decay sink with analytic derivatives; other thermal channels remain to be supplied |

## Licence

MIT. The opacity, equation-of-state and conductivity data are redistributed under their
own terms; see `data/*/README.md`.
