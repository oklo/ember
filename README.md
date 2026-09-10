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

The completed **0.1 solar-mass model reaches 3.600 trillion years** from a
specified static main-sequence model. It now has a central region stable against convection and an outer convective
envelope. Hydrogen burning continues. **Hydrogen exhaustion and white-dwarf cooling have not been reached.**
The current goal is to evolve the helium remnant to **100 K effective temperature**
using thermodynamics, heat transport and atmospheres valid in the conditions
the star reaches.

| Latest completed model | Value |
|---|---:|
| Initial composition | XH = 0.7, He3 = 0, Z = 0.02; GS98 metals |
| Mass mesh | 512 points; fixed baryonic mass of 0.1 Msun |
| Central hydrogen mass fraction | 0.09839 |
| Radius | 0.1419 Rsun |
| Luminosity | 0.002085 Lsun |
| Effective temperature | 3274 K |
| Central temperature | 8.670 million K |
| Convective mass fraction | 0.5073 |

The [September 10 working paper](docs/reports/2026-09-10/ember_status_and_future.pdf)
([source, figures and reproduction instructions](docs/reports/2026-09-10/README.md))
compares Ember with the published Laughlin, Bodenheimer and Adams (1997)
model, using values stated in the paper and approximate curves read from its
figures. The current F77 reconstruction is a separate supplementary comparison.
The paper also describes computing time and the physics still needed for the
remnant. The
[current development plan](docs/COLD_REMNANT.md) and latest
[handoff entry](HANDOFF.md) record subsequent work and job status.

The completed run uses **64 composition tables for the equation of state**,
refined hydrogen-poor TOPS opacity, and **120 atmosphere models** that resolve
absorption and scattering over frequency. Convection accounts for composition
gradients; plasma-neutrino losses enter the energy equation. The new 512-point calculation
contains 2354 states from the initial model through 3.560 trillion years,
followed by 277 accepted steps to 3.600 trillion years.

A refined table of **120 atmospheres extends to hydrogen fraction 0.1** at
warm surface temperatures, including new points at hydrogen fraction 0.15.
At an independent atmosphere near that composition, the matching-temperature
discrepancy falls from 2.134% to 0.6468%. Independent checks at hydrogen
fractions 0.125 and 0.175 differ by 0.7971% and 1.089%. These are local
interpolation checks, not a bound on lifetime accuracy.

At 3.560 trillion years, both the 512-point and 1024-point calculations have stable central regions,
containing 37.23% and 37.77% of their mass, respectively, and extending to about
45% of their radius. At the same age their luminosities differ by 0.07813%,
but central hydrogen differs by 1.752%. The earlier stable shell has reached
the center. Further mesh and time-step checks are needed to determine the
transition age accurately. The earlier checkpoint counter limit was corrected
without loosening accuracy tolerances. See the
[completed comparison](docs/results/evolution_transition_3560gyr_v1.json).

The next continuation reached 3.609 trillion years before encountering the
interior opacity family's lower hydrogen boundary. Additional TOPS source
calculations are being requested; no table is extrapolated to continue the run.
Separate calculations with fourfold tighter time tolerances and with 1024 mesh
points are checking the transition, using the new two-worker mode.

**Performance:** `--step-workers 2` runs independent time-step estimates on two
CPU workers, and conduction calculations now reuse repeated source lookups.
A checked 512-point benchmark uses 32.76% less elapsed time, with complete
output files unchanged byte for byte. The default remains one worker. See
[the measurements and CPU/GPU assessment](docs/COMPUTATIONAL_COST.md).

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
