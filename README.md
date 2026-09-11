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

The **0.1 solar-mass model has reached 3.818 trillion years**. Its core
transports heat without convection and is surrounded by a convective envelope.
Hydrogen burning continues. The calculation stops at its atmosphere table's
**4000 K** limit. Its complete history, final structure, input receipt and
transport classification pass the recorded checks.

| Latest accepted state | Value |
|---|---:|
| Initial composition | XH = 0.7, He3 = 0, Z = 0.02; GS98 metals |
| Mass mesh | 512 points; fixed baryonic mass of 0.1 Msun |
| Central hydrogen mass fraction | X = 0.004580 |
| Surface hydrogen mass fraction | X = 0.1731 |
| Radius | 0.1379 Rsun |
| Luminosity | 0.004387 Lsun |
| Effective temperature | 4000 K |
| Central temperature | 11.05 million K |
| Convective mass fraction | 0.1170 |

[History and endpoint checks](docs/results/evolution_atmosphere_limit_3818gyr_v1.json).

The goal remains an evolved helium remnant at **100 K effective temperature**,
using material properties, heat transport and atmospheres valid along its
trajectory. The [current development plan](docs/COLD_REMNANT.md) describes
remaining physics; the [handoff](HANDOFF.md) records calculations and inputs.

The [September 10 working paper](docs/reports/2026-09-10/ember_status_and_future.pdf)
([source and figures](docs/reports/2026-09-10/README.md)) now presents the
track through **3.795 trillion years**. All its Ember curves refer to that
calculation. It compares Ember
primarily with Laughlin, Bodenheimer and Adams (1997), using stated values and
approximate curves read from their figures. The reconstructed F77 results are
a separate, faint supplementary comparison.
The photospheric figure compares the LBA97 and Ember opacity maps on identical
axes and a shared approximate opacity scale, with both 0.1-solar-mass tracks.
It marks the grains missing from Ember's gas-only background.

The ATOMIC opacity family for the hot interior extends to **zero hydrogen**. Twelve independent
composition checks differ by at most **0.05257%** on tested hot-profile states;
this is an opacity interpolation check, not a stellar lifetime error bound.
All previous hot data and the cool opacity inputs remain unchanged. The plotted
calculation uses a separately checked **72-table EOS** and **144-model
atmosphere extension through 3800 K**. All 144 source models have EOS support.
Two independent 3700 K atmospheres differ from interpolation by at most
0.7399% in matching temperature and 0.5707% in gas pressure; these are local
checks, not a lifetime error bound. Depth checks pass, and no condensates appear
in the twelve added structures. See the
[atmosphere acceptance](docs/results/nongrey_t3800_acceptance_v1.json) and
[EOS data notes](data/eos/README.md).

The **156-model extension through 4000 K** also passes its runtime/EOS,
intermediate-temperature, depth and condensation checks and supplies the
latest 3.818-trillion-year calculation. [Acceptance record](docs/results/nongrey_t4000_acceptance_v1.json).

The **172-model extension through 4400 K** is now checked. Its independent
matching-state differences are below **0.5987%**, eight depth comparisons pass,
and no condensates occur in the sixteen new structures. All 156 existing states
remain unchanged. The saved star is continuing with this table, using a tested
restart procedure that preserves its structure and composition.
[Atmosphere acceptance](docs/results/nongrey_t4400_acceptance_v1.json),
[restart comparison](docs/results/atmosphere_extension_restart_v1.json).

Numerical controls through 3.600 trillion years find central-hydrogen changes
of 0.1849% with a fourfold tighter time control and 1.683% with twice as many
mass points. These measurements apply to that age and those fixed inputs.
[Accuracy comparison](docs/results/evolution_accuracy_3600gyr_v1.json).

**Production timing:** the partial trajectory through 3.818 trillion years
used **41.19 CPU minutes** and **25.09 awake elapsed minutes**, with 7147 accepted
steps and 269 rejected attempts. Atmosphere generation is separate, reusable
work. The complete cooling-track cost remains unmeasured. Reuse of complete nuclear calculations reduced CPU time
by **42.21% in a short benchmark**. The completed longer comparison through
3.757 trillion years used **37.31% less CPU time**, with identical entire output
files: **162.4 versus 101.8 CPU minutes**. Background load varied; this is one
long pair, not a timing guarantee. Sharing the electron calculation within each nuclear response
saves a further **14.97%** in a short comparison. Both changes are installed;
both reported stellar calculations include them. See
[the timing notes](docs/COMPUTATIONAL_COST.md) and
[the latest history and checks](docs/reports/2026-09-10/evolution_latest_provenance.json).

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
