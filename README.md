# ember

A stellar evolution code for the lowest-mass stars — the ones that burn
hydrogen for trillions of years, turn blue instead of red when their fuel runs
out, and end as helium white dwarfs that outlive everything else in the galaxy.

`ember` is a ground-up C++23 rewrite of the FORTRAN line developed in
[oklo/Henyey](https://github.com/oklo/Henyey), which reconstructed the code of
Laughlin, Bodenheimer & Adams (1997). That repository is the historical
artifact and stays as it is. This one is where the physics moves forward: the
reconstruction's job was to reproduce 1997, and it did; `ember`'s job is to be
right by today's standards.

## Status

Built for the lowest-mass stars first; the interfaces are shaped so
that the ultra-cold evolution of *massive* remnants - a solar remnant, or
something near the Chandrasekhar mass - can be added as further terms rather
than as a rewrite. See `docs/ROADMAP.md`.

The code now completes a **coupled 10-Gyr evolution experiment for a
0.1 Msun star** with composition-dependent FreeEOS potentials and
AESOPUS/TOPS opacity, explicit pp burning, instantaneous convective mixing,
and an AMES-COND non-grey boundary. At 4096 points the final model has
**R=.128915 Rsun, L=.00091942 Lsun, Teff=2799.32 K**;
hydrogen decreases from .7 to .697402 and He3 grows to .00259806.
Eighteen test suites pass, including source-data, conservation, coupled
thermal evolution, timestep refinement and the independent Lane–Emden benchmark.

This is a bounded early-evolution experiment. Time starts from a specified
static composition with zero He3, rather than formation. The consistent EOS
still treats metals as helium; opacity uses a nominal isotope approximation;
the non-grey atmosphere retains a frozen solar-mixture T/P boundary with
explicit small-composition-change limits. The run does not reach hydrogen
exhaustion or white-dwarf cooling. See [docs/EVOLUTION.md](docs/EVOLUTION.md)
for equations, numerical checks, source coverage and reproduction.
[docs/FREEEOS.md](docs/FREEEOS.md) documents the underlying potential and its
source-fit uncertainties. The retained [CMS19 experiment](docs/CMS19.md)
keeps its static-only energy guard.

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

`-ffast-math` is deliberately *not* used: it licenses the compiler to assume no
NaN or infinity, and a stellar model legitimately probes states where a table
returns one. Those should surface, not be optimised away.

## Physics sources

| Ingredient | Source |
|---|---|
| Composition | Asplund, Amarsi & Grevesse (2021) |
| Low-T opacity | AESOPUS 2.1 gas (Marigo et al. 2024), log R to 6; Ferguson (2005) also available with grains |
| High-T opacity | LANL TOPS ATOMIC at X=.6/.65/.7/.75, Z=.02; OPAL GS98 subset also available; strict smooth blends |
| Equation of state | FreeEOS 3.0 EOS1 data represented by one C2 Helmholtz potential, variable H with He3 number-density mapping; CMS19 static-only and analytic ideal+FD alternatives |
| Conduction | Cassisi et al. (2007) *(pending)* |
| Convection | Böhm–Vitense MLT, optically thick; Schwarzschild criterion |
| Nuclear rates | Reduced pp network with explicit He3, atomic mass-defect heating and classical screening; modern rate/screening audit pending |
| Atmosphere | AMES-COND-2000 non-grey tau=100 states via MESA, explicit solar-mixture proxy; Eddington grey alternative |

## Licence

MIT. The opacity and equation-of-state tables are redistributed under their
own terms; see `data/*/README.md`.
