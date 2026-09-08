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

The physical modules and damped Henyey solver now converge an experimental
**0.1 Msun static equilibrium** with molecular/partially ionized CMS19 H/He
thermodynamics, AESOPUS/TOPS radiative opacity, pp heating, mixing-length
convection and a grey atmosphere. Fine meshes settle near **R=.1286 Rsun,
L=.000976 Lsun, Teff=2845 K**, with energy and independent virial checks.
Thirteen test suites pass, including the independent Lane–Emden benchmark.

This is not yet a calibrated physical stellar model. The source EOS has
substantial local pressure/entropy consistency errors and an internal-energy
join defect. CMS19 is therefore **restricted to static calculations**;
positive-dt energy equations reject it. Physical atmosphere data, consistent
caloric/composition physics, mixing and time-step control remain pending.
No evolutionary run has been made. See [docs/EQUILIBRIUM.md](docs/EQUILIBRIUM.md)
for reproduction, numerical checks and limitations, and
[docs/CMS19.md](docs/CMS19.md) for the EOS audit.

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
| High-T opacity | LANL TOPS ATOMIC at X=.7, Z=.02; OPAL GS98 subset also available; strict smooth blends |
| Equation of state | CMS19 pressure/entropy, experimental static only; analytic ideal+FD also available |
| Conduction | Cassisi et al. (2007) *(pending)* |
| Convection | Böhm–Vitense MLT, optically thick; Schwarzschild criterion |
| Nuclear rates | JINA REACLIB / NACRE II *(pending)* |
| Atmosphere | Eddington grey implemented; table reader ready, physical grids pending |

## Licence

MIT. The opacity and equation-of-state tables are redistributed under their
own terms; see `data/*/README.md`.
