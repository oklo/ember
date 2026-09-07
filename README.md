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

Early. The equation of state, composition handling and build are in and
tested; the solver and the remaining physics modules are being written. See
`docs/ROADMAP.md`.

## Design

**Analytic derivatives everywhere.** The Fortran ancestor formed its Jacobian
by centred differences — five equation-of-state calls per zone per Newton
iteration. Every physics module here returns its own derivatives, which costs
accuracy nothing and time a factor of several.

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
| Low-T opacity | Ferguson et al. (2005), grains included |
| High-T opacity | OPAL / OPLIB *(pending)* |
| Equation of state | Chabrier, Mazevet & Soubiran (2019) *(pending)*; ideal+FD implemented |
| Conduction | Cassisi et al. (2007) *(pending)* |
| Nuclear rates | JINA REACLIB / NACRE II *(pending)* |
| Atmosphere | model-atmosphere boundary *(pending)* |

## Licence

MIT. The opacity and equation-of-state tables are redistributed under their
own terms; see `data/*/README.md`.
