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

The [September 9 scientific report](docs/reports/2026-09-09/ember_status_and_future.pdf)
([LaTeX and recovery artifacts](docs/reports/2026-09-09/README.md)) covers the
completed track, input validation limits and the remaining lifetime/remnant
physics. [The next-session prompt](docs/NEXT_SESSION_PROMPT.md) and the newest
entry in [HANDOFF.md](HANDOFF.md) describe where to resume.

The new opt-in **GS98 metal-bearing EOS and Ledoux transport** are described
in [docs/FORWARD_EVOLUTION.md](docs/FORWARD_EVOLUTION.md). That 512-point
model reaches **2.85 trillion years** at R=.148761 Rsun, L=.00185757 Lsun and
Teff=3106.84 K, with XH=.30289 and fully convective structure. The extended
72-cell gas atmosphere family supplies the boundary. At two trillion years,
an independent repeat reproduces the entire stellar output byte for byte;
doubling the mesh changes luminosity by +.06163%, and fourfold tighter time
tolerances change it by -.00361%. Semiconvective and thermohaline composition
transport are implemented. [Native restarts](docs/RESTART.md) reproduce the
subsequent short test trajectory exactly. The checkpoint-enabled executable
also repeats the entire 2.5-trillion-year track byte for byte, then successfully
continues its saved state through 2.75 to 2.85 trillion years. The settled-grain condensate
radiation experiment remains in progress.
The earlier validated controls and their numerical refinements follow.
The remaining requirements for complete lifetimes and a mass–metallicity
survey are in [docs/LIFETIME_SURVEY.md](docs/LIFETIME_SURVEY.md).

Built for the lowest-mass stars first; the interfaces are shaped so
that the ultra-cold evolution of *massive* remnants - a solar remnant, or
something near the Chandrasekhar mass - can be added as further terms rather
than as a rewrite. See `docs/ROADMAP.md`.

The code now completes a **coupled one-trillion-year evolution experiment
for a 0.1 Msun star**, including explicit He3 burning, instantaneous convective
mixing and thermal evolution. The model remains fully convective, with
hydrogen reduced from .7 to about .533 by baryonic mass. Twenty-four test
suites pass, including independent source queries and conservation checks.
At 2048 points the final model has **R=.143515 Rsun, L=.00134953 Lsun,
Teff=2920.28 K**. Mesh and timestep comparisons accompany the reference.

The extended model uses eight FreeEOS composition potentials, broader
AESOPUS/TOPS X/Z opacity families with helium isotope number-density mapping,
electron conduction, and a convective composition correction anchored to the
solar AMES-COND atmosphere. An optional **helium-rich non-grey grid** now
supplies 48 independently computed radiative/convective gas atmospheres,
interpolated in hydrogen, helium-3, effective temperature and gravity.
Its source checks and remaining approximations are documented in
[docs/NONGREY.md](docs/NONGREY.md). With that grid, a 512-point model reaches
one trillion years at **R=.143578 Rsun, L=.00133023 Lsun, Teff=2909.13 K**,
remaining fully convective. The stellar run takes 10.2 minutes on an M4 Max;
the expensive atmosphere calculation is performed beforehand. Metals, isotope
cross sections and screening retain documented approximations. Time starts
from a specified static composition, excluding formation; hydrogen exhaustion,
the blueward turn and white-dwarf cooling remain future milestones.

See [docs/EXTENDED_EVOLUTION.md](docs/EXTENDED_EVOLUTION.md) for reproduction,
source coverage, numerical comparisons and the remaining physical inputs.
[docs/EVOLUTION.md](docs/EVOLUTION.md) gives the evolution equations and
historical 10/20-Gyr results; [docs/NUCLEAR.md](docs/NUCLEAR.md) describes the
SFII/SVH burning prescription. [docs/FREEEOS.md](docs/FREEEOS.md) documents
the potential and its source-fit uncertainties. The retained
[CMS19 experiment](docs/CMS19.md) keeps its static-only energy guard.

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
| High-T opacity | LANL TOPS ATOMIC at seven X values from .3 to .75, Z=.01/.02/.03; OPAL GS98 subset also available; strict smooth blends |
| Equation of state | FreeEOS 3.0 EOS1 represented by one C2 Helmholtz potential; GS98 metal-bearing H/He3 family, isotope corrections; historical H/He proxy and static alternatives retained |
| Conduction | Ioffe pure-ion tables; Cassisi et al. (2007, 2021), Blouin et al. (2020); approximate mixture resistivities |
| Convection | Böhm–Vitense MLT; Schwarzschild or full-EOS Ledoux buoyancy; semiconvective and thermohaline composition transport |
| Nuclear rates | Reduced pp network, SFII S-factor quadrature, finite-degeneracy SVH screening, atomic mass-defect heating; legacy prescription retained |
| Atmosphere | TLUSTY208/SYNSPEC54 composition-dependent non-grey gas grid at tau=100; AMES-COND anchor plus grey convective composition correction remains the default control |

## Licence

MIT. The opacity, equation-of-state and conductivity data are redistributed under their
own terms; see `data/*/README.md`.
