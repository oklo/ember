# Composition-dependent physics and first coupled evolution

The 0.1 solar-mass experiment now evolves composition and thermal structure
together. A 4096-point run reaches **10 billion years after its specified
initial static composition**, with H1 falling from .7 to .69740182 and He3
rising from zero to .00259806. This clock does not include formation or a
pre-main-sequence calculation, and the initial zero-He3 model is an explicit
initial condition. It is an early main-sequence experiment with bounded
mixture approximations, not a complete stellar lifetime.

## Reproduce

```
cmake --build build
ctest --test-dir build --output-on-failure
mkdir -p out
build/apps/ember-evolve 4096 1e10 1e7 > out/evolution-4096-10gyr.json
build/apps/ember-evolve 1024 1e10 1e7 .25 > out/evolution-1024-10gyr-tight.json
python3 scripts/verify_composition_data.py
```

Arguments are mesh points, duration in years, initial timestep in years,
and an optional multiplier on all timestep error tolerances (default 1).
The program emits progress on stderr and JSON on stdout; success requires
reaching the requested duration. A failed step preserves the last accepted
model and age. Default arguments give a smaller 512-point, 100-Myr run.
The compact [reference track](results/evolution_m010_10gyr.json) and
[convergence comparisons](results/evolution_convergence.json) are versioned;
full profiles are in ignored `out/`.

## Abundance and energy conventions

`Composition` records its abundance basis. Existing static configurations
retain atomic mass fractions. Evolution uses baryonic abundances
`X_i = A_integer,i * Y_i`, with integer mass numbers 1, 3, 4, etc. Thus
`sum(X_i)=1` expresses conserved baryons. Rates consume three H1 per pp
reaction including the fast deuteron capture, and explicitly evolve He3
through ppI and the reduced ppII branch. All nuclear mass defects use the
physical **atomic masses**, including when the abundance basis is baryonic:

```
eps_deposited + eps_nuclear_neutrinos = -c^2 sum[(A_atomic/A_integer) dX_i/dt]
```

No abundance clipping or renormalization hides a nuclear mass defect.
The evolution model's density and fixed enclosed mass are baryonic;
Newtonian gravity uses that conserved mass. It does not feed the small
binding-energy change back into gravitational mass. For this reason a
nominal baryonic .1 Msun, X=.7 model is not precisely the same initial
particle inventory as the historical atomic-mass .1 Msun, X=.7 model.

The thermal energy equation is evaluated with the actual old and new
compositions:

```
dL/dm = eps_deposited - [u_new-u_old + P_new*(1/rho_new-1/rho_old)]/dt
```

`u` excludes nuclear rest mass but includes the EOS's composition-dependent
thermal, ionization and interaction energies. Nuclear neutrinos have already
been subtracted from deposited heating. Trapezoidal nodal mass weights are
shared by energy, burning and mixing; the unresolved central sphere adds
its mass to the first node. The reported luminosity balance is an audit of
this **discrete local first law**, not an independent global gravitational
energy integral or a claim of timestep accuracy.

## Composition-dependent EOS

`CompositionHelmholtzEos` combines four complete material Helmholtz
potentials generated from FreeEOS 3.0.0 EOS1 at source hydrogen atomic mass
fractions **.60, .65, .70 and .75**. Each uses the original .0125-dex
temperature/density grid and strict source-stencil stability masks.
Potentials, rather than independently tabulated pressure and energy, are
interpolated linearly in composition. Thermodynamic identities therefore
hold at every supported mixture; composition derivatives can jump at knots.
The runtime exposes analytic pressure and energy derivatives when H1 or He3
replaces He4 at fixed T, density and metals.

He3 changes the ion/electron number density. Let `h=XH/wH`, `n3=X3/w3`,
`n4=(X4+Z)/w4`, where `w` is the selected abundance weight. With physical
atomic masses `aH,a4`, the equivalent source state is
`s=aH*h+a4*(n3+n4)`, `Xsource=aH*h/s`, `rho_source=s*rho`.
The material `F/T` is `s` times the interpolated source potential at that
state. An analytic ideal helium isotope mixing/translation entropy term
distinguishes He3 and He4. Tiny isotope shifts of electronic levels and
nuclear-spin entropy constants are omitted. Radiation is added once at
the actual model density. Metals remain the explicit He4 proxy.

Twenty-eight independent FreeEOS queries use elemental mole abundances
`eps_H=XH`, `eps_He=X3/3+(1-XH-X3)/4` directly, at compositions and thermal
states between source knots. Maximum sampled differences are .00872% in
P/E and .0540% in cv/cp/adiabat. These queries validate the number-density
mapping and composition interpolation; they do not independently validate
the isotope entropy term or the metal approximation. The wider original
source-fit audit still finds 1–2% response differences at some molecular
joins; see [FREEEOS.md](FREEEOS.md).

All 723,460 source evaluations across the four grids converged. Their raw
responses are archived, and all four potentials reproduce byte for byte.
Checksums are in `data/eos/sources/freeeos300_composition_manifest.json`.
The independent reference can be regenerated with:

```
python3 scripts/generate_freeeos_composition_reference.py /path/to/probe /tmp/composition-reference.dat
```

The external probe build and source archive checksum are in FREEEOS.md.
Normal builds and reimport need no Fortran or network.

## Opacity and atmosphere coverage

AESOPUS already spans hydrogen abundance. New TOPS calculations supply
the same four X planes at fixed GS98 Z=.02. Every returned H, He and metal
fraction is verified; every server-substituted density is excluded.
Two common rectangles retain 4×50×63 low-temperature cells and 4×36×71
high-temperature cells. Interpolation uses native log density, with the
existing temperature blends. `dln(kappa)/dX` is analytic through both blends.
Sixty independent original source cells test all four planes at runtime.
Requests, exact normalized text and SHA-256 hashes are archived under
`data/opacity/sources/`; reimport is entirely offline.

The evolution driver explicitly selects `NominalAbundanceOpacity`:
nominal baryonic X/Z are used in atomic-mixture opacity tables, with He3
treated as He4. This omits isotope opacity effects and the small
atomic/baryonic mass-convention correction to X, Z and density. The wrapper
rejects He3 above .005. The abundance caps are conservative operational
limits for this experiment, **not calibrated opacity-error bounds**.

`FrozenCompositionAtmosphere` keeps the original GN93 solar COND tau=100
temperature/pressure boundary, while recomputing density from the actual
composition-dependent EOS. It rejects `abs(Xsurface-.7)>.005`, He3>.005,
changed metals, and the existing Teff/gravity source bounds. Thus this run
does not claim a composition-dependent non-grey atmosphere. The underlying
fixed-composition reader still rejects changed abundances or abundance basis.

## Coupled burning and convective mixing

At each timestep the solver finds connected Schwarzschild-unstable regions
using the same midpoint radiative and adiabatic gradients as the structure
equations. Each connected region solves its integrated pp abundance
equations by backward Euler, assuming instantaneous homogeneous mixing.
Stable isolated nodes burn locally. The abundance Newton solve uses analytic
reaction and screening composition derivatives and a positivity-preserving
line search. Metals are inert but mix conservatively.

Burn/mix and thermal structure are iterated until both converge and the
convective partition stops changing. This is a converged block iteration,
not a single unevaluated operator-split update. The 0.1 Msun test remains
fully mixed. Disconnected-region tests verify conservation without mixing
across stable boundaries. This does not validate moving radiative-core
boundaries: Ledoux stability, diffusion, semiconvection and mesh adaptation
remain pending.

The timestep controller compares one full step with two half steps and
accepts the two-half-step solution. Default error scales are 1e-5 in log
structure, 1e-8 in each absolute active-species abundance, and 1e-4 in
relative surface luminosity. These are local estimates, not global error
bounds. A maximum .001 abundance change also guards each implicit step.

The nuclear audit corrected an older screening normalization error:
the classical Debye charge density is `ne + sum(ni Zi^2)`; multiplying
ion-averaged charges by electron density introduced an extra mean ionic
charge. The corrected expression and composition derivatives are tested
independently. See the [Debye expression in Phyu et al. (2020)](https://academic.oup.com/ptep/article/2020/9/093D01/5899683).
The retained classical weak-screening model and exp(2) cap still require
qualification for intermediate coupling and partially degenerate electrons;
the pp coefficients and reduced ppII treatment also need a modern rate audit.
Old static numerical references predate this screening correction.

## Numerical checks and limits

Eighteen CTest suites pass. Added checks cover abundance-basis guards,
EOS/opacity composition derivatives, independent source data, baryon and
nuclear energy conservation, stiff positive burning, disconnected mixing,
coupled thermal evolution, failed-step rollback and timestep refinement.
Halving fixed backward-Euler steps over 1 Gyr reduces the radius difference
by a factor .397 in the tested 512-point model.

| Points | Initial R/Rsun | R/Rsun at 10 Gyr | L/Lsun at 10 Gyr | He3 at 10 Gyr |
|---:|---:|---:|---:|---:|
| 1024 | .12884233 | .12891572 | .0009194957 | .0025985795 |
| 2048 | .12882960 | .12890799 | .0009197530 | .0025996158 |
| 4096 | .12884559 | .12891498 | .0009194175 | .0025980645 |

At 4096 points, Teff=2799.32 K, Tc=4.55637e6 K and rhoc=356.961 g/cm³.
The reference uses 23 accepted macrosteps, each storing two implicit half
steps. Its largest recorded last-half-step luminosity imbalance is 4.70e-9
relative, and nuclear mass-defect imbalance is 1.08e-8. Each recorded
diagnostic belongs to that half step, not the whole macrostep integral.

The 2048-to-4096 final differences are .00542% in R, .0365% in L and
1.55e-6 in absolute He3. Mesh changes remain nonmonotonic. Tightening all
time-error scales fourfold at 1024 points changes final R by 1.57e-8
relative, L by 1.02e-7 relative and He3 by 6.21e-8 absolute.
These numerical differences are much smaller than the presently unmeasured
physical uncertainties in mixtures, atmospheres, screening and convection.

Before substantial hydrogen depletion: extend source composition coverage,
improve physical metals/isotopes, replace the frozen atmosphere with an
appropriate composition grid, qualify nuclear rates/screening, and add
mesh adaptation and composition-gradient transport. Later helium-WD cooling
also needs conduction, thermal neutrinos, diffusion and dense-ion/solid
physics. This milestone does not yet provide those phases.
