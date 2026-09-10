# Coupled evolution: equations and historical checkpoints

The latest metal-bearing EOS, shared elemental inventory, Ledoux MLT and
conservative burning/diffusion solver are documented in
[FORWARD_EVOLUTION.md](FORWARD_EVOLUTION.md). The historical experiments
below retain their original physical selections.

**Current status:** the shared GS98 metal mixture, Ledoux transport and
extended 72-cell non-grey gas atmosphere family have reached 2.85 trillion years.
Mesh and timestep refinements are complete through two trillion years:
doubling the mesh changes luminosity by +.06163%, fourfold tighter timestep
tolerances change it by -.00361%, and a fresh-start repeat gives byte-identical
stellar output. [Native restarts](RESTART.md) have passed an exact short-trajectory
comparison; the checkpoint-enabled executable repeats the full 2.5-trillion-year
track byte for byte and has continued that state through 2.75 to 2.85 trillion years.
Coupled condensate experiments remain in progress. See
[FORWARD_EVOLUTION.md](FORWARD_EVOLUTION.md) for this checkpoint and
[EXTENDED_EVOLUTION.md](EXTENDED_EVOLUTION.md) for the earlier trillion-year
calculation. The abundance/energy and coupled solver equations below still
apply. The four-plane sources and frozen
boundary described here belong to the retained `early` transport option.
The optional helium-rich non-grey atmosphere grid and its independent source
checks are described in [NONGREY.md](NONGREY.md).

## Historical 20-Gyr checkpoint

The previous 0.1 solar-mass experiment reached **20 billion years after its
specified initial static composition** using Solar Fusion II rate integrals
and finite-degeneracy Salpeter–Van Horn screening. At 4096 points,
H1 falls from .7 to .6953 and He3 rises from zero to .004740.
The star remains inside the existing EOS, opacity and frozen-atmosphere
composition bounds. This clock excludes formation and pre-main-sequence
evolution. It remains an initial main-sequence experiment, not a complete
stellar lifetime.

The new nuclear prescription replaces the evolution driver's capped classical
screening; the library's legacy default preserves static benchmarks. Its
source audit, assumptions, derivatives and independent reference integrals
are in [NUCLEAR.md](NUCLEAR.md). SVH is a smooth approximate interpolation,
with ideal Fermi electron susceptibility; it is not derived from the FreeEOS
interaction potential. This remains a significant physical uncertainty.

## Reproduce

```
cmake --build build
ctest --test-dir build --output-on-failure
mkdir -p out
build/apps/ember-evolve 4096 2e10 1e7 1 sfii-svh early > out/evolution-sfii-svh-4096-20gyr.json
build/apps/ember-evolve 1024 2e10 1e7 .25 sfii-svh early > out/evolution-sfii-svh-1024-20gyr-tight.json
python3 scripts/verify_composition_data.py
```

Arguments are mesh points, duration in years, initial timestep in years,
an optional multiplier on all timestep error tolerances (default 1), and
an optional nuclear prescription (default `sfii-svh`). The other choices
are `sfii-debye`, `sfii-legacy-screening` and `legacy`; the first two are
screening sensitivity controls. The final `early` transport argument selects
the historical EOS/opacity/atmosphere used for these references. Each run starts from its own equilibrium
at the specified initial composition. To reproduce the previous checkpoint:

```
build/apps/ember-evolve 4096 1e10 1e7 1 legacy early > out/evolution-legacy-4096-10gyr.json
```

The program emits progress on stderr and JSON on stdout; success requires
reaching the requested duration. A failed step preserves the last accepted
model and age. Default arguments give a 512-point, 100-Myr run.
The compact [20-Gyr reference](results/evolution_m010_20gyr.json),
[convergence comparisons](results/evolution_20gyr_convergence.json) and
[nuclear profile audit](results/nuclear_m010_20gyr.json) are versioned;
full profiles are in ignored `out/`. The older
[10-Gyr reference](results/evolution_m010_10gyr.json) and
[comparisons](results/evolution_convergence.json) retain their original
nuclear prescription and are historical results.

## Historical 20-Gyr numerical checks

Nineteen CTest suites pass. The new tests compare 18 bare-rate/response
queries and 10 screening queries with independent numerical integrals,
check thermal and composition derivatives in both abundance bases, and
verify limits and conservation. Coupled evolution uses the new prescription;
halving fixed backward-Euler steps over 1 Gyr reduces the radius difference
by a factor .442 at 512 points. A composition-boundary failure also returns
the unchanged previous model.

| Points | Initial R/Rsun | R/Rsun at 20 Gyr | L/Lsun at 20 Gyr | He3 at 20 Gyr |
|---:|---:|---:|---:|---:|
| 1024 | 0.1249 | 0.1251 | 0.0008403 | 0.004738 |
| 2048 | 0.1249 | 0.1251 | 0.0008406 | 0.004741 |
| 4096 | 0.1249 | 0.1251 | 0.0008408 | 0.004740 |

At 4096 points, Teff=2779 K, Tc=4.626e+6 K and rhoc=390.7 g/cm³.
The reference uses 41 accepted macrosteps, with no rejected macrosteps.
Each stores two implicit half steps. The largest recorded last-half-step
luminosity and nuclear rest-mass imbalances are 4.96e-9 and 1.65e-8 relative.
These audit the local discrete equations, not global timestep accuracy.

The 2048-to-4096 final differences are .00456% in R, .0254% in L and
3.82e-7 in absolute He3. Tightening time tolerances fourfold at 1024 points
changes final R by 4.33e-8 relative, L by 2.62e-7 and He3 by 1.31e-7 absolute.
These checks do not yet establish an asymptotic spatial convergence order.
At the final central state the electron susceptibility relative to its
classical value is .4311. The pp screening exponent is .2801 under SVH
versus .5160 under the previous classical prescription, at the same state.
The reduced ppII/pp rate ratio is at most 3.03e-7 on this profile.
Numerical differences are substantially smaller than the change of screening
prescription; that comparison does not calibrate the remaining physical error.

The atmosphere and opacity abundance guards are **unchanged**. At 20 Gyr
only about .00026 remains before the .005 surface hydrogen/He3 limits.
This limit motivated the separate extended composition sources and atmosphere
calculation documented in EXTENDED_EVOLUTION.md; the old guards remain intact.

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

The historical `early` driver option selects `NominalAbundanceOpacity`:
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
The capped classical prescription remains available as `legacy`; the new
evolution default and its remaining limitations are documented in NUCLEAR.md.
Old static numerical references predate this charge-normalization correction.

## Historical 10-Gyr checks and remaining limits

The original 10-Gyr checkpoint used legacy pp fits and classical capped
screening. Its eighteen CTest suites passed at that checkpoint. Those checks
cover abundance-basis guards,
EOS/opacity composition derivatives, independent source data, baryon and
nuclear energy conservation, stiff positive burning, disconnected mixing,
coupled thermal evolution, failed-step rollback and timestep refinement.
Halving fixed backward-Euler steps over 1 Gyr reduces the radius difference
by a factor .397 in the tested 512-point model.

| Points | Initial R/Rsun | R/Rsun at 10 Gyr | L/Lsun at 10 Gyr | He3 at 10 Gyr |
|---:|---:|---:|---:|---:|
| 1024 | .1288 | .1289 | .0009195 | .002599 |
| 2048 | .1288 | .1289 | .0009198 | .002600 |
| 4096 | .1288 | .1289 | .0009194 | .002598 |

At 4096 points, Teff=2799 K, Tc=4.556e+6 K and rhoc=357.0 g/cm³.
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

The extended calculation supplies broader composition sources, a convective
atmosphere correction and electron conduction. Physical metals/isotopes,
non-grey composition coverage and nuclear screening remain approximations.
Mesh adaptation and composition-gradient transport are needed before trusting
moving radiative-core or shell boundaries. Thermal neutrinos, diffusion and
dense-ion/solid physics remain necessary for later helium-WD cooling. Neither
checkpoint yet provides those phases.
