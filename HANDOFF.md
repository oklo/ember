# ember — handoff

Updated 2026-09-07 after composition-dependent EOS/opacity and coupled
burning/mixing. Read this, then **`docs/EVOLUTION.md`**, `docs/FREEEOS.md`
and the data READMEs. Use `git status` and `git log` for current history.
GitHub pushes are authorized. The previous static/non-grey checkpoint is
`6cc4fb4`; its numerical values predate the screening correction below.

**The 0.1 Msun experiment now evolves for 10 billion years with coupled pp
burning, instantaneous convective mixing and thermal structure.** At 4096
points: R=.12891498 Rsun, L=.00091941748 Lsun, Teff=2799.32 K,
Tc=4.55637e6 K, rhoc=356.9611 g/cm³. H1 falls from .7 to .69740182;
He3 rises from zero to .00259806. Eighteen CTest suites pass.

**This is bounded initial evolution.** The age origin is the specified
static composition, not formation. New evolution uses conserved **baryonic**
mass and abundances; old static drivers retain atomic mass fractions.
Nuclear rest mass defects are separately included in energy, with no
renormalization that hides their release. Newtonian gravitational mass
remains the conserved baryonic mass. Do not silently equate the old and
new mass/composition conventions.

Four FreeEOS potential planes cover source X=.6/.65/.7/.75, with explicit
He3 number-density mapping and ideal isotope entropy. Metals remain the
He4 proxy. New TOPS data cover the same X planes at GS98 Z=.02, rejecting
every source-substituted density. Evolution explicitly uses nominal X/Z
and He3-as-He4 opacity, capped at He3=.005. COND remains a **frozen GN93
T/P boundary** with actual-EOS density inversion, limited to
|Xsurface-.7|<=.005 and He3<=.005. These limits are operational assumptions,
not calibrated physical-error bounds. No composition-dependent atmosphere
source grid or white-dwarf cooling run exists yet.

The new nuclear audit corrected the old classical Debye screening charge
normalization (an extra mean ionic charge). This changes even static runs;
`docs/results/equilibrium_m010_screening_corrected.json` records the updated
atomic-basis reference. Classical weak screening and the exp(2) cap still
need qualification for intermediate coupling/degenerate electrons.
FreeEOS molecular-fit joins still limit sampled response accuracy to about
1–2%. Exact discrete conservation is not a physical accuracy estimate.

The 10-Gyr 4096-point track uses 23 accepted macrosteps. Largest recorded
last-half-step luminosity and nuclear mass-defect imbalances are 4.70e-9
and 1.08e-8 relative. Mesh differences remain nonmonotonic; timestep and
mesh comparisons are in `docs/results/evolution_convergence.json`.

---

## 1. What this is, and why it exists

`ember` is a ground-up C++23 stellar evolution code for the **lowest-mass
stars** — objects that burn hydrogen for trillions of years, turn *blue*
rather than red when their fuel runs out, and end as helium white dwarfs that
outlive everything else in the galaxy.

It is the successor to a FORTRAN line living at
`/Users/greglaughlin/Projects/low_mass_stars` (GitHub: `oklo/Henyey`), which
reconstructed the code of **Laughlin, Bodenheimer & Adams (1997)**. That
reconstruction succeeded and is finished; it stays as a historical artifact
with its period-styled website. `ember` is where the physics moves forward.

**Motivation from the user (Greg Laughlin):** he is discussing with **Fred
Adams** a full-scale update of their **1997 Rev Mod Phys** paper. That update
will want the ultra-cold evolution of *massive* white dwarfs — a solar
remnant, possibly objects near the Chandrasekhar mass. Those are **not to be
built yet**, but the structure must accept them smoothly. Section 6 covers
how.

**The milestone that defines "working":** a 0.1 M☉ star evolved end to end —
Hayashi track, trillions of years of hydrogen burning, the blueward turn, and
down the helium-white-dwarf cooling track below 10⁻⁶ L☉ — in one run.
**GitHub push is now authorized.** On 2026-09-07 the user explicitly asked
to push the current checkpoint and continue development, superseding the
earlier requirement to wait for the end-to-end run. The scientific milestone
is unchanged. The user then explicitly confirmed creation of the private
repository **https://github.com/oklo/ember**. `master` tracks `origin/master`;
use `git log` and the configured remote for current history.

---

## 2. Build and test

```
cd /Users/greglaughlin/Projects/ember
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

Host is an **Apple M4 Max**, 10 performance + 4 efficiency cores. `clang`
accepts `-mcpu=apple-m4`. Accelerate is linked and available, while the
current small Henyey blocks use a portable pivoted kernel. CMake and Ninja
were installed via Homebrew this session.

`-ffast-math` is **deliberately not used** — it licenses the compiler to
assume no NaN/Inf, and a stellar model legitimately probes states where a
table returns one. Those must surface.

---

## 3. Current state

Eighteen suites pass: eos, opacity, dense_opacity, nuclear, structure,
convection, atmosphere, henyey, relaxation, stellar_equilibrium, cms19,
tops, low_mass_equilibrium, helmholtz, nongrey_equilibrium,
composition_physics, burning_mixing and coupled_evolution. The independent n=3 radiative-polytrope
benchmark and the older ionized-EOS 0.5 Msun benchmark are retained.

| Module | Implementation | State |
|---|---|---|
| Constants/composition | `constants.hpp`, `composition.{hpp,cpp}` | CODATA 2018, nominal solar units, 8 species, AAG21 metal helper |
| Analytic EOS | `eos_components.cpp`, `eos_composite.hpp`, `fermi.{hpp,cpp}` | ions, radiation, relativistic FD electrons; component Hessians |
| CMS19 EOS | `eos_cms19.{hpp,cpp}`, `src/jet2.hpp` | original pure-H/He TP density/entropy, additive-volume mixture, actual interpolation Hessians, strict fluid support; **static only** |
| Helmholtz EOS | `eos_helmholtz`, `eos_composition` | C2 FreeEOS potentials at four X values, He3 number-density mapping, analytic thermal/composition derivatives; explicit metal/isotope approximations |
| EOS interface | `eos.hpp` | transport response derivatives, optional density bounds, virtual PT inversion, explicit internal-energy availability |
| Opacity | `opacity_table.{hpp,cpp}`, source wrappers | strict monotone tables in native log R or log rho, fixed Z, single- or multiple-X support |
| Low-T opacity | `opacity_aesopus.hpp`, `opacity_ferguson.hpp` | AESOPUS 2.1 gas log R to 6; Ferguson with grains also available |
| Hot opacity | `opacity_opal.hpp`, `opacity_tops.hpp` | OPAL log R to 1; TOPS two un-clamped rectangles at X=.6/.65/.7/.75, Z=.02 |
| Blending | `opacity_blend.{hpp,cpp}` | smooth ln-kappa blends, analytic derivatives, strict overlap intersection |
| Nuclear | `nuclear_pp.cpp` | pp chains, explicit He3, atomic/baryonic bases, corrected classical screening and analytic composition/thermal derivatives |
| Structure | `structure.{hpp,cpp}` | four residuals and analytic Jacobian; independent numerical checks; unsupported positive-dt energy rejected |
| Convection | `convection.{hpp,cpp}` | BV58 MLT, Schwarzschild criterion, bounded cubic and analytic response |
| Atmosphere | `atmosphere_grey.cpp` | radiative Eddington T(tau), variable opacity, adaptive integration and analytic sensitivities; EOS/opacity bounds intersect |
| Atmosphere table | `atmosphere_table.{hpp,cpp}` | AMES-COND tau=100 physical states; untouched MESA cells; strict bounds and explicit solar-mixture proxy; synthetic tests retained |
| Boundaries | `boundary.{hpp,cpp}` | regular unresolved central sphere and interchangeable atmosphere surface, analytic derivatives |
| Henyey/relaxation | `henyey.cpp`, `relaxation.cpp` | pivoted blocks, iterative refinement, damped Newton with residual AND undamped correction checks |
| Seeds/apps | `examples/stellar_seed.hpp`, `apps/equilibrium.cpp` | n=3 or n=1.5 Lane–Emden seed; rebuild mass after envelope adjustment; structured diagnostics |
| Evolution | `evolution.{hpp,cpp}`, `apps/evolve.cpp` | coupled backward-Euler burning/mixing/thermal structure; adaptive step-doubling on fixed mass mesh; bounded 10-Gyr experiment |
| Massive WD hooks | `losses.hpp`, `conduction.hpp` | declared, unimplemented; not the current task |

Data travel with the code. Source hashes, strict coverage and reproduction
commands are in `data/eos/README.md`, `data/atmosphere/README.md` and
`data/opacity/README.md`. The FreeEOS source builder/generator/importer and
COND importer reproduce the new data; all 180,865 direct source evaluations
were reproduced byte for byte from a clean FreeEOS build. Normal builds
need no Fortran or network. The earlier six CMS19/opacity imports were reproduced byte for byte with stdlib Python scripts:

- `scripts/import_opal.py`: 9,880 original GS98 Z=.02 cells, log T=4..7.1,
  log R=-8..1, 10 X planes. The original source is ragged outside this box.
- `scripts/import_aesopus.py`: 149,810 original AESOPUS 2.1 **gas** cells,
  log T=2..4.5, log R=-8..6, 10 X planes. Archive filenames retain 2.0 but
  headers identify 2.1; selected tables do not include grains.
- `scripts/import_cms19.py`: 121×441 original (log rho, log S) pairs per
  pure component on the source TP grid. Runtime masks reject unphysical
  corners. No source energy values are imported or repaired.
- `scripts/import_tops.py`: 3150 low-rectangle and 2556 high-rectangle cells
  (overlapping), native rho grid, only X=.7, Z=.02. Every server-substituted
  density is excluded. Exact returned text and request are versioned under
  `data/opacity/sources/`; a fresh service request is not required to build.

The generic `TabulatedOpacity::Range` now names its density fields
`logD_min/max` and includes the `DensityAxis`. log R requires the -3 term
in the temperature derivative; native log rho does not. All strict source
bounds remain enforced. TOPS itself blends rectangles over log T=5.6..5.7;
the stellar driver blends AESOPUS to its selected hot opacity over 4.4..4.5.

---

## 4. Reproduce and continue

### Current coupled evolution reference

```
build/apps/ember-evolve 4096 1e10 1e7 > out/evolution-4096-10gyr.json
build/apps/ember-evolve 1024 1e10 1e7 .25 > out/evolution-1024-10gyr-tight.json
python3 scripts/verify_composition_data.py
```

Positions are points, duration years, initial step years, optional tolerance
multiplier. JSON stdout includes accepted history and full final profile;
stderr reports progress. Each macrostep compares one full versus two half
steps and retains the two-half-step solution. Failed steps roll back; do
not accept the last unconverged composition or structure. Convection uses
connected Schwarzschild regions and instantaneous mass-conserving mixing,
iterated with the thermal solve until both converge. See EVOLUTION.md.

Raw data for all four EOS composition planes and all TOPS planes are
versioned with checksum manifests. The EOS and opacity imports reproduce
byte for byte offline. The original external FreeEOS probe is available
in `/tmp/ember-freeeos-reproduce/probe` on the development host; build
instructions are in FREEEOS.md. Independent reference generation:

```
python3 scripts/generate_freeeos_composition_reference.py /tmp/ember-freeeos-reproduce/probe /tmp/composition-reference.dat
```

### Historical FreeEOS/COND static reference at `6cc4fb4`

```
build/apps/ember-equilibrium 4096 .1 .15 --eos freeeos --hot-opacity tops --seed-index 1.5 --atmosphere cond-solar-proxy > out/equilibrium-m010-nongrey-reference.json
```

The following values belong to `6cc4fb4` before screening correction;
the same command now gives R=.12935124, L=.00094042956, Teff=2810.43 K.
The old checkpoint had nine accepted updates, residual 6.60e-13, correction 6.76e-12, nuclear
balance about 2e-14 relative and independent virial error -6.396e-7.
The 1024/2048/4096 models agree to <.014% in R and <.078% in L between
successive meshes; virial error falls by four each doubling, while R/L
changes are nonmonotonic. Versioned metadata and the separate EOS/atmosphere
comparison are in `docs/results/equilibrium_m010_nongrey.json` and
`docs/results/eos_atmosphere_comparison.json`.

FreeEOS source options are EOS1 `(3,1,-2)`, H=.7/He=.3, all metals zero.
The C2 potential grid is .0125 dex in log T/log Q with strict stability
masks. Read `docs/FREEEOS.md` before modifying it: the H2 partition-function
join at 9000 K means local source identities alone are insufficient.
Energy consistency is exact for the implemented potential; the source
response agreement has a separate, sampled 1–2% precision qualification.
A fresh external FreeEOS process per isotherm avoids a source cache failure
at dense-to-dilute row resets; never accept nonzero source `info`.

COND uses only untouched Teff=1800..3300 K, log g=3.5..6 MESA cells, matched
at Rosseland tau=100. `--tau-top` is grey-only. The seed must use the
atmosphere's local T, not Teff. Do not relabel its native GN93 mixture as an
exact match to the interior, or silently extrapolate outside this rectangle.

### Retained CMS19/grey reference

```
mkdir -p out
build/apps/ember-equilibrium 4096 .1 .15 --eos cms19 --hot-opacity tops --seed-index 1.5 > out/equilibrium-m010-reference.json
```

Use 2048 points for a smaller working model. Positional arguments are
points, mass/Msun, seed radius/Rsun. Radius is free to relax. `--tau-top`
can be set explicitly; default .001 for CMS19, 1e-6 for the ionized EOS.
All composition is fixed at X=.7, Y=.28, Z=.02, He3=0. CMS19 explicitly
represents metals as helium (effective Y=.3). MLT alpha=1.9 is uncalibrated.
The surface is a radiative grey atmosphere matched at tau=2/3. No grain
opacity, conduction, composition mixing or age advancement is included.

At 4096 points: 12 updates, residual 4.10338e-10, undamped correction
3.06322e-10, R=.12860617 Rsun, L=.00097562953 Lsun, Teff=2844.5719 K,
Tc=4.5663102e6 K, rhoc=366.87522 g/cm³. Nuclear luminosity agrees with
surface L to ~2e-15 relative; independent virial error is -6.5848e-7.
The unresolved center is 1.23e-10 of the mass. All profiles are ordered.

| Points | R/Rsun | L/Lsun | Teff (K) | Absolute virial error |
|---:|---:|---:|---:|---:|
| 128 | .12911938 | .0009629408 | 2829.64 | 6.856e-4 |
| 256 | .12935265 | .0009420455 | 2811.62 | 1.699e-4 |
| 512 | .12885921 | .0009637381 | 2833.08 | 4.229e-5 |
| 1024 | .12862818 | .0009746770 | 2843.63 | 1.055e-5 |
| 2048 | .12861667 | .0009751461 | 2844.10 | 2.635e-6 |
| 4096 | .12860617 | .0009756295 | 2844.57 | 6.585e-7 |

Coarse thermal results are nonmonotonic. Fine meshes agree to ~.1% in L
and .02% in R, but R/L do not yet show clean second-order convergence;
virial error does. Do not call a tiny nonlinear residual spatial accuracy.
Doubling the atmosphere starting column to .002 changes the independent
256-point R/L by .000402%/.00259%; a same-mass-mesh 2048-point perturbation
also passes in `test_low_mass_equilibrium`. This does not validate grey
atmospheric physics. Full current JSON/log is in ignored `out/`; compact
reference metadata are versioned in `docs/results/equilibrium_m010.json`.

### Retained CMS19 source audit — its energy guard still applies

`docs/CMS19.md` gives the implementation and independent source audit.
CMS19 response derivatives follow the **actual interpolated entropy**.
The Maxwell identity D=P*delta/(rho*T*cp*grad_ad)-1 is not artificially
forced to zero. Max |D|=.24178 in the 4096-point star at T=27304 K,
rho=.12394 g/cm³, m/M=.99999310; mass RMS=.0043487. Small outer mass is
not evidence of small influence on R or L.

`python3 scripts/audit_cms19_energy.py /tmp/ember-eos2019.tar.gz` reproduces
both defects directly from original source columns:

- Pure H at log T=4.45, log P[GPa]=1.35 has D=-.26498098, before ember
  interpolation. This is a source-version consistency issue as well as an
  interpolation accuracy question.
- Pure He at rho=1 g/cm³ drops in internal energy from 8.88587e13 to
  8.22697e13 erg/g between T=891251 and 1e6 K. Secant dU/dT=-6.05890e7
  erg/g/K while listed entropy cv=+9.72478e7.

Do not fabricate internal energy, silently repair table cells, loosen the
Maxwell check, or replace an entropy-derived adiabat just to make an
identity pass. Exploratory U-TS/free-energy and pressure-integrated-potential
reconstructions developed negative heat capacities at source joins and
were rejected. They exist only as ignored scratch, not library physics.
`Cms19Eos` has no caloric state or He3 support; energy and unavailable
abundances are NaN and dt>0 zones/center throw explicitly.

The runtime conservatively masks whole 4×4 stencils by source density,
fluid/quantum limits, and log T=3.2..7.3. It is a selected computational
subset, not proof of accuracy everywhere inside. Grey tau_top=1e-6 falls
below its density support; defaulting CMS19 to .001 is explicit and tested.
Never extrapolate a surface integration below the EOS density floor.

### Next actions, in order

1. Qualify pp coefficients and screening for the current partially degenerate
   core; replace the classical capped approximation with a suitable consistent
   prescription. The charge-normalization bug is fixed, not the full screening
   uncertainty. Keep rate normalization and atomic/baryonic conventions explicit.
2. Extend source composition coverage toward hydrogen depletion and improve
   physical metals/isotope EOS treatment and source-fit joins. Do not clamp
   the family at its current equivalent-source X limits.
3. Obtain a composition-dependent non-grey atmosphere or quantify a justified
   extension. The present frozen COND T/P wrapper must stop at its explicit
   abundance caps; deleting those guards does not supply new atmosphere physics.
4. Adaptive mesh, Ledoux/diffusive mixing and moving convective-boundary
   validation before following a radiative core or burning shell. Repeat
   age/composition-aware literature/observational comparisons with the new track.
5. Conduction, thermal neutrinos and dense-ion/cooling physics for the original
   end-to-end 0.1 Msun milestone. Massive WD physics remains future work.
   GitHub pushes remain authorized, as recorded in §1.

### Retained old runs and practical source notes

```
build/apps/ember-polytrope 256 > out/polytrope-256.json
build/apps/ember-equilibrium 512 .5 .6 > out/equilibrium-m050.json
build/apps/ember-equilibrium 128 .1 .2 > out/old-ionized-m010-failure.json
```

The first two converge. The last uses the old ionized EOS/OPAL/n=3 defaults
and still fails at the OPAL density edge; it is retained as an explicit
failure diagnostic. The .5 Msun case has R=.94219963 Rsun, L=.019877124
Lsun, Teff=2232.77 K; its inflated cool envelope is not realistic. The n=3
stellar seed also proved poor for the new CMS19 .1 Msun model; n=1.5
converges. This changes only initialization, never the solved equations.

Original EOS archive URL:
`https://perso.ens-lyon.fr/gilles.chabrier/DirEOS/DirEOS2019.tar.gz`.
The source README and 2021 README recommend the 2019 IVL tables for stars,
2021 interactions for brown dwarfs; never use 2021 effective H as pure H.
Temporary originals are `/tmp/ember-eos2019.tar.gz`,
`/tmp/ember-eos2021.tar.gz`, `/tmp/ember-aesopus21-gs98.zip`, and
`/tmp/ember-GS98hz.gz`. Hashes and commands are documented with the data.
They are not needed for normal compilation or tests.

TOPS: `/results` can return a stale prepared calculation despite new
parameters. Submit through `/submit`, then verify returned composition and
grid dimensions before import. The final archived request uses the public
form identifiers, a blank mixture name, and the exact 21-element GS98 mass
mixture. Its returned warnings identify clamped densities; every such
pair is excluded. Request-specific curl `-k` was needed for the local
certificate-chain failure. All data needed for TOPS reimport are versioned.

---

## 5. Hard-won lessons — these cost real days in the FORTRAN line

Do not rediscover these.

1. **Fixed-format input is a trap.** `NPRIN=100000` is six digits in an `I5`
   field; it shifted every later field and silently ate the Hayashi-start
   flag. A whole session was spent misdiagnosing the resulting physics as an
   "MLT regression". *`ember` must never column-count input.*
2. **Fixed-width output is the same trap.** A model counter in `I4` printed
   `****` past 9999 and corrupted downstream analysis — I misreported where
   runs ended because of it. *Structured output only.*
3. **An equation form can be valid only in a limit.** Scaling the radiative
   equation by the convective efficiency φ = ∇/∇_rad is sound *only where φ is
   order unity*. Deep in a giant envelope φ ≈ 1e-6 — perfectly efficient
   convection — and multiplying the matrix row by a millionth decoupled those
   zones' temperature from the luminosity. **This was the cool-giant wall that
   stopped every giant for weeks.** Fix: superadiabatic zones use the plain
   gradient form, well conditioned at any efficiency. *Where a formulation has
   a domain, assert it.*
4. **Constants tuned for the Sun break elsewhere.** A fitting point at a fixed
   mass fraction drifts *above the photosphere* when a star swells to 300 R☉
   (measured: T(N)/Teff = 0.86, i.e. cooler than the effective temperature).
   Anything that is "a good value for the Sun" is suspect.
5. **Silent table extrapolation is lethal.** Below its floor the 1983 opacity
   table returned κ ~ 1e-9; the envelope went transparent and cool giants lost
   their Hayashi limit entirely, expanding to 314 R☉ and breaking. *In `ember`,
   `FergusonOpacity` throws outside its table. Keep it that way.*
6. **Write tests as physics statements, not numbers.** Every test caught
   something real this session:
   - electron normalisation short by √π (caught by `n = 2e^η/λ³`)
   - NaN under strong degeneracy: at β ~ 1e-4 the FD derivative integrand is a
     shell of width ~β that a single quadrature panel never samples
   - a wrong ppII branch ratio (caught by the mass defect, **not** by
     "mass fractions sum to zero" in the atomic-mass convention, where
     Σ dX/dt = −(ε+ε_ν)/c². In the new baryonic convention their sum is
     zero and atomic/integer mass ratios give the released energy)
   - **and twice the *test* was wrong, not the code**: a point labelled "the
     ideal limit" with 13% of its pressure in radiation, and a shell built by
     one-point integration checked against a centred-difference equation.
     Suspect the test too.

Two more lessons from connecting the atmosphere in ember:

- **Differentiate the actual interpolant.** Interpolating opacity derivatives
  separately missed the dependence of the temperature slope limiter on the
  input opacities. The full Ferguson atmosphere test saw a ~0.6% mismatch in
  its Teff pressure response; a direct opacity regression also failed. Both
  pass after propagating derivatives through the limiter. No opacity values
  were intentionally changed.
- **An inversion must report failure.** A small pressure residual can conceal
  a bad density when radiation dominates, and an exhausted iteration must not
  return a density as a success. `rho_from_PT` now checks both and throws.

Two more from assembling the analytic zone Jacobian:

- **Accurate values do not guarantee accurate second derivatives.** The
  original 64-node Fermi shell produced plausible pressures but roughly 0.2%
  errors in degenerate response derivatives, unchanged when the finite
  difference step shrank. Splitting the shell at its center resolved them.
  Use centered number-weighted kernels for the density-constrained electron
  Hessians, and preserve `f(1-f)` even when `f` rounds to one. The derivation
  is in `docs/JACOBIAN.md`; do not replace it with large cancelling terms.
- **Differentiate the implemented heating rate.** Approximate reaction
  energy weights and omitted screening derivatives did not differentiate the
  pp mass-defect heating value. Its derivatives now use the same masses,
  neutrino losses, and screening cap as its value. The cap is still only the
  existing weak-screening approximation, not new dense-matter physics.

From the first complete boundary-value solve:

- **Newton convergence is not spatial accuracy.** The coarse polytrope can
  converge to a tiny residual while its radius differs by percent levels
  from the continuum solution. Keep the independent Lane–Emden mesh-refinement
  test and the check that the unresolved central sphere shrinks correctly.
- **Damping is not convergence.** A small accepted step may only reflect a
  line-search restriction. Require a small undamped remaining correction
  as well as a small residual. Physics-domain errors may reject a trial;
  they never authorize table extrapolation or an implicit fallback.

From connecting physical opacity to a stellar trial:

- **Changing a seed's density changes its enclosed mass.** Keeping the old
  Lane–Emden mass mesh after adjusting envelope density produced folded
  surface layers. Rebuild the mesh consistently and inspect monotonicity.
- **Stored points do not establish atmospheric coverage.** The earlier
  Ferguson trial had every stored point below log R=1, yet its atmosphere
  reached the edge. AESOPUS has now supplied that missing cool density data;
  OPAL's hot density limit remains when that source is selected. The new
  TOPS import supplies the .1 Msun hot dense coverage. Report the failing module explicitly.
- **A checked correction can need refinement.** Finer stellar meshes exposed
  a linear backward-error failure absent in the small benchmarks. Correcting
  the original linear residual fixes this without weakening acceptance.
  Independent virial integration then verifies spatial hydrostatic accuracy.

---

## 6. Hooks for massive white dwarfs (the Adams & Laughlin update)

Deliberately not built. Nothing about them should require rearranging what
exists. **The electrons are already exact at Chandrasekhar-mass densities** —
what such a star additionally needs is *ion* physics, which is exactly why the
EOS is a sum of `EosComponent`s:

| Need | Where it goes |
|---|---|
| Coulomb energy of the ion lattice | new `EosComponent` |
| Crystallisation, latent heat, Debye solid | new `EosComponent` |
| C/O phase separation on freezing | new `EosComponent` |
| Plasmon/pair neutrino losses | `NeutrinoLosses` (declared, null) |
| Electron conduction | `Conduction` + `CombinedOpacity` (declared) |
| C/O interiors | `Composition` already carries C12, O16 |

Noted but **not** designed: general relativity in the structure equations,
which the last percent of mass before the Chandrasekhar limit would want.

---

## 7. State of the FORTRAN line (context, not work)

`/Users/greglaughlin/Projects/low_mass_stars`, branch `main`, pushed through
`6689b6d`. This session fixed three walls there:

- `cd26677` — the φ-scaling wall (see §5 item 3)
- `1807d06` — adaptive fitting point (§5 item 4)
- `6689b6d` — Ferguson 2005 opacities replacing AJR83 (§5 item 5)

Solar calibration under Ferguson: **X = 0.737, α = 1.93 → L = 1.004,
R = 1.000, Teff = 5791**. (α rose because Ferguson is ~3× more opaque near
10,000 K.) Note these tracks are **no longer LBA97's** — `IOPC=2` preserves
the 1997 physics exactly for reproducing the paper.

Results: 0.30 and 0.35 M☉ complete their lives and **do not** ignite helium
(peak Tc 4.5e7 and 6.0e7). No star has reached 1e8 K, so on these models the
**flash mass is above 0.55 M☉** — still unmeasured.

**Known regression:** with Ferguson opacities the 0.10 M☉ flagship no longer
reaches the 1e-6 L☉ stop; it stalls at L = 5.2e-6, Teff = 1704 K. Open
question worth the user's judgement: *should* Ferguson's grain opacities apply
to a high-gravity helium-white-dwarf atmosphere at all?

There may still be a background 0.55 M☉ run (`fg055`) in the session
scratchpad; it was holding its Hayashi line correctly at Teff ≈ 3400 K.

---

## 8. Practical gotchas in this environment

- **The shell's cwd resets** between tool calls. Use absolute paths or `cd`
  inside every command.
- **`pkill -f "henyey77 <"` matches nothing** — the `<` is shell syntax, not
  part of the command line. Two processes once wrote the same output file and
  produced an interleaved, backwards-running track. Use `pkill henyey77`.
- The session scratchpad is **cleared between sessions**; anything worth
  keeping goes in the repo.
- Long runs: launch with `run_in_background`, watch with a `Monitor` whose
  filter matches **every** terminal state, not just success.
