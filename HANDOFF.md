# ember — handoff

Written 2026-09-07 at the end of a long session. Read this first; it is meant
to be the only thing you need.

Updated 2026-09-07 after implementing central boundary conditions, Henyey
block elimination, and complete-model relaxation against a Lane–Emden
benchmark. The latest checkpoint is described below; use git for its hash.

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

Eight test suites, all passing (EOS, opacity, nuclear, structure, convection,
atmosphere, Henyey, relaxation). A complete **controlled radiative polytrope**
now converges; a physical stellar equilibrium and an evolutionary run remain
pending. The benchmark has constant heating and opacity and an explicitly
artificial atmosphere. It must not be reported as the 0.1 M☉ scientific run.

| Module | File | Status |
|---|---|---|
| Constants | `include/ember/constants.hpp` | CODATA 2018, IAU 2015 nominal solar |
| Composition | `composition.{hpp,cpp}` | 8 species, AAG21 mixture, scaled to (X,Z) |
| EOS interface | `eos.hpp` | `EosResponse` includes cp, delta, grad_ad derivatives; missing response implementations throw; `rho_from_PT` reports failure |
| EOS components | `eos_component.hpp`, `eos_components.cpp` | ions, radiation, relativistic FD electrons, with component Hessians |
| Composite EOS | `eos_composite.hpp` | additive; == monolithic to 1e-12 over 5 regimes |
| Fermi–Dirac | `src/fermi.{hpp,cpp}` | relativistic, centered density-constrained Hessians; shell split at the Fermi surface |
| Opacity | `opacity.hpp`, `opacity_ferguson.{hpp,cpp}` | Ferguson 2005, monotone, edges throw |
| Interpolation | `interp.{hpp,cpp}` | monotone Hermite; parameter derivatives follow active limiter branches |
| Nuclear | `nuclear.hpp`, `nuclear_pp.cpp` | pp chains, He3 explicit; heating derivatives include the same mass defect, neutrinos, and screening as the value |
| Model | `model.hpp` | (ln r, ln ρ, ln T, L); solver mesh starts at positive m[0] and ends at M |
| Structure | `structure.{hpp,cpp}` | 4 residuals + analytic Jacobian; value-only path and numerical Jacobian retained for tests |
| Convection | `convection.{hpp,cpp}` | BV58 MLT, Schwarzschild criterion, analytic gradient partials; wired into transport |
| Atmosphere | `atmosphere.hpp`, `atmosphere_grey.cpp` | Eddington grey, varying opacity, adaptive integration and analytic sensitivities |
| Atmosphere tables | `atmosphere_table.{hpp,cpp}` | strict reader and bilinear interpolation; **no production grid yet** |
| Surface boundary | `boundary.{hpp,cpp}` | two surface residuals and analytic 2x4 Jacobian, assembled in relaxation |
| Central boundary | `boundary.{hpp,cpp}` | regular unresolved central sphere: r³=3m/(4πρ), L=m*eps_total; analytic derivatives |
| Henyey | `henyey.{hpp,cpp}` | pivoted block elimination, linear work/storage; checks original linear equations |
| Relaxation | `relaxation.{hpp,cpp}` | damped Newton, fixed mesh/composition, explicit failure, residual AND correction convergence |
| Continuum benchmark | `examples/radiative_polytrope.hpp`, `apps/polytrope.cpp` | independent n=3 Lane–Emden reference; JSON output via ember-polytrope |
| Losses | `losses.hpp` | **declared, null** — plasmon neutrinos for massive WDs |
| Conduction | `conduction.hpp` | **declared, unimplemented** — Cassisi 2007 |

Data: `data/opacity/ferguson_gs98_z020.dat` (versioned with the code
deliberately — a result is reproducible only if its numbers travel with it).

Benchmark: `build/apps/ember-polytrope 256 > out/polytrope-256.json` after
`mkdir -p out`. At 64/128/256 points, radius errors are 3.140%/0.8147%/0.2049%,
with 6/7/7 accepted Newton updates. This checks second-order spatial
convergence independently of the tiny nonlinear residual. See
`docs/HENYEY.md` for the equations, units, failure behavior, and limitations.

---

## 4. Immediate next steps, in order

**Next active dependency:** high-temperature opacity from item 8, to permit
an actual stellar interior with the analytic EOS and pp heating. The Henyey
machinery is now tested, but Ferguson alone cannot reach interior
temperatures. A physical equilibrium and a suitable initial model are the
next scientific solver validation; the polytrope is only controlled physics.

1. **Done: mixing-length convection** in the transport equation, plus the
   Schwarzschild criterion. The bounded cubic retains small gradient
   differences in both limits; the transport row remains in plain gradient
   form. Tests cover flux conservation, element cooling, gradient derivatives,
   and row conditioning with ∇/∇_rad < 1e-6. See `docs/CONVECTION.md` for the
   equations and coefficient convention. This is optically thick interior
   MLT; composition mixing and optically thin losses remain separate work.
2. **Atmosphere boundary condition.** Tabulated model atmospheres
   (PHOENIX/BT-Settl for M dwarfs and BDs) still need **physical data imported**.
   Grey fallback, table reader, and surface residuals are implemented. The
   BT-Settl AGSS2009 STRUCTURES archive timed out over HTTP/HTTPS this session;
   spectrum files are not pressure-temperature structures. The only added
   atmosphere table is a labeled synthetic test fixture under `tests/data`.
   See `docs/ATMOSPHERE.md` for format, sources, and integration conventions.
   The grey fallback uses radiative T(tau), not a convective atmosphere.
3. **Done: analytic Jacobian assembly.** `EosResponse` now supplies the
   transport derivatives from component Hessians. Zone assembly propagates
   module derivatives through all four equations; `zone_equations` and
   `zone_residual_numerical` remain independent of the transport-response
   path for tests. Checks cover contraction, degenerate helium, finite and
   efficient convection, and zero/inward luminosity. See `docs/JACOBIAN.md`.
   The old left-endpoint gravitational-energy difference remains; its
   discretization order is not solved by differentiating it. Positive `dt`
   now requires a previous model on the identical mass mesh.
4. **Done: Henyey block elimination and central boundary conditions.**
   Six active rows eliminate four variables with row pivoting; two
   constraints propagate outward and close with the surface equations.
   `relax` uses fixed physical units and a backtracking line search. It
   requires both residual and undamped correction convergence, and returns
   the last accepted model with an explicit failure status on exhaustion.
   Positive `dt` supports a fixed previous model, but no age/composition
   advance or time-step controller is implemented. A small positive inner
   mass replaces the impossible `ln r` at m=0; check that sphere's size.
5. **Adaptive mesh** — refine at burning shells, coarsen in isothermal cores.
   Must carry a *density* term in the smoothness measure (see §5 item 4).
6. **Time-step control** by estimated error, not iteration-count heuristics.
7. **CMS19 / Chabrier–Debras tabulated H/He EOS**, blended to the analytic
   form outside the table.
8. **Conduction** (Cassisi 2007) + `CombinedOpacity`; **OPAL/OPLIB** above
   31,600 K where Ferguson stops.
9. **CNO out of equilibrium** (runs once on the pre-MS and never again).
10. 0.1 M☉ end to end. GitHub pushes are authorized before this milestone.

The default α = 1.9 is **not calibrated for ember**. Do not import the
Fortran solar calibration as if the EOS and atmosphere were identical.

Physical equilibrium work can use the grey fallback while the physical
atmosphere grid is sourced, once opacity spans the interior. Do not mark
tabulated atmosphere physics complete merely because its reader passes
synthetic tests.

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
     "mass fractions sum to zero" — they must *not*; they fall at
     Σ dX/dt = −(ε+ε_ν)/c², and testing against zero hides the error)
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
