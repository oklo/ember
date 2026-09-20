# Roadmap

The active objective is physically supported evolution of an initially
0.1-solar-mass star through hydrogen exhaustion, helium-remnant cooling past
**Teff=100 K**, and disappearance under explicit proton/nucleon-decay scenarios.
Formation belongs in the full lifetime account; the current numerical clock
starts from a static main-sequence model. Nuclear and cooling ages are outputs,
while unmeasured particle lifetimes and future environments are assumptions.
Compare those assumptions and retain a stable-baryon alternative.

**September 11 status:** the latest computed model reaches **3.890 trillion
years**, with central **X=0.0008157**, at an interior opacity density limit.
Its full endpoint audit is pending; the checked continuation reaches
**3.881 trillion years**. Dense opacity is the immediate constraint. See
[COLD_REMNANT.md](COLD_REMNANT.md) for the current acceptance order,
[LIFETIME_SURVEY.md](LIFETIME_SURVEY.md) for the wider programme, and the
[September 11 working paper](reports/2026-09-11/ember_status_and_future.pdf).
The checklist below retains architecture history with current qualifications.

## Done

- [x] Repository, C++23 build (CMake + Ninja), host-tuned on Apple silicon
- [x] Constants (CODATA 2018, IAU 2015 nominal solar values)
- [x] Composition: 8 species, AAG21 solar mixture, scaled to (X, Z)
- [x] Equation of state interface with analytic derivatives
- [x] Ideal + radiation + non-relativistic Fermi-Dirac electrons
- [x] Limit tests: ideal, radiation-dominated, Chandrasekhar, Maxwell relation

- [x] Relativistic Fermi-Dirac, panelled at the Fermi surface
- [x] Limit tests: ultra-relativistic (4/3 law) and the mildly relativistic
      white-dwarf core between the two laws
- [x] Ferguson (2005) low-temperature opacities, monotone-interpolated, with
      table edges as errors rather than extrapolation
- [x] Composable equation of state, proven equal to the monolithic one across
      five regimes
- [x] Hooks for cold, massive remnants (see below)
- [x] pp chains with He3 followed explicitly, energy derived from the mass
      defect of the code's own nuclide masses rather than a separate Q table
- [x] Mixing-length convection and the Schwarzschild criterion, with a
      bounded cubic solve, analytic gradient partials, and a regression for
      the transport-row conditioning at grad/grad_rad below 1e-6
- [x] Grey atmosphere with varying opacity, radiation pressure, adaptive
      integration, and analytic sensitivities to Teff and gravity
- [x] Atmosphere table reader and surface residuals with an analytic Jacobian
- [x] Analytic zone Jacobian, including EOS transport responses and screened
      mass-defect heating derivatives; numerical assembly retained for tests
- [x] Regular central boundary conditions, pivoted Henyey block elimination,
      and damped Newton relaxation on a fixed mesh
- [x] Complete radiative-polytrope benchmark against an independent
      Lane–Emden solution, including spatial convergence and failure checks
- [x] OPAL GS98 fixed-Z rectangular subset and smooth Ferguson blend, with
      derivative checks and strict density/composition coverage
- [x] AESOPUS 2.1 GS98 gas opacity to log R=6, original cells and strict bounds
- [x] Experimental static stellar driver with structured success/failure diagnostics
- [x] Fixed-composition 0.5 Msun equilibrium with the current approximate EOS
      and grey atmosphere; three-mesh convergence, energy and virial checks
- [x] Henyey iterative refinement retaining the original backward-error threshold
- [x] CMS19 original H/He pressure/entropy cells, actual interpolation Hessians,
      strict fluid masks and explicit static-only energy guard
- [x] LANL TOPS hot dense opacity at X=.7, Z=.02; native-density interpolation,
      excluding every server-substituted density
- [x] Experimental 0.1 Msun static equilibrium through 4096 points; fine-mesh
      agreement, nuclear balance, independent virial and atmosphere-depth checks
- [x] Consistent fixed-composition C2 Helmholtz EOS from FreeEOS 3.0 EOS1;
      caloric and transport derivatives independently tested; source-fit join
      regularization and sampled accuracy explicitly audited in `FREEEOS.md`
- [x] Physical AMES-COND non-grey tau=100 boundary via untouched MESA cells,
      explicit solar-mixture proxy, and a combined 4096-point stellar solve

- [x] Composition-dependent Helmholtz potential family, H1/He3 responses,
      explicit baryonic abundance convention; metals/isotope approximations qualified
- [x] Composition-dependent TOPS source planes, strict common density support
- [x] Coupled implicit pp burning, conservative instantaneous convective mixing,
      thermal first law and adaptive step-doubling; bounded 10-Gyr experiment

## Beyond the first target: massive white dwarfs

The 0.1 Msun star is the milestone, but the same machinery is wanted for the
ultra-cold evolution of heavier remnants - a solar remnant, and structures
approaching the Chandrasekhar mass - in support of an update to Adams &
Laughlin (1997). Those are deliberately *not* built yet, but nothing about
them should require rearranging what is here. The interfaces they need already
exist:

| What a cold massive remnant needs | Where it goes |
|---|---|
| Coulomb energy of the ion lattice | a new `EosComponent` |
| Crystallisation, latent heat, Debye solid | a new `EosComponent` |
| C/O phase separation on freezing | a new `EosComponent` |
| Plasmon/pair neutrino losses | `NeutrinoLosses`; HRW plasma implemented, other channels pending |
| Electron conduction (Cassisi 2007, Blouin 2020) | `Conduction` + `CombinedOpacity` |
| Carbon/oxygen interiors | `Composition` already carries C12 and O16 |
| Arbitrary relativity at high density | already done; the electrons are exact |

The electron gas is already adequate for a Chandrasekhar-mass star: what such
an object additionally requires is the *ion* physics, which is why the
equation of state is a sum of components and not one closed expression.
Anything needing general relativity in the structure equations (the last
percent of mass before the limit) would be a new `Structure` term and is noted
but not designed.

## Next

- [ ] Coulomb corrections and crystallisation (Potekhin & Chabrier)
- [x] GS98 metal potential and refined composition coverage through zero H,
      within the original thermal domain; cold dense EOS still pending
- [x] Hydrogen-poor TOPS refinements, X/Z source blends and Cassisi-family
      conduction; source and mixture approximations remain under assessment
- [ ] CNO out of equilibrium (needed on the pre-main-sequence, where it runs
      once and never again); Chugunov screening for dense matter
- [ ] Match modern model-atmosphere composition to the interior; quantify
      the thin-atmosphere approximation and atmospheric-grid interpolation
- [x] Model on a Lagrangian mass mesh in (ln r, ln rho, ln T, L)
- [x] The four structure equations as zone residuals, with a Jacobian checked
      against the residual it differentiates
- [ ] Quantitatively credible 0.1 Msun equilibrium: reduce mixture/source-fit
      approximations and repeat age/composition-aware stellar and observational
      comparisons. Current numerical and bulk checks are in `EQUILIBRIUM.md`
- [ ] Adaptive mesh (the "temporary points" idea, done properly)
- [x] Time stepping with error control from one full versus two half steps
- [x] Ledoux buoyancy and semiconvective/thermohaline composition transport
- [ ] Microscopic diffusion/settling, its energetics, and converged moving boundaries
- [ ] Hydrogen exhaustion and helium-remnant cooling through 100 K
- [ ] Quantum solid and partially ionized/neutral envelope physics
- [ ] Coupled grains, element depletion and settling in cold atmospheres
- [ ] Residual and pycnuclear burning with process-specific rate bounds
- [ ] Gas, interstellar-object and planetary-debris accretion and retained heat
- [ ] Heating from WD encounters that remain unbound, including crystallized stars
- [ ] Evolving galactic environments, ejection and stochastic heating histories
- [ ] Conditional dark-matter capture, thermalization, annihilation and evaporation
- [ ] Proton/nucleon-decay channels, nuclear daughter products and energy deposition
- [ ] Gravitational decay through virtual black holes; distinguish collective tunneling
- [ ] Decay-driven mass loss, loss of degeneracy, neutral matter and optical thinning
- [ ] Dispersal and the final discrete-particle regime with an explicit end criterion
- [ ] Full formation-to-disappearance timeline in the working paper

These investigations are specified in [ULTRACOLD_PHYSICS.md](ULTRACOLD_PHYSICS.md).
The [approximation remit](APPROXIMATION_REMIT.md) determines when detailed
structure, reduced models or conditional estimates are appropriate. The
[Fable assignment](research/fable/ASSIGNMENT.md) covers the independent
encounter survey and literature work.
Their importance is set by the calculated structure and luminosity, not a
universal 100 K switch. Environmental and particle assumptions may alter or
prevent a particular temperature crossing.

## Lessons carried over from the Fortran line

These cost real time there and the design here is meant to make them
impossible:

1. **Fixed-format input.** A six-digit number in a five-column card field
   silently ate a flag and produced a physically wrong run that looked fine.
   Configuration here is typed and validated, never column-counted.
2. **Fixed-width output.** A model counter in an `I4` field printed `****`
   past 9999 and corrupted downstream analysis. Output here is structured.
3. **Equation forms valid only in a limit.** Scaling the radiative equation by
   the convective efficiency is correct only where that efficiency is order
   unity; at 1e-6 it decoupled the temperature from the luminosity and killed
   every giant. Where a formulation has a domain, the code should assert it.
4. **Fixed structural constants.** A fitting point at a fixed mass fraction
   drifts above the photosphere when a star swells to 300 solar radii. Anything
   that is "a good value for the Sun" is suspect.
5. **Extrapolated tables.** Below its floor the 1983 opacity table returned
   kappa ~ 1e-9 and the envelope became transparent, removing the Hayashi
   limit entirely. Table edges must be errors, not silent extrapolation.
