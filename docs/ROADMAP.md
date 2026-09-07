# Roadmap

The target that defines "working": a 0.1 solar-mass star evolved end to end,
from the Hayashi track through twelve trillion years of hydrogen burning, the
blueward turn, and down the helium-white-dwarf cooling track to below
10^-6 Lsun — in one run, with no hand-holding.

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
- [x] Atmosphere table reader and surface residuals with an analytic Jacobian;
      physical atmosphere data remain pending
- [x] Analytic zone Jacobian, including EOS transport responses and screened
      mass-defect heating derivatives; numerical assembly retained for tests
- [x] Regular central boundary conditions, pivoted Henyey block elimination,
      and damped Newton relaxation on a fixed mesh
- [x] Complete radiative-polytrope benchmark against an independent
      Lane–Emden solution, including spatial convergence and failure checks

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
| Plasmon/pair neutrino losses | `NeutrinoLosses`, declared, currently null |
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
- [ ] CMS19 / Chabrier-Debras tabulated H/He equation of state, blended to the
      analytic form outside the table
- [ ] Opacity: Ferguson 2005 (tables already assembled in the Fortran line),
      OPAL/OPLIB above, conduction from Cassisi 2007
- [ ] CNO out of equilibrium (needed on the pre-main-sequence, where it runs
      once and never again); Chugunov screening for dense matter
- [ ] Import and validate physical model-atmosphere grids (BT-Settl archive
      timed out; the existing atmosphere table is synthetic test data only)
- [x] Model on a Lagrangian mass mesh in (ln r, ln rho, ln T, L)
- [x] The four structure equations as zone residuals, with a Jacobian checked
      against the residual it differentiates
- [ ] Physical stellar equilibrium with actual pp heating and an explicitly
      selected atmosphere; requires high-temperature opacity coverage
- [ ] Adaptive mesh (the "temporary points" idea, done properly)
- [ ] Time stepping with error control rather than iteration-count heuristics
- [ ] 0.1 Msun end to end

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
