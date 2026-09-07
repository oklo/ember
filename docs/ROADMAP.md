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

## Next

- [ ] Coulomb corrections and crystallisation (Potekhin & Chabrier)
- [ ] CMS19 / Chabrier-Debras tabulated H/He equation of state, blended to the
      analytic form outside the table
- [ ] Opacity: Ferguson 2005 (tables already assembled in the Fortran line),
      OPAL/OPLIB above, conduction from Cassisi 2007
- [ ] Nuclear: pp chains with modern rates, He3 tracked; CNO out of equilibrium
- [ ] Atmosphere: tabulated model atmospheres, with a grey fallback
- [ ] Henyey solver: block elimination in (ln r, ln rho, ln T, L)
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
