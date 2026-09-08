# Static stellar equilibrium

The current **0.1 Msun fixed-composition equilibrium** uses a consistent
FreeEOS-based Helmholtz potential, AESOPUS/TOPS radiative opacity, pp
heating, MLT, and an AMES-COND non-grey boundary at Rosseland tau=100.
At 4096 points it gives **R=.12973875 Rsun, L=.00094857566 Lsun,
Teff=2812.29 K**. This is an experimental static model with declared
composition and source-fit approximations, not an age-matched calibration
or an evolutionary calculation.

## Reproduce the current reference

```
cmake --build build
mkdir -p out
build/apps/ember-equilibrium 4096 .1 .15 --eos freeeos --hot-opacity tops --seed-index 1.5 --atmosphere cond-solar-proxy > out/equilibrium-m010-nongrey-reference.json
ctest --test-dir build --output-on-failure
```

Positional arguments are mesh points, mass/Msun and seed radius/Rsun;
radius is free to relax. The n=1.5 Lane–Emden seed is an initial guess,
not an equilibrium constraint. Its envelope uses the atmosphere's local
matching temperature, which exceeds Teff at tau=100. Enclosed mass is
rebuilt after adjusting the envelope; relaxation uses that fixed mass mesh.
Use 2048 points for a smaller working model.

The old driver defaults remain 128, .1, .2, ionized EOS, OPAL, n=3 and grey
atmosphere; that particular trial fails at the OPAL density edge. Use the
explicit command above. `--atmosphere grey` selects the radiative Eddington
alternative at tau=2/3; `--tau-top` then defaults to .001 for CMS19/FreeEOS
and 1e-6 for the analytic ionized EOS. COND rejects `--tau-top`.

Exit status zero requires convergence. stdout is JSON; stderr contains
short diagnostics. Failures retain their reason and last accepted profile
when available. Diagnostics include the residual, undamped correction,
Newton history, nuclear integral, independent virial check, EOS Maxwell
identity, and full profile. A JSON file alone never establishes convergence.
Compact current metadata are in
[results/equilibrium_m010_nongrey.json](results/equilibrium_m010_nongrey.json).

## Separate effects of EOS and atmosphere

These 4096-point models share mass, initial composition, opacity and
MLT alpha=1.9. Each is independently seeded and relaxed. The FreeEOS
models use the .0125-dex potential grid. The grey runs use tau_top=.001.

| EOS | Atmosphere | R/Rsun | L/Lsun | Teff (K) | Maximum local Maxwell defect |
|---|---|---:|---:|---:|---:|
| CMS19 | Grey, tau=2/3 | .12860617 | .0009756295 | 2844.57 | .24178 |
| CMS19 | COND, tau=100 | .12938578 | .0009370485 | 2807.53 | .24124 |
| FreeEOS potential | Grey, tau=2/3 | .12857762 | .0010069939 | 2867.48 | 4.44e-16 |
| FreeEOS potential | COND, tau=100 | .12973875 | .0009485757 | 2812.29 | 4.44e-16 |

The combined change increases radius by .88%, decreases luminosity by
2.77% and cools Teff by 32.28 K relative to CMS19/grey. These differences
include replacing the EOS physical model, its interpolation and the
atmospheric boundary; they are not a measurement of one source error.
The machine-precision identity measures internal consistency, not physical
accuracy. [FREEEOS.md](FREEEOS.md) gives independent caloric/transport
checks and the source-response discrepancies introduced near fit joins.
The original [CMS19 source audit](CMS19.md) and energy guard are unchanged.

Against BHAC15's .1 Msun, 5-Gyr point, the new model is approximately
4.6% larger and 8.8% brighter, with Teff 1.3 K higher. Its inner state is
Tc=4.58580e6 K and rhoc=350.601 g/cm³. The temperature agreement has
improved, while the radius offset has grown. Unmatched abundances,
atmosphere models and age prevent interpreting these as calibrated errors.
The previous observed-star mass runs are clearly retained as CMS19/grey
in [LITERATURE_COMPARISON.md](LITERATURE_COMPARISON.md).

## Current mesh and integral checks

| Points | R/Rsun | L/Lsun | Teff (K) | Absolute relative virial error |
|---:|---:|---:|---:|---:|
| 1024 | .12974188 | .0009485454 | 2812.23 | 1.025e-5 |
| 2048 | .12972454 | .0009492765 | 2812.96 | 2.560e-6 |
| 4096 | .12973875 | .0009485757 | 2812.29 | 6.396e-7 |

The 4096-point solve takes nine accepted updates, with residual 6.60e-13
and undamped correction 6.76e-12. The pp integral and surface luminosity
agree to about 2e-14 relative; profiles are ordered with positive heat
capacities and finite internal energy. The unresolved center contains
1.23e-10 of the stellar mass.

As before, virial closure is independently integrated from the profile:

```
3*integral(P/rho dm) - 4*pi*R^3*P_surface = integral(G*m/r dm).
```

Its error decreases by approximately four per mesh doubling. Radius
changes by less than .014% and luminosity by less than .078% between
successive fine meshes. Their changes remain nonmonotonic, so this is
fine-mesh agreement, not demonstrated second-order convergence of R/L.
The separate .025-to-.0125-dex EOS-grid comparison changes the 1024-point
radius/luminosity by about .0030%/.021%; source-fit joins prevent claiming
asymptotic EOS-grid convergence from those two resolutions.

## Remaining physical limitations

Composition is fixed at X=.7, Y=.28, Z=.02, He3=0. The EOS represents
metals as helium (effective Y=.3), opacity uses the actual fixed X/Z with
GS98 metals, and resolved nuclear metals use the AAG21 helper. The COND
GN93 solar mixture and helium abundance are not matched to the interior;
using it requires an explicit solar-mixture-proxy option. No parameters
were adjusted to fit another model or an observed star.

The EOS potential regularizes small FreeEOS source-fit discontinuities.
At the audited off-grid points, source differences reach 1–2% in some
responses despite much smaller pressure differences. These are sampled
approximation errors, not uniform physical uncertainty bounds. Exact
thermodynamic consistency alone does not remove that qualification.

The non-grey boundary supplies actual material pressure and temperature
at tau=100 from structures, but retains the thin-atmosphere approximation:
its radius is identified with the photospheric radius and omitted mass is
neglected. The interior opacity has no grains. MLT alpha=1.9 is uncalibrated;
conduction, mixing, composition advance and time-step control are absent.
Current TOPS opacity supports only X=.7,Z=.02. The new EOS supplies a
validated fixed-composition caloric state, but He3/composition changes
remain unsupported, so it is not yet a burning-evolution model.

Next work should match the physical mixtures, improve source-fit joins and
qualify modern atmosphere/geometry effects, then add composition responses,
mixing, adaptive mesh and time-step control. These are still required for
the end-to-end low-mass evolutionary milestone.

## Verification and provenance

Fifteen CTest suites pass. `helmholtz` checks independent differences of
E/S/P, transport responses, PT inversion, source values at off-grid states,
C2 continuity, and static/time-dependent zone Jacobians. `atmosphere`
checks original COND values, derivatives, bounds and explicit mixture
selection. `nongrey_equilibrium` solves 1024/2048/4096-point stars, checking
profiles, caloric stability, nuclear balance, mesh agreement and virial
convergence. CMS19, opacity, the older stellar cases and independent
Lane–Emden continuum tests remain in place.

All 180,865 direct FreeEOS evaluations were reproduced byte for byte from
a fresh source build on this host. Both new runtime tables were reproduced
byte for byte from their versioned numerical inputs. Normal C++ builds need
neither Fortran nor network. Source versions, checksums, conversion details
and limitations are in [FREEEOS.md](FREEEOS.md),
[EOS provenance](../data/eos/README.md) and
[atmosphere provenance](../data/atmosphere/README.md).

## Retained CMS19/grey mesh experiment (2026-09-07)

These older CMS19/grey runs start independently from a .15 Rsun, n=1.5 seed with identical
physics, tau_top=.001 and surface matching at tau=2/3.

| Points | R/Rsun | L/Lsun | Teff (K) | Scaled residual | Absolute relative virial error |
|---:|---:|---:|---:|---:|---:|
| 128 | .12911938 | .0009629408 | 2829.64 | 5.31e-10 | 6.856e-4 |
| 256 | .12935265 | .0009420455 | 2811.62 | 1.67e-10 | 1.699e-4 |
| 512 | .12885921 | .0009637381 | 2833.08 | 1.04e-10 | 4.229e-5 |
| 1024 | .12862818 | .0009746770 | 2843.63 | 1.90e-11 | 1.055e-5 |
| 2048 | .12861667 | .0009751461 | 2844.10 | 1.68e-10 | 2.635e-6 |
| 4096 | .12860617 | .0009756295 | 2844.57 | 4.10e-10 | 6.585e-7 |

The 4096-point result has Tc=4.56631e6 K, rhoc=366.875 g/cm³ and
R=8.94713e9 cm. It converges in 12 accepted updates with an undamped
correction of 3.06e-10. Surface L=3.73471e30 erg/s agrees with the pp
heating integral to about 2e-15 relative. Radius increases, density and
temperature decrease, and luminosity stays positive.

The virial diagnostic is independently integrated from the stored points:

```
3*integral(P/rho dm) - 4*pi*R^3*P_surface = integral(G*m/r dm).
```

Both integrals use trapezoidal mass quadrature with a uniform unresolved
central sphere. The relative error drops approximately fourfold with each
mesh doubling. Nonlinear tolerances remain residual 1e-9 and undamped
correction 1e-8; linear backward-error acceptance remains 1e-10. No table
bounds or convergence requirements were loosened.

Coarse R and L are **not monotonic with resolution**. Sharp envelope/EOS
features need substantially more than 256 points. The 1024-to-2048 and
2048-to-4096 luminosity changes are .0481% and .0496%; the corresponding
radius changes are .00895% and .00817%. Thus fine models agree to about
.1% in L and .02% in R, but R/L are not yet in demonstrated clean
second-order convergence. Adaptive envelope resolution remains worthwhile;
a small virial error alone does not guarantee the thermal profile's accuracy.

Doubling tau_top from .001 to .002 changes the independent 256-point
seed-and-solve result by .000402% in radius and .00259% in luminosity.
The regression also repeats this perturbation on the **identical**
2048-point converged mass mesh, bounding R/L changes below .02%.
This assesses the finite starting column, not the accuracy of a grey
radiative T(tau) law for this atmosphere.

## Retained ionized-EOS benchmark at 0.5 Msun

Reproduce with `build/apps/ember-equilibrium 512 .5 .6`. This uses the
ionized EOS, AESOPUS/OPAL and n=3 seed. All listed runs start independently
from a .6 Rsun seed.

| Points | Updates | R/Rsun | L/Lsun | Teff (K) | Scaled residual | Relative virial error |
|---|---:|---:|---:|---:|---:|---:|
| 128 | 18 | 0.94635323 | 0.019650339 | 2221.4802 | 5.49e-11 | 1.179e-3 |
| 256 | 18 | 0.94342814 | 0.019831831 | 2230.0411 | 8.47e-13 | 2.927e-4 |
| 512 | 19 | 0.94219963 | 0.019877124 | 2232.7675 | 1.34e-10 | 7.291e-5 |

At 512 points the inner temperature is 8.24398e6 K and density 88.8026 g/cm³.
Radius increases outward, density and temperature decrease, and luminosity
is positive. Surface luminosity equals the trapezoidal pp heating integral
to about 1e-14 relative. Radius changes by 0.130% and luminosity by 0.228%
from 256 to 512 points. The independently integrated virial theorem,
`3*integral(P/rho dm) - 4*pi*R^3*P_surface = integral(G*m/r dm)`, closes
with an error decreasing by about four per mesh doubling. This checks
hydrostatic spatial accuracy independently of the nonlinear residual.

The finer initial meshes originally failed the Henyey original-equation
check. Up to three iterative-refinement solves now correct `J*dy+f` using
the same coefficients. The original componentwise backward-error threshold
of 1e-10 remains unchanged, as do nonlinear tolerances of 1e-9 residual
and 1e-8 undamped correction. A rejected correction reports its measured
backward error. No table bounds or convergence requirements were loosened.

These results establish a converged solution of the implemented equations.
They do not establish realistic M-dwarf radii or temperatures. The ionized
EOS and radiative grey atmosphere are substantial envelope approximations;
the very cool, inflated result is a reason to replace those approximations
before making observational comparisons or evolving the model.
