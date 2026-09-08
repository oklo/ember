# Static stellar equilibrium

An **experimental 0.1 Msun fixed-composition equilibrium now converges**
with CMS19 H/He pressure/entropy, AESOPUS/TOPS radiative opacity, pp heating,
MLT and an Eddington grey atmosphere. It passes fine-mesh agreement, nuclear
energy balance and an independent virial check. Its EOS has a substantial
local thermodynamic consistency defect, and its atmosphere is still grey.
**This is a numerical milestone toward a physical stellar model, not an
observational calibration or an evolutionary calculation.**

A [comparison with published main-sequence grids](LITERATURE_COMPARISON.md)
finds R about 3.7% larger, L about 12% brighter and Teff about 34 K hotter
than BHAC15 at .1 Msun and 5 Gyr. This is a comparison of bulk properties,
with different composition and atmosphere physics, not a calibration.

## Reproduce the 0.1 Msun reference

```
cmake --build build
mkdir -p out
build/apps/ember-equilibrium 4096 0.1 0.15 --eos cms19 --hot-opacity tops --seed-index 1.5 > out/equilibrium-m010-reference.json
ctest --test-dir build --output-on-failure
```

The 2048-point run is also a useful working model. Positional arguments are
mesh points, mass/Msun and seed radius/Rsun; radius is free to relax. The
flags select the EOS, hot opacity, polytropic seed index and optional
`--tau-top`. CMS19 defaults to tau_top=.001; ionized EOS defaults to 1e-6.
The older driver defaults remain 128, .1, .2, ionized EOS, OPAL and n=3;
that particular trial still fails at the OPAL density edge. Use the explicit
command above for the new model. `--help` describes the interface.

Exit status is zero only on convergence. stdout is JSON, stderr contains
short diagnostics. A failed relaxation returns its last accepted profile
and failure reason. Invalid input or initial physics returns a JSON failure
without a profile. Output includes the residual and undamped correction,
Newton history, pp heating integral, independent virial error, EOS Maxwell
defect and its location, and the mass/radius/density/temperature/luminosity
profile. A JSON file's existence alone never establishes convergence.

Composition is held fixed at X=.7, Y=.28, Z=.02, He3=0. This is equilibrium
at the specified composition, not chemical equilibrium or a star at a
specified age. The CMS19 EOS explicitly approximates metals as helium
(effective Y=.3); GS98 opacity keeps the actual X and Z. Resolved nuclear
metals use the AAG21 helper. MLT alpha=1.9 is uncalibrated. No mixing,
composition advancement, conduction or grain opacity is included.

The n=1.5 Lane–Emden seed supplies a convective initial guess, with pp
luminosity and outer temperature/density adjusted toward the atmosphere.
Enclosed mass is then rebuilt from the changed density. Relaxation solves
on that fixed mass mesh. The seed is not an equilibrium constraint. The
unresolved central sphere contains about 1.23e-10 of the total mass.

## Mesh and integral checks (2026-09-07)

All runs start independently from a .15 Rsun, n=1.5 seed with identical
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

## What remains physically unresolved

CMS19 removes the fully ionized envelope approximation, but pressure and
entropy in this source version are not an exact thermodynamic pair. Define
D=P*delta/(rho*T*cp*grad_ad)-1. At 4096 points, max |D|=.24178 near
T=27304 K, rho=.12394 g/cm³ and m/M=.99999310; the mass-weighted RMS is
.0043487. The small mass in those layers does not bound their effect on
stellar radius or luminosity. An original pure-H source node already has
a 26.5% defect without any ember interpolation. This limitation is visible
in JSON and explicitly tested; [CMS19.md](CMS19.md) gives the audit.

The source internal-energy columns have an independent join defect. Energy
is unavailable in this implementation, and **positive-dt structure calls
reject it**. Thermodynamically consistent caloric physics and support for
He3/composition changes are required before attempting evolution.

TOPS supplies actual hot dense opacity where OPAL ended at log R=1. Its
server response included some substituted densities; the importer excludes
every one of those cells, retaining two complete original rectangles.
Native log-rho interpolation and smooth blends preserve strict support.
This TOPS request has only X=.7, Z=.02, so it cannot follow fuel depletion.
See [opacity provenance](../data/opacity/README.md).

Next work should resolve EOS pressure/entropy/energy consistency and import
physical atmosphere structures, then reassess the equilibrium before adding
composition mixing, adaptive mesh and time-step control. A credible surface
boundary needs pressure-temperature structures, not spectrum files. No
physical atmosphere grid has yet been imported; the table reader's fixture
is explicitly synthetic.

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

## Verification

Thirteen CTest suites pass. `cms19` checks independent original pure and
mixture cells, molecular/neutral/ionizing limits, actual interpolation
responses, full zone Jacobians, strict domain failures, finite atmosphere
sensitivity and rejection of unvalidated time-dependent energy. It keeps
the known dense ionization consistency defect visible.

`tops` checks all 5706 imported rectangular cells, independent native
values, both native-density and blend derivatives, dilute scattering,
fixed composition and exclusion of substituted-density regions.
`low_mass_equilibrium` independently solves 512/1024/2048 meshes, checks
R/L agreement, ordered profiles, pp balance, virial convergence and the
atmosphere starting-depth perturbation. The older 0.5 Msun integration test
and independent radiative Lane–Emden continuum benchmark remain in place.

All six newly imported EOS/opacity files were reproduced byte for byte from
their checksum-pinned originals. CLI invalid-input and table-domain failures
were checked for nonzero exit status and valid JSON. The compact reference
metadata are retained in [results/equilibrium_m010.json](results/equilibrium_m010.json);
the command above regenerates the full profile and Newton history.
