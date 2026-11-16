# Helium-rich non-grey atmospheres

The continuous lifetime calculation combines composition-dependent gas grids
at an interior matching optical depth of 100. Its selected files and physical
limits are described in [the common driver](LIFETIME_DRIVER.md) and
[the hydrogen-rich atmosphere inputs](../data/atmosphere/lifetime_hydrogen_envelope/README.md).
This document describes the helium-rich source calculations and the common
source-generation and acceptance checks. Gas-opacity coverage does not by
itself supply the full atmosphere physics needed for white-dwarf cooling.
See [the current trajectory](COLD_REMNANT.md) and [EOS data](../data/eos/README.md).

`CompositionAtmosphereGrid` supplies an atmosphere boundary from independent
LTE radiative-transfer calculations on a rectangular grid in baryonic H1,
He3, effective temperature and surface gravity. The five Ember metal
abundances and Rosseland matching depth are fixed for a family. This backend
is distinct from the differential grey correction described in
[ATMOSPHERE.md](ATMOSPHERE.md).

The source pipeline uses [TLUSTY208/SYNSPEC54](https://arxiv.org/abs/2104.02829)
and the molecular data distributed with
[Synple](https://github.com/callendeprieto/synple/tree/703549b60031c5f6dd2446ad933b129c070a6560).
It solves hydrostatic balance, wavelength-dependent radiative transfer and
radiative/convective equilibrium. SYNSPEC supplies absorption per unit
baryonic mass. TLUSTY adds electron and Rayleigh scattering; they must not
also be included in the absorption table.

## Physical specification

The recipe includes atomic lines, the general solar molecular line list,
TiO and water, the H-, H2+, He-, CH, OH and H2- continua, and collision-induced
absorption from H2-H2, H2-He, H2-H and H-He. Molecular equilibrium is enabled
below 10,000 K. The source's partition functions restrict relevant molecular
temperatures to 1000–16,000 K; the material grid stays inside that
interval. Wavelengths are in vacuum, spanning 900 Å to 30 microns. The adopted
microturbulence is 1 km/s and mixing-length parameter is 1.9.

The distributed CIA tables have separate temperature limits: H2-H2 and
H2-He span 1000–7000 K, H2-H spans 1000–2500 K, and H-He spans
1000–10,000 K. The source holds coefficients constant above their upper
temperature endpoints and supplies negligible absorption outside their
wavenumber intervals. These retained source prescriptions are additional
approximations, particularly for H2-H in dissociation layers. The source
H2- continuum spans 3505–1.519e5 Å; it contributes zero outside that range.

The atmosphere is a plane-parallel gas in LTE, with no irradiation,
condensates, rainout, magnetic support or spots. These choices are physical
approximations. In particular, cool upper layers can form condensates even
when Teff exceeds their condensation temperature. This gas calculation must
not be advertised as reproducing the dust chemistry of the COND anchor.
The MLT parameter is uncalibrated; equal alpha values in two different MLT
implementations do not establish identical convective efficiency.

Metal opacity uses the same declared GS98 proxy as the extended interior
opacity family. It does not resolve Ember's five-species metal inventory
element by element. The atmosphere EOS includes molecular metal chemistry;
the selected interior FreeEOS family uses the same GS98 element pattern,
with a documented trace potassium omission. Differences between the two
thermodynamic treatments are measured at the atmosphere matching states.

Helium-3 is important in this star. A source calculation using He4 masses
throughout would give the wrong number of helium atoms per gram. Both source
programs therefore load an explicit element-mass file. With atomic mass unit
`mu`, H has mass `mu`, and the mean helium mass in units of `mu` is

```
AHe = (X3 + X4) / (X3/3 + X4/4).
NH/rho = XH/mu
NHe/rho = (X3/3 + X4/4)/mu.
```

Other elements use declared representative integer isotope masses. The
patched source weights multiply its native `HMASS`, so the density used for
opacity, the molecular EOS and hydrostatic balance is consistently baryonic.
There is no gravity rescaling. Electronic cross sections, isotope line
shifts and collision-induced absorption remain approximated by the source
helium data; an effective helium mass in Doppler widths is not an
isotope-resolved line profile. Source physical constants are retained for
reproducibility; their Stefan–Boltzmann constant differs from Ember's by
about 0.022%.

## Numerical balance depth

The optional specification field `energy_balance_tau_division` sets TLUSTY's
`TAUDIV`: the depth below which differential flux balance supplements local
radiative equilibrium. It defaults to 0.01, preserving existing input files.
The value must lie between the top and matching optical depths. The importer
checks it against the saved input, and changing it changes the run identity.
The correction, flux, chemistry and matching-depth acceptance checks are unchanged.

At 4800 K and log g = 6.2, and at 4700 K and log g = 6.0, setting this depth to
0.0001 allowed previously oscillating pure-hydrogen columns to converge.
The independent 4800 K, log g = 6.0 control changed matching temperature by
0.0001852% and gas pressure by 0.001634%.
[Source checks and assumptions](https://github.com/oklo/ember/blob/e8f27e4f1eefdac7f4a25aa3c22a642475afb705/docs/results/highg_energy_balance_sept27_v1.json).
This is a numerical control for source calculations, not an automatic change
to a running star's atmosphere.

## Reproduction

Normal C++ builds require neither a Fortran compiler nor network access.
Generating a source family requires Python 3.11 or later, gfortran, several
GB of temporary storage, and a substantial calculation. The preparation
step retains and verifies pinned downloads and records compiler, executable,
data and patch hashes. The public FTP login in the script is the one
published by the source author for these line lists.

```
python3 scripts/prepare_nongrey_sources.py /tmp/ember-atmosphere-source
python3 scripts/generate_nongrey_grid.py \
  /tmp/ember-atmosphere-source/prepared.json \
  data/atmosphere/sources/nongrey_candidate_specification.json \
  /tmp/ember-atmosphere-family --jobs 3 \
  --initial-models data/atmosphere/sources/nongrey_gs98_z020
```

`--offline` on preparation uses the verified download cache. Generation
supports `--opacity-only` and `--plane N`. Source runs are reused only when
their input/executable fingerprints and output hashes match. Changed
specifications or provenance require a new work directory. A Fortran zero
exit code alone is insufficient: `STOP` can return zero after a failure.
Each material isotherm runs in a fresh source process and is checkpointed
independently. Assembly copies the original cells into the full source
table without interpolation. Completed isotherms survive a failed later
calculation; no incomplete atmosphere family is imported.
The 21-temperature material grid extends down to 1076 K. Its explicitly
stored temperature values preserve the earlier 19 isotherms exactly.
`--reuse-opacity WORK` finds matching isotherm inputs even when their indices
change after an extension, then verifies the complete run fingerprint.
`--initial-models WORK` can use completed structures from a matching family
as initial guesses; the final model is solved again with the selected source.
It also reads the compressed completed structures in the installed source
archive, verifying the original receipt against decompressed output bytes.
The command above uses those converged profiles to avoid repeating the
difficult initial convection search. Omitting this option generates coarse
starting models, which can need additional initialization work before every
cell converges.
It also accepts interrupted structures accompanied by `checkpoint.json`:
the reader checks their original input physics, executable/input fingerprint,
file hashes, grid labels and positive structure. Such a checkpoint makes no
final flux or convergence claim and cannot supply an imported grid cell.
The [archived checkpoint regression](https://github.com/oklo/ember/blob/e8f27e4f1eefdac7f4a25aa3c22a642475afb705/docs/results/nongrey_initial_checkpoint.json)
demonstrates this distinction and rejects altered or mislabeled inputs.
Opacity preparation uses at most four simultaneous source jobs; atmosphere
generation supports up to eight. Rejections are reported per model and never
filled into the grid.
`python3 scripts/nongrey_status.py WORK` reports saved isotherms, validated
atmospheres and the last-started source runs. A saved PID is diagnostic
information and is not treated as proof that a process is still running.

A complete accepted work directory can be archived without its large opacity
binaries and reimported offline:

```
python3 scripts/archive_nongrey_grid.py /tmp/ember-atmosphere-family \
  data/atmosphere/sources/nongrey_gs98_z020
python3 scripts/import_nongrey_grid.py \
  data/atmosphere/sources/nongrey_gs98_z020/manifest.json \
  data/atmosphere/nongrey_gs98_z020_tau100.dat
```

The archive checks source output receipts, each original opacity table's
composition and axes, and the table hash loaded by every atmosphere. It
retains the recipe, source provenance, isotherm inputs and diagnostics,
atmosphere inputs and outputs, and a hash inventory. Archiving cannot promote
an incomplete or rejected family.

The native dense synthesis executable is retained as `dense_synspec` in the
preparation receipt. The default sampling variant evaluates the original
line and continuum routines at the prescribed opacity wavelengths, using
at most 16 samples per continuum interpolation interval. It retains the
source line-selection threshold and broadening prescriptions. All selected
lines must fit in memory; it fails if the source starts chunking a list.
To fit the pinned inputs in the macOS executable address space, the extra
H2/He molecular broadening arrays are reserved for the third molecular list
(water); a nine-column list in another position is rejected.

The source patches also correct array bounds and output bookkeeping exposed
by the cool calculations: inconsistent abundance-array dimensions,
out-of-range helium-profile and CIA accesses, previous-line cache indices,
the H2- wavelength limit, and a final-row electron-density overwrite.
They preserve the dense algorithm for comparison. Bounds checking is enabled
in both builds. Source outputs retain full precision in the final depth
profile, so precision loss in printed temperatures cannot dominate the
boundary interpolation.
TLUSTY's atomic partition-function machinery supports elements through Zn.
All nonzero GS98 proxy elements lie in this range; heavier trace placeholders
remain in the chemistry bookkeeping without activating those atomic
queries. A source patch preserves their explicitly requested abundance
instead of substituting the default inactive-element abundance.
They also repair TLUSTY's He/H2 Rayleigh setup, which used an uninitialized
frequency, and bound both temperature indices after a Newton trial leaves
the opacity table. These changes are included in the pinned patch receipt.
In Rybicki mode the source disables the separate particle-conservation
unknown and updates density through `ELDENS`. The original molecular routine
then restored the old density and molecular weight whenever that unknown
was disabled. The patch retains the chemical density for Rybicki updates.
An independent final evaluation of density from chemical equilibrium at
each reported T/Pgas must agree with the structure to 0.2%. A pre-fix pilot
retained an atomic molecular weight and differed from FreeEOS by 54% at
the join; it is rejected and cannot supply physical grid cells.

The molecular-equilibrium stopping tolerance is tightened from `1e-3` to
`1e-8`. The original tolerance made entropy derivatives insufficiently
accurate for deep convection: a 300-depth model with tiny temperature
corrections still reported a 0.34% flux error. With tighter chemistry this
falls to 0.0012%. The importer verifies the loaded chemistry tolerance.
A source specification may declare `molecular_equilibrium_tolerance` below
its default of `1e-8`; the importer requires that value in the source log.
A nonideal atmosphere that solves density to `2e-9` uses `1e-11` for its inner
chemical equilibrium. The tighter inner solve avoids noise in the density
inversion. Flux, hydrostatic and source-coverage checks remain unchanged.
The final convection derivative uses `DERT=1e-5`; the thermodynamic entropy
differences retain the source's 1% temperature/pressure intervals. These are
separate numerical controls. SYNSPEC's opacity calculations retain their
native chemical stopping tolerance; they do not differentiate the entropy.

Zero-metal inputs encode absent elements as `1e-99`, because the external
programs interpret a zero abundance as solar. Both source patches exclude
molecules containing an element at this absence encoding (`<=1e-90`).
Otherwise, RUSSEL's minimum atomic pressure can produce substantial carbon
molecules in cold hydrogen–helium gas despite the declared absence of carbon.
The exclusion is independent of temperature and density; ordinary finite-metal
abundances are unchanged. [Native source checks](https://github.com/oklo/ember/blob/e8f27e4f1eefdac7f4a25aa3c22a642475afb705/docs/results/absent_element_chemistry.json)
recover the expected neutral H2/He particle mass and preserve warmer opacity
values to within 5.723e-5 in log opacity. This correction does not validate
finite-metal chemistry, condensates or arbitrary low-temperature coverage.

`prepare_nongrey_sources.py --cold-algebra` avoids overflow in cold LTE
opacity calculations. Boltzmann factors are combined before exponentiation;
Planck functions use negative exponents. Molecular-line strength is retained
for both line selection and wing extent. Matched spectra at 307.7, 1000 and
8000 K are unchanged in the tested H/He cases; 207 and 221 K cases are finite.
The 207 K control also passes with invalid-operation and overflow traps.
This numerical repair adds no opacity process or lower-temperature CIA data.
[Source checks](https://github.com/oklo/ember/blob/e8f27e4f1eefdac7f4a25aa3c22a642475afb705/docs/results/cold_opacity_algebra_sept30_v1.json).

The earlier `DERT=.001` produces a repeating small correction cycle in a
3000 K/logg5.15 solar-composition model. Reducing the derivative interval
to `1e-5` converges the retained structure in three iterations: maximum
temperature correction 4.4e-7, maximum local flux error 8.3e-6. Neither
the correction criterion nor the independent flux limit is relaxed.

An independent [convection Jacobian check](https://github.com/oklo/ember/blob/e8f27e4f1eefdac7f4a25aa3c22a642475afb705/docs/results/nongrey_convection_jacobian.json)
exposed two errors in RYBENE's interior integral-flux stencil: the diagonal
used the wrong pair of face derivatives, and the upper-face gradient
derivative omitted the radiative-loss factor used by the lower face.
The patch assembles the derivative of the actual flux difference and applies
the same radiative-loss factor at both faces. Both face temperature
derivatives now use the declared `DERT`. The convection flux law and energy
residual are unchanged. With the native source routines and an analytic
ideal molecular EOS, all nine checked matrix entries agree with centered
finite differences to 1.3e-6 relative; the original diagonal has the wrong
sign in every fixture. Run `scripts/check_nongrey_convection.py SOURCE WORK`
to repeat this source-level check with gfortran. The historical source patch
is retained with the validation files. At 2800 K/logg4.9/XH=.7/X3=0, the prior source needs Ng acceleration
and 28 iterations; the corrected source converges without acceleration in
eight iterations. Both final profiles pass the flux and chemistry checks;
their matching T/Pgas differ by only +0.000165%/-0.000786%. The installed
family has been reconverged with the corrected source.

The [native opacity derivative check](https://github.com/oklo/ember/blob/e8f27e4f1eefdac7f4a25aa3c22a642475afb705/docs/results/nongrey_opacity_derivative.json)
finds two further OPACTR errors. The forward opacity differences were divided
by the temperature at the last depth, rather than at the depth being
differentiated. The perturbation also left molecular weight and gas pressure
changed. The repair uses the local temperature, restores the saved particle
density/gas pressure, recomputes the base chemical state and restores the
original structural density, electron density and mean weight. An analytic
temperature-dependent EOS/opacity fixture exposes 50% and 25% derivative
errors and a 0.2% pressure drift before repair; all six corrected queries
agree with their prescribed secants and restore the state to 1.5e-14.
`scripts/check_nongrey_opacity_derivative.py SOURCE WORK` repeats the check.
Full solar and helium-rich control atmospheres pass the independent physical
checks; their matching states differ from earlier converged results by less
than one part per million. Earlier source patches remain archived so the
preceding numerical studies retain their original provenance.

The cool solver limits each temperature update to a factor of 1.03 and uses
`ILGDER=1`, so the initial model and final convection diagnostic use the
same logarithmic gradients as the Rybicki equations. With the inconsistent
default diagnostic, a pre-chemical-fix 2800 K pilot reported a 42% flux error despite small
iteration corrections; the consistent 100-depth calculation reports 0.24%.
The 0.2% final acceptance limit is unchanged. Optional `initial_depths` and
`initial_frequencies` settings solve a coarse starting atmosphere, then
interpolate only that initial guess onto a finer column mesh. The fine
calculation must independently converge at its requested resolution, and
its exact initial structure is archived and included in the run fingerprint.
New source runs also save a structure checkpoint after each iteration;
these initial guesses cannot substitute for the required final diagnostics.

The hottest high-gravity XH=.7/X3=.12 cell needs a convective initial
profile. Copying the completed X3=0 model directly can leave deep layers
initially stable and produce large radiative temperature corrections.
The archived `solar-he3-hot-convective-initialization` control enables the
source's `ICONRE=3, ICONRS=1, IDEEPC=3, CRFLIM=-1` refinement for the first
three iterations; `IMUCON=200` disables its later experimental branch.
It then converges the unchanged equilibrium equations in eight further
iterations, with flux error 6.3e-6 and maximum correction 3.3e-7. Its completed
target-composition structure supplies the initial guess for a separate
final cell using the canonical family parameters. This refinement changes
initialization, not the accepted flux law or convergence limits.

The archived [reference depth study](https://github.com/oklo/ember/blob/e8f27e4f1eefdac7f4a25aa3c22a642475afb705/docs/results/nongrey_source_benchmark.json)
uses the author's 5500 K, log g=4.5 input and original opacity table. Small
iteration corrections at 70 and 100 depths coexist with flux errors of
1.73% and 0.754%; both are rejected. At 200 depths the maximum flux error is
0.196%. This historical study predates the He/H2 Rayleigh correction and
uses the original gradient diagnostic. The capacity-only change reproduced
its T/P state exactly; that statement does not describe the later physical
Rayleigh fix.
This is a source regression, not validation of a cool helium atmosphere.

The [dense/sampled opacity comparison](https://github.com/oklo/ember/blob/e8f27e4f1eefdac7f4a25aa3c22a642475afb705/docs/results/nongrey_opacity_sampling.json)
covers 2800/6000 K, densities `1e-8`/`1e-5` g/cm³, XH=.55 and X3=.1.
Across those four states, Rosseland absorption means differ by at most
0.74%, and Planck absorption means by at most 0.60%. Individual spectral
samples can differ much more, particularly near narrow lines and continuum
edges. These averages do not validate the emergent spectrum or establish an
atmospheric error bar. Original binary tables, source inputs, logs and hashes
are archived under `data/atmosphere/sources/nongrey_validation/`; use
`scripts/compare_nongrey_opacity.py` to repeat the comparison.

The [material interpolation audit](https://github.com/oklo/ember/blob/e8f27e4f1eefdac7f4a25aa3c22a642475afb705/docs/results/nongrey_material_interpolation.json)
compares direct source opacities at four independent T/rho states with
TLUSTY's bilinear interpolation in logarithmic opacity, temperature and
density. At 2890/4290 K and densities `1.14e-5`/`1.47e-4` g/cm³, the
interpolated Rosseland absorption means exceed the direct values by up to
4.52%; Planck absorption means are lower by up to 10.81%. This is a larger
numerical uncertainty than the wavelength and depth checks. Halving the temperature spacing and reducing the density log10 spacing
from 0.556 to 0.394 over the held-out atmosphere's actual support changes
its matching T/Pgas by -0.00067%/+0.332%. The finer model independently
passes the flux, chemical-density and full-profile support checks. This
local boundary comparison does not bound errors elsewhere in the family. Archived reference stencils copy original source cells
exactly; the comparison can be repeated with
`scripts/compare_nongrey_opacity.py --interpolate-reference`.

## Computational cost

Atmospheres are computed once and reused across stellar timesteps. Chemistry,
radiative transfer and convergence dominate source cost. Opacity generation
is limited to four concurrent jobs because of the line-list memory footprint;
atmosphere generation supports up to eight. Use fewer workers when sharing
resources with stellar calculations. `scripts/profile_nongrey_source.py`
measures source routines in a separate instrumented build.

## Acceptance and interpolation

The importer requires all composition/Teff/gravity cells, archived source
inputs and outputs, and matching checksums. It checks the requested mixture,
baryonic element masses, loaded-mass marker, transfer mode, convection
parameter and depth/frequency settings. Each accepted structure must have:

- a complete final profile and independent final correction vector;
- maximum relative correction at most `1e-5`;
- total local flux within 0.2% of the requested flux at every depth;
- emergent bolometric flux consistent with the requested Teff;
- gas pressure consistent with the specified gravity to 1%;
- density consistent with an independent chemical-equilibrium evaluation
  at the final temperature and gas pressure to 0.2%;
- its entire converged T/rho profile inside the original opacity rectangle;
- matching optical depth strictly bracketed by source depth points.

Intermediate source iterations can probe outside the opacity rectangle.
TLUSTY's internal clipping/extrapolation during those iterations is not a
proof of source coverage; the independent final-profile check is mandatory.
The acceptance tolerances measure solver consistency, not physical accuracy.
Opacity sampling, material-grid interpolation, depth resolution and
composition interpolation need separate convergence checks.

The imported state is log T and log gas pressure interpolated in log
Rosseland depth. Runtime interpolation uses convex multilinear weights in
`(XH, X3, log10 Teff, log10 g)`. Analytic Teff, gravity and composition
derivatives come from that same interpolant, including cross terms.
Radiation pressure `a T^4/3` is added once; density is obtained from the
caller's EOS at the actual composition. Derivatives are one-sided at grid
knots, as appropriate for a multilinear interpolant.

There is no composition clipping, extrapolation, missing-cell interpolation
or fallback atmosphere. A grid's documented mixture approximations require
explicit selection. The evolution application accepts:

```
build/apps/ember-evolve 512 1e8 1e7 1 sfii-svh wd \
  nongrey:/absolute/path/to/atmosphere.dat
```

The non-grey evolution mode chooses its trial radius and effective
temperature from the interior of the supplied Teff/gravity rectangle.
It scales only the starting luminosity guess; relaxation still solves the
original energy equations. The previous 0.15 Rsun polytropic guess has
Teff≈2564 K and can lie below a cool-dwarf grid starting at 2600 K.

Synthetic C++ fixtures test the reader and its derivatives. They are not
physical atmosphere data. Numerical agreement at sampled matching states
does not establish a physical uncertainty bound over an entire cooling track.
