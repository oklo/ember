# One evolution program

`ember-evolve --lifetime` connects the common stellar solver to a contracting
Hayashi starting model, explicit initial D and pp/CN burning, plasma neutrino
losses, face luminosities and the finite metal-transport interface. It supports
the wholly convective limit and an explicitly selected hot diffusion law with a
mixed cool envelope. The atmosphere selection keeps the same matching depth
through the supplied composition and gravity coverage. The complete Hayashi-to-white-dwarf
calculation is not yet established.

```
ember-evolve --lifetime CONFIG TARGET_YEARS NEW_OUTPUT_DIRECTORY \
  [--restart CHECKPOINT] [--restart-step-years Y] [--max-steps N] [--cpu-seconds S]
```

A restart normally retains the saved next trial step. Optional
`--restart-step-years Y` changes only that trial interval, capped by
`maximum_step_years`; the saved structure, age and composition are unchanged.
This is useful after an input-coverage stop has driven the saved interval too
small for a well-conditioned thermal solve. Every trial still passes the normal
time-error and conservation checks. The requested and effective intervals are
recorded in `execution.json`.

The optional [integrated outer envelope](ENVELOPE.md) uses the same boundary
solver from Hayashi contraction onward and includes its mass in conservation.

The configuration is a text file containing one `key "value"` per line.
Version 1 requires the following keys; paths are relative to the configuration:

- `eos`, `opacity_low`, `opacity_warm`, `opacity_bridge`, `opacity_hot`,
  `conduction`, `atmosphere`, `collisions`, `composition`.
- `mass_Msun`, `initial_radius_Rsun`, `initial_Teff_K`,
  `initial_entropy_loss`, `points`, `zone_threads`.
- `screened_minimum_T_K`, `initial_step_years`, `maximum_step_years`,
  `structure_tolerance`, `species_tolerance`, `energy_tolerance`,
  `coupling_abundance_tolerance`, `inventory_abundance_tolerance`, and
  `version` (set to `1`).

For screened-core transport, optional `screened_heat_lower_T_K` separates the
start of the heat-law transition from the lower temperature limit of the
microscopic species law. It defaults to `screened_minimum_T_K` and must lie
between that limit and `screened_heat_upper_T_K`. This permits independently
checked species coverage to grow while retaining the same heat prescription.
It does not establish ionization or collision-table support; those still need
physical assessment. Optional `quantum_screening_zeta_max` selects the
[nuclear quantum correction](NUCLEAR.md); its default is zero.

The optional `structure_prediction "linear"` uses the recent evolution to
estimate the next radius, density, temperature and luminosity before solving.
The accepted composition and thermal history remain the physical starting
state. Every full step and half step retains the same convergence, conservation
and time-accuracy checks. If the predicted solve fails, the program retries
from the previous structure. The default is `"none"`; after a restart, the
prediction history is rebuilt from newly accepted states.

The physical-metal burning/diffusion solve also reuses the latest composition
from the coupled iteration as its starting guess. Together with structure
prediction, this reduced a matched 1 Gyr interval from **297.8 to 143.5 CPU
seconds**, with negligible structural and fuel differences. This is one timing
comparison in a gradual phase, not a benchmark for the complete track.
[Comparison](results/solver_starting_guesses_sept27_v1.json).

Optional `collision_taylor_radius "1e-4"` reuses collision responses with a
first-order Taylor expansion around the last exact evaluation at each face.
The radius bounds changes in ln T, ln density, H, helium-3, total metals and
ln screening length; population changes also stay within 1% of their previous
values. A changed active species set, source table, or nonpositive response
requires an exact evaluation. Source-table bounds are checked on every query.
The approximation changes transport coefficients, not their conservative
exchange between neighboring cells. Its radius is recorded in the checkpoint;
the default is zero (exact evaluation), and the supported maximum is 1e-4.

With both starting-guess improvements, this reduced the same 1 Gyr comparison
from **143.5 to 81.76 CPU seconds**: **3.642 times** faster than the original
**297.8 seconds**, retaining full-step/two-half-step checks. Relative global
changes were at most **1.022e-9**. A separate four-interval check recomputed
every reuse exactly, measuring a maximum coefficient-group error of
**1.195e-8**. These are gradual-phase measurements, not a complete-track or
flash-convergence benchmark. Optional `collision_verify_reuse "1"` enables
that diagnostic; `collision_reuse.json` records counts and measured errors.
[Collision comparison](results/collision_reuse_sept27_v1.json).

The recovered high-hydrogen atmosphere families overlap over **X = 0.99–0.995**.
The lower family remains valid through X = 0.995; the next reference requires
Z no greater than **0.0001**. This overlap supports the current track as metals
settle, while preserving the measured columns and optical depth 100. A local
comparison found matching-temperature differences below **0.5354%** and
gas-pressure differences below **0.6945%** between the overlapping descriptions.
The accepted stellar continuation preserves its starting physical state and
passes the existing timestep and conservation checks.
[Atmosphere overlap](results/hydrogen_atmosphere_join_sept27_v1.json).

The initial composition file contains the nine baryonic mass fractions in
`Composition` order. The CN inventory starts with the declared GS98 isotope
mixture and evolves thereafter. Mass and starting-state parameters are inputs;
the main sequence is an outcome. The entropy-loss source constructs only the
initial state. Subsequent evolution uses the internal-energy change and work
term, without an imposed luminosity or composition reset.

During primordial-D burning, the selected early approximation requires an
exactly homogeneous, wholly convective star. Every solved full/half interval
checks its estimated mixing heat and convective travel time. The sum of the
absolute species heat contributions must remain below **0.01000%** of local
luminosity, and the estimated D nonuniformity below **1e-8** in mass fraction.
The controls supporting these choices are in
[the initial transport assessment](results/pms_initial_transport_bound_sept26_v1.json).
These limits bound this particular mixing approximation; they do not validate
a microscopic kinetic-heat or partially ionized collision prescription.

After initial D reaches zero through the burn solver, total H/He3/metal
composition enthalpy enters the heat equation explicitly. No abundance is
manually reset. The default `transport "whole_convective"` requires the star
to remain wholly convective and homogeneous. A mixing-gradient
estimate must remain below **1e-8** in absolute mass fraction; a generous
microscopic-drift estimate below **1e-6**, and its kinetic-heat estimate below
**5%** of local luminosity. The latter two are order-of-magnitude proxies,
not validated partially ionized collision coefficients. Bounded stellar
controls vary conduction and residual heat to assess their structural effect.
[Transport assessment](results/convective_transport_integration_sept26_v1.json).

To permit hot radiative regions, select:

```text
transport "screened_core"
screened_minimum_T_K "2000000"
screened_heat_upper_T_K "3000000"
maximum_relative_mixing_gradient ".01"
```

Microscopic H/He/metal fluxes then act between mixed regions and through the
radiative core. All metals share a mass velocity with separate collision charge
contributions. Both endpoints of a species boundary must meet the microscopic
law's domain. Its fully stripped metal approximation remains conditional;
the temperature threshold alone does not demonstrate complete ionization.
There is no zero-flux replacement in an unsupported cool radiative layer.

The heat response joins material conduction and total composition enthalpy in
the cool envelope to the full screened heat law over the stated temperature
interval. Analytic derivatives include the changing join weight. This does not
attenuate microscopic species exchange. Each convective region is checked for
homogeneity and the gradient needed to carry its composition flux; the configured
allowance is a fraction of its H, He3 or total-metal abundance. The same cool
drift and kinetic-heat estimates described above remain in force. A model that
needs greater abundance gradients requires finite mixing. These checks assess
the instantaneous-convection approximation, not every uncertainty in diffusion.

The join agrees with the existing stellar transport implementation at saved
faces. A short stellar comparison and an independent radiative-core/mixed-envelope
control pass, including explicit rejection of an unsupported cool boundary.
[Integration evidence](results/envelope_transport_integration_sept26_v1.json).
For phases where burning and mixing times become comparable, the same program
can select `convective_mixing "finite_implicit"` or `"finite_lagged"` with
`transport "screened_core"`. The default is `"instantaneous"`. Finite implicit
mixing recomputes the mass conductance during the coupled solve; finite lagged
mixing holds its value from the starting model and requires saved composition
heat rates for a stratified restart.

The optional `instantaneous_mixing_below_T_K` retains instantaneous mixing on
convective faces with either endpoint below that temperature. Zero, the default,
leaves all faces finite. The hot microscopic law cannot support such a choice
in a cool envelope; an explicit cool mixing approximation is then needed. Its
gradient and drift/heat checks remain in force. Radiative faces are never merged,
and every exposed interface must satisfy the microscopic provider's domain.
The numerical mixing temperature is independent of that physical domain limit.
Both the physical convection regions and numerical composition blocks are
checked on the returned structure, after its composition and heat update.

Initial D still uses the assessed wholly convective approximation. Once D is
exhausted, each interval selects the configured mixing treatment from its own
starting model, including the second half of a time-accuracy comparison.
The mixing mode and temperature enter restart identity. D exhaustion itself
is not a criterion for losing convection or requiring finite mixing.

A **12.73 Myr** comparison at **3.609 Tyr** passes every full/half-interval
audit. Finite mixing changes the helium-3 inventory by **5.107e-15** relative
to instantaneous mixing in this slow phase. The default program reproduces
three intervals and the final physical checkpoint exactly. This establishes
integration for the tested states, not convergence through a flash.
[Driver and transport checks](results/finite_mixing_driver_sept27_v1.json).

Atmosphere overlaps, the bounded approximation for small metal changes and the
measured response to settling are described in [ATMOSPHERE.md](ATMOSPHERE.md).
`atmosphere_metal_chain` selects the measured intervals in the same driver;
the contraction boundary remains exactly selected at the initial Z. These
selections, all response files, the transport choice, join temperatures and
mixing allowance enter restart identity.

`atmosphere_hydrogen_interval` adds the retained hydrogen-rich reference
families. `atmosphere_hydrogen_envelope` adds explicit trace-helium and nearly
pure-hydrogen gas sources, joined over declared composition and gravity
intervals at the same optical depth of 100. The initial boundary remains
unchanged. All source bytes and interval coordinates enter restart identity,
and the runtime packager includes their complete input set.
The cool, high-gravity source contains bounded inferred rows; using that family
does not establish independently solved coverage throughout the late track.
[Source limits](../data/atmosphere/lifetime_hydrogen_envelope/README.md),
[coverage and exact continuation checks](results/hydrogen_envelope_integration_sept27_v1.json).

The volume-face thermal gradient now uses the same logarithmic-temperature
mean with and without a material-heat provider. A zero heat/conduction
provider leaves structure equations and convective mixing unchanged. Older
nodal calculations retain their original discretization.

One full interval is compared with two half intervals. Accepted models use
the latter. Structure, physical isotope abundances, surface luminosity and
nuclear energy control the timestep. Separate global isotope and mass-defect
budgets and the first law check every interval. The nuclear-error scale is the
larger of nuclear and surface power, so vanishing nuclear power does not impose
an inappropriate relative-precision requirement.
If any full or half interval fails a conservation check, the driver rejects
that trial and halves the timestep. It preserves the accepted model and heat
rates and stops after eight consecutive rejections. This retry policy changes
neither the budgets nor which intervals must pass them. The independent
controller can still request an immediate stop on audit failure.
The nonlinear composition tolerance and global inventory budget are independent.
The former limits a local correction; the latter checks mass-weighted global
conservation. A stratified radiative-core control converges with a **1e-12**
local tolerance while passing the unchanged **1e-14** inventory budget and
energy checks; forcing the local correction to **1e-15** stalls its line search.
The global budget need not be larger than the local correction tolerance.
Neither setting replaces the full-step/two-half-step accuracy check.
The configured absolute species time-error allowance may be at most **1e-4**.
During gradual radiative-interior evolution, a **1 Gyr** comparison near
**3.626 Tyr** supports that setting: it uses **20** accepted intervals instead
of **39** at 1e-5 and **18.03%** less CPU. Surface temperature changes by less
than **0.001 K** and the helium-3 inventory by **0.000184%**. The largest local
helium-3 difference is **4.027e-7** in mass fraction at the envelope boundary.
No cell is omitted from the time-error estimate. Structure accuracy remains
1e-4, nuclear-power accuracy .005 and the global inventory budget 1e-14.
This comparison supports the gradual phase, not time resolution during a flash.
[Matched-age comparison](results/envelope_time_accuracy_sept27_v1.json).
The preceding 1e-5 setting has an independent
[500 Myr comparison](results/radiative_core_time_accuracy_sept27_v2.json).
The continuous star now uses a tighter **1e-13** stratified solve after a matched
**11.76-Myr** comparison at radiative-core onset. This reduces audit retries
without changing the global budget; luminosities differ by **0.00004235%**.
[Comparison and exact continuation](results/radiative_core_solve_accuracy_sept27_v1.json).
When instantaneous convection mixes the whole star, the driver tightens the
composition correction to the smaller of the selected tolerance and **1e-15**.
This inexpensive one-region solve needs that precision during initial
deuterium burning: a **1000-year** test with **1e-12** throughout fails the
unchanged inventory and mass-power checks. The tighter setting follows the
actual mixing regions at each coupling iteration, including a newly formed
radiative core; it does not depend on an assigned evolutionary phase or age.
The general stratified tolerance applies once more than one region is present.
Both settings are bound into restart identity, and all interval audits remain
unchanged.

The optional `opacity_hydrogen_response` family selects the existing bounded
composition extensions through `RadiativeOpacity`. Its bounds are explicit:
`opacity_minimum_Z`, `opacity_maximum_Z`, and `opacity_maximum_X`. The prepared
selection uses **0**, **0.16**, and **1**, respectively. Low-temperature tables
interpolate at fixed hydrogen share; the final **1e-7** in hydrogen at zero
metals uses the measured source slope. Warm/bridge hydrogen and low/high-metal
extensions retain their earlier linear approximations. Hot hydrogen opacity
uses OPLIB composition ratios, joined to the TOPS composition slope at high
density over **log R = 1.25–1.45**, bounded by **log R = 1.8** and
**log T = 5.6–6.1**. These are explicit composition approximations, not new
opacity calculations; temperature and density coverage is still enforced.
The selected linear high-metal method does not load the unused second metal
table family. The opacity contributes radiation only; the heat transport
implementation supplies conduction. Selection and all source bytes are bound
into restart identity. See `docs/results/opacity_lifetime_extension_sept27_v1.json`
for saved-profile checks and the earlier physical sensitivity tests.

`history.jsonl` retains every accepted model's scalar diagnostics;
`attempts.jsonl` retains acceptance and conservation checks. The history includes
central density, hydrogen and helium-3, and the locations of maximum helium-3
abundance and specific nuclear power. Here `q` means enclosed mass divided by
stellar mass. `burning_half_max_*` describes the contiguous cells surrounding
the strongest specific nuclear power, above half its maximum: cell count and
inner/outer mass faces. Separate burning peaks are not combined. These fields
help identify a narrowing shell; they do not establish mass convergence. If
nuclear power is zero throughout, the peak location, width and cell count are
zero. The helium-3 peak location is zero when helium-3 is absent. The diagnostics
reuse the nuclear rates already evaluated for the luminosity integral.
Only the initial,
latest and final full checkpoints are retained by default. Checkpoints include
CN/D abundances, the luminosity convention and, when applicable, material-heat
rates. Restart checks bind the executable, physical settings, family coordinates
and every selected table's contents. Configuration and table paths, comments,
thread counts, run duration, timestep ceilings and execution budgets may change.
Changed accuracy settings or table contents are rejected before material tables
are loaded. `execution.json` records the controls for each invocation. An exact
replay requires retaining the same effective timestep limits; changing them is
a new numerical control, not a promise of identical saved steps.
Older checkpoints remain readable with their frozen executable; the new driver
does not silently reinterpret a checkpoint from a different executable or
physics selection. `report.json` distinguishes the requested age from a planned stop or a
physical/source-domain failure. A bounded run finishing is not a claim that the
full evolutionary track has been completed.

To avoid repeatedly parsing the EOS text planes, pack the selected family once:

```sh
build/apps/ember-pack-eos data/production/eos.dat data/production/eos.bin
```

Set `eos` in a copied configuration to the binary file, using a path relative to
that configuration. The bundle contains the same source values, masks and
logarithmic coordinates as the text input, in little-endian IEEE binary64.
Interpolation and the construction of composition derivatives are unchanged.
It is a self-contained input, checked by content; the loader does not trust
modification times or silently substitute a cache. Keep the original source
manifest and provenance for regeneration. An existing output is never overwritten.

The optional `eos_low_metal_interpolation "quadratic"` uses the first three
metal planes below the first nonzero metal node. Between the first and second
nonzero nodes, a quintic weight returns to the cubic interpolant with continuous
first and second derivatives. For the current family these boundaries are
Z = 0.005 and 0.02. The default remains `"cubic"`. This setting is part of the
restart identity.

All thermal and composition derivatives come from the same free-energy
interpolation; the analytic ideal mixing terms and actual stellar abundances
are retained. This avoids consulting a distant metal-rich source at trace Z.
It does not extend the temperature/density domain of the required low-metal
planes. A matched 10 Myr cooling comparison changes luminosity by 3.637e-6
and He3 inventory by 2.001e-7; the convective boundary is identical. Metal
enthalpy derivatives can change appreciably even when their contribution at
trace abundance is small. The comparison and derivative checks are recorded in
[the EOS assessment](results/low_metal_eos_sept27_v1.json).

The retained 1-Gyr contraction reproduces all 1020 physical history rows and
both complete physical checkpoints exactly with binary loading. In one matched
M4 check, CPU time through construction of the initial relaxed model fell from
6.748 to 2.591 s, including input checks and loading. This is a startup improvement;
it does not imply the same speedup for a long evolutionary calculation.


### Optional reuse of nearby transport states

`eos_taylor_radius "1e-4"`, `buoyancy_reuse_spacing "1e-4"`, and
`screening_reuse_spacing "1e-4"` select first-order response reuse. Each defaults
to zero (exact evaluation), and the driver limits each to 1e-4. The selected
values enter the restart identity. `verify_response_reuse "1"` recomputes every
EOS and buoyancy reuse exactly and writes its measured errors to
`response_reuse.json`; it is intended for short comparisons.

Transport EOS reuse retains the source masks for every requested derivative
plane. The composition radius is relative to each species, including reference
helium-4. Active-species changes force exact evaluation. Buoyancy and nuclear
screening use fixed logarithmic grid anchors. Buoyancy uses exact inversion
near a source boundary. Configure these process-wide options before starting
workers; configuring buoyancy reuse clears the previous calculation's anchors.
Per-point transport EOS anchors may depend on thread scheduling within the
measured approximation error.

Convection-region tests now share the configured zone threads; the connected
regions are assembled in mesh order. The driver permits up to 16 threads;
eight are selected on the M4 Max after the additional parallel work below. This does not change the mixing prescription.

A matched 1-Gyr segment starting at 3.718 Tyr, with four threads and every
full/two-half timestep check retained, used **48.97 CPU seconds** and
**20.34 wall seconds**, versus **78.08** and **41.37** for the preceding code.
Both accepted 28 intervals. Relative changes in global stellar quantities were
below **2e-12**; all isotope and energy audits passed. The concurrent timing
pair is an estimate of the gain, not a hardware-independent benchmark.
Over four further checked intervals, more than 460,000 exact comparisons bounded EOS
potential/enthalpy group errors by **2.782e-8**; the maximum absolute error in
the buoyancy term was **1.634e-6**. These checks concern gradual shell burning,
not convergence through a flash. Sampled timestep checking remains unselected.
[Comparison and verification](results/response_reuse_sept27_v1.json).


### Coupled-solve accuracy and additional parallel work

The optional settings `coupling_stop_tolerance "1e-10"`,
`material_heat_tolerance "1e-7"`, `verification_residual_tolerance "1e-8"`
and `verification_correction_tolerance "1e-7"` separate the outer iteration
from the inner species solve and its integrated conservation bounds. The
returned composition still passes a structure residual and correction check.
Omitting these settings retains the previous stopping criteria. Selected values
enter the restart identity; the driver requires them to remain smaller than
the corresponding time-discretization tolerances.

Face mixing coefficients, transported-heat residuals and nuclear sources used
in flux reconstruction now share the zone threads. Ordered conservation sums
remain serial. On the same Hayashi-origin star at **3.841 Tyr**, a **5 Gyr**
comparison took **31.53 wall seconds** with eight threads, versus **54.95 seconds**
with four threads and the preceding code: **1.743 times** faster. CPU use was
**154.7** versus **142.0 seconds**; the elapsed-time gain here comes from using
more cores. Both retained 75 intervals and every full/two-half check.
Relative global differences were at most **1.068e-12** and the helium-3 profile
difference, normalized by its peak, was **5.861e-11**.

Five integration tests cover ordinary, finite and lagged convection, unchanged
conservation bounds, parallel face coefficients, actual convective-region sets
and invalid solver settings. The production histories record boundary counts;
they do not independently establish identical boundary positions. This
comparison covers gradual shell burning; sampled timestep control and flash
convergence remain separate work.
[Measurements](results/parallel_coupling_sept27_v1.json).


### Refining a species solve before rejection

A small Newton correction and a closed global inventory can still leave a
local species balance outside its existing tolerance. In that case, the
solver refines only the species solution, starting from its returned answer,
then reconstructs and checks the fluxes again. Failure still rejects the step;
no abundance, flux or conservation requirement is adjusted afterward.

Over a matched **1 Gyr** interval near **3.592 Tyr**, this removes three
local-continuity rejections and reduces CPU use by **5.699%**. The altered
step sequence changes luminosity by **2.712e-6** relative and total He3 by
**6.851e-6**. Six transport and mixing tests pass.
[Comparison and test scope](results/species_continuity_refinement_sept28_v1.json).

### Composition response in the structure iteration

The optional `linearized_burning "1"` includes the change in nuclear heating
caused by the local end-of-step composition responding to temperature and
density. Three auxiliary implicit burns estimate that response once per solve.
The structure iteration includes the resulting linear terms, which vanish as
the coupled solution converges. Direct temperature and density derivatives
already supplied by the nuclear network are not counted twice.

This is an approximate Jacobian. It omits transport feedback and applies no
response to instantaneously mixed regions. With finite mixing, its local
response likewise omits mixing during the auxiliary burns. Every returned
assisted solution passes the original structure equations, physical convection
test, species conservation and energy checks. An unsuccessful auxiliary burn
falls back to the existing iteration. The option is off by default and enters
the restart identity when enabled. It works with `structure_prediction "linear"`.

For the same Hayashi-origin star near **3.997 Tyr**, a **1 Gyr** comparison
with two threads per run used **15.47 CPU seconds** versus **149.5 seconds**
with the current production solver. Measured elapsed times were **10.14** and
**79.28 seconds**: **7.815 times** faster in this phase, with **9.663 times**
less CPU. Both runs retained every full-step/two-half-step check and the same
acceptance thresholds. The improved iteration accepted eight longer intervals,
versus 49; rejections fell from 30 to one.

The maximum relative global difference was **0.01382%**. The largest local
ln-density difference was **2.994e-4**, and the helium-3 profile difference
normalized by its peak was **1.102e-4**. Convective mass and boundary counts
agreed. Five tests cover the differentiated structure equations and ordinary,
finite and lagged mixing. This single comparison measures gradual shell
burning under shared machine load; it does not establish flash convergence
or a whole-lifetime speed-up.
[Comparison](results/linearized_burning_sept27_v1.json).

## Dense hydrogen opacity in a convective envelope

`opacity_dense_hydrogen_maximum_logR` selects the upper density coordinate for
an existing hydrogen-composition approximation; its default is 1.8 and its
permitted interval is 1.8--2.5. Here R is rho/(T/1e6)^3 in CGS units. A
nondefault selection is stored in restart identity. The original temperature
bounds, source-table checks and smooth join are retained. This option does
not extrapolate temperature or density beyond the underlying source tables.

The calculation uses the opacity slope between the retained X=0.70 and
X=0.75 planes to continue to the actual hydrogen abundance. Near the cooling
model's density limit, at 415400 K, convection carries 99.95% of the thermal
luminosity. The selected upper bound is 2.5. A matched 1 Myr comparison with
half and twice the approximated opacity changes global quantities by at most
4.260e-12 fractionally; all time and conservation checks pass. Completed
atmosphere coverage limits the duration of this comparison. The source data
support the larger density interval, but the composition continuation remains
an approximation and needs reassessment if these layers become radiative.
[Measured response](results/dense_envelope_density_extension_sept27_v1.json).


## Second-order accepted states

`richardson_extrapolation "1"` enables an optional second-order estimate from
`2 * (two half-steps) - (full step)`. The default is off. Every interval retains
the full/two-half time-error check and the three original conservation audits;
there is no increase in the allowed time error. Structure variables, abundances
and stored species heat rates use the same combination.

The candidate must retain nonnegative abundances and the same convective and
instantaneously mixed regions as the starting model, full step and both half
steps. Each assessment uses that model's stored species heat rates. Atmosphere,
EOS and opacity support, hydrostatic balance, mass continuity, temperature
transport and boundary conditions are checked directly. The dimensionless
algebraic allowance is 1e-5; the global energy-measure allowance is 0.001 of the
larger radiated or nuclear energy over the interval. Species accounting must
meet the configured species time accuracy. This global energy measure differs
slightly from the local discrete first law and does not replace its audits.
The backward-Euler thermal equation is not imposed on a second-order estimate.
Any failed assessment uses the ordinary two-half-step result. Initial-deuterium
evolution also uses that fallback. Integrator selection and assessment limits
are recorded in checkpoint identity.

A matched 1 Gyr comparison near 3.997 Tyr, using the current physics and two
threads, compared both methods with a 5 Myr-step reference. The reference ends
0.03955 yr short because of accumulated absolute-age rounding; this negligible
offset is included in the receipt. Richardson reduces the helium-3 peak error
from 0.005233% to 0.0004588%, while the largest local structure difference is
0.04542%. It costs 20.61 retained CPU seconds versus 14.64 for the ordinary
method. It remains off in the production trajectory: this comparison supports
its accuracy but does not show a speed gain or validate rapid burning. The
controller tests cover missing assessment, declined and throwing assessments,
and acceptance; a real initial-deuterium step confirms the application fallback.
[Comparison](results/richardson_current_physics_sept27_v1.json).

## Abundance changes during instantaneous mixing

`abundance_cap_after_mixing "1"` measures the abundance cap from the initial
composition conservatively averaged over each instantaneous mixed region.
It prevents newly incorporated material from imposing a finite composition
jump that cannot be reduced by shortening the timestep. Radiative cells keep
the pointwise comparison; finite and lagged finite mixing retain the original
cap regardless of this setting. The default is off. The selected cap reference
is recorded in restart identity.

Every full and half interval still passes the usual composition, structure
and energy error checks. Linear structure prediction is skipped after a
preceding abundance jump larger than the cap, because instantaneous mixing
is not proportional to elapsed time. This only changes the initial guess.
A resolved envelope-boundary comparison and finite-mixing regression checks
are recorded in [the numerical comparison](results/instantaneous_mixing_step_cap_sept27_v1.json).

The optional `opacity_dense_hydrogen_maximum_logT` extends the declared
source-slope hydrogen-opacity approximation from its default upper log T of
6.1 to at most 6.3. It preserves the source temperature/density checks and
the selected log R limit. It is recorded in restart identity. A 20 Myr
half/double-opacity comparison supports this choice in the efficiently
convective envelope of the current remnant; it does not establish accurate
radiative opacity for arbitrary pure-hydrogen layers. See the
[temperature-domain comparison](results/dense_envelope_temperature_extension_sept27_v1.json).


### Dense gas opacity in a cool envelope

`opacity_cold_dense` selects the compact
[`cold_dense_hydrogen.dat`](../data/opacity/cold_dense_hydrogen.dat) gas table
for dense, nearly metal-free hydrogen envelopes. It replaces the low-temperature
source smoothly over log R = 5.6–5.9, where R = rho/(T/10^6)^3. Temperature
joins span 3000–3500 K and 10000–12000 K; composition joins limit its use to
hydrogen-rich mixtures with Z below 1e-10. Source temperature, density and
composition bounds remain enforced. The integer-mass source convention is
converted explicitly, with derivatives that preserve extinction per length.

The table includes gas absorption and electron scattering; nonideal chemistry,
grain opacity and direct H3+ lines remain incomplete. Its use is supported by
the envelope's weak sensitivity to opacity: multipliers 0.1 and 10 change
luminosity by at most 0.0001184% over a 140 Myr cooling comparison, ending at
1902 K. `opacity_cold_dense_scale` selects those controls; the default is 1.
Both the file hash and scale enter restart identity. These checks support
this cooling segment, not a final cooling age.
[Source details](../data/opacity/cold_dense_hydrogen.json) ·
[Numerical checks](results/cold_dense_gas_opacity_sept30_v1.json).

### Radiative opacity in the conductive interior

`opacity_conductive_interior "1"` permits a bounded density continuation of
radiative opacity in cool, dense cells. It measures the source opacity slope
near 8500 g/cm³ and joins smoothly over 8500–9500 g/cm³. The original source
is recovered over 3.4–3.6 MK. Temperature and composition still require
supported source queries at the anchor; the option does not supply missing
opacity data. The assessed lower temperature limit is 800 kK.

The continuation is allowed only when a factor-ten uncertainty in radiation
changes total heat conductivity by at most 0.1%. The check uses half the
tabulated electron conductivity, conservatively below the microscopic value
measured along the tested cooling segment. This margin must be reassessed
as the physical conditions change. `opacity_conductive_scale` accepts values
from 0.1 to 10 for sensitivity calculations; its default is 1. Both selections
are recorded in restart identity. `conductive_opacity.json` records the number
of continued evaluations and largest conditional transport uncertainty.

The 150 Myr cooling comparison changes luminosity by at most 0.008869% for
the two extreme opacity scales. This supports the approximation in that
segment, conditional on its declared uncertainty; it is not validation of
the opacity source outside its tabulated domain.
[Numerical and stellar checks](results/conductive_interior_opacity_sept28_v1.json).

The separate hydrogen-envelope option, `opacity_conductive_envelope "1"`,
continues the measured density slope where radiation is a small contribution
alongside electron conduction. Its density blend ends at 190 g/cm³, before
one hot-source corner at 199.5 g/cm³. Temperature, composition and source
checks remain active. Restart identity records
`hydrogen_density_continuation.v4`. A smooth blend recovers the original
radiative source between 320 and 300 kK, where the conducting channel turns
off. The source's coverage limits still apply. Missing conduction therefore
cannot justify radiative extrapolation.

The join and its derivatives pass native checks; the largest measured change
in total diffusive conductivity within the supported cold overlap is
1.808e-8. A matched 10 Myr stellar comparison leaves the reported global
quantities unchanged. Cooling then continues to the 2000 K atmosphere edge.
These checks support the opacity join, not the complete low-temperature
physics or a final cooling age.
[Native and stellar checks](results/conductive_opacity_cold_join_sept30_v1.json).

`opacity_conductive_envelope_uncertainty` selects the maximum fractional change
in total diffusive conductivity for a factor-100 reduction of the continued
radiative opacity. The default is 0.001; values up to 0.01 permit explicit
sensitivity studies. Any nondefault value enters restart identity. This is
an assumed opacity range, not a measured physical uncertainty. Select a
larger allowance only after checking its effect on the stellar calculation.
The envelope opacity multiplier accepts 0.01–100 for those comparisons.

In a matched 1 Gyr comparison, the two extreme opacity multipliers change
luminosity by at most 0.04833% and effective temperature by 0.01166%.
The largest exercised local transport bound is 0.001142. An allowance of
0.003 supports this segment; it is not a general cooling-age error bound.
The corrected density blend changes luminosity by 0.005844% in its separate
1 Gyr comparison. All energy and time-step checks are unchanged.
[Numerical checks](results/conductive_envelope_opacity_sept30_v1.json).

In the assessed H-rich cooling-envelope layers, direct pressure-ionization
checks support `screened_minimum_T_K "400000"`. This changes a domain check,
not the collision law. A matched 10 Myr continuation leaves global quantities
unchanged; cooler evolution still requires covered EOS and opacity inputs.
[Domain assessment](results/cold_transport_domain_sept30_v1.json).
