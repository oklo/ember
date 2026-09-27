# One evolution program

`ember-evolve --lifetime` connects the common stellar solver to a contracting
Hayashi starting model, explicit initial D and pp/CN burning, plasma neutrino
losses, face luminosities and the finite metal-transport interface. It supports
the wholly convective limit and an explicitly selected hot diffusion law with a
mixed cool envelope. Continuous late-life atmosphere selection remains under
integration; the complete Hayashi-to-white-dwarf calculation is not yet established.

```
ember-evolve --lifetime CONFIG TARGET_YEARS NEW_OUTPUT_DIRECTORY \
  [--restart CHECKPOINT] [--max-steps N] [--cpu-seconds S]
```

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
The configured absolute species time-error allowance may be at most **1e-5**.
A matched **500 Myr** comparison after the centre becomes radiative uses
**17** accepted intervals at 1e-5 versus **45** at 1e-6, with maximum relative
temperature difference **2.475e-6** and total helium-3 difference **2.766e-5**.
The radiative-core continuation selects 1e-5; the structure allowance remains
1e-4, nuclear-power allowance .005 and global inventory budget 1e-14.
This comparison does not establish accuracy during a flash.
[Matched-age control](results/radiative_core_time_accuracy_sept27_v2.json).
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
`attempts.jsonl` retains acceptance and conservation checks. Only the initial,
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

The retained 1-Gyr contraction reproduces all 1020 physical history rows and
both complete physical checkpoints exactly with binary loading. In one matched
M4 check, CPU time through construction of the initial relaxed model fell from
6.748 to 2.591 s, including input checks and loading. This is a startup improvement;
it does not imply the same speedup for a long evolutionary calculation.
