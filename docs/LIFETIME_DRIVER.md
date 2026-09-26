# One evolution program

`ember-evolve --lifetime` connects the common stellar solver to a contracting
Hayashi starting model, explicit initial D and pp/CN burning, plasma neutrino
losses, face luminosities and the finite metal-transport interface. It is under
integration. The wholly convective limit is assessed through main-sequence arrival.
Radiative microscopic transport and complete atmosphere coverage must still
be supplied before this is the continuous production calculation.

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

After initial D reaches zero through the burn solver, the same instantaneous
convective mixing remains selected. Total H/He3/metal composition enthalpy
then enters the heat equation explicitly. No abundance is manually reset.
The star must remain wholly convective and homogeneous. A mixing-gradient
estimate must remain below **1e-8** in absolute mass fraction; a generous
microscopic-drift estimate below **1e-6**, and its kinetic-heat estimate below
**5%** of local luminosity. The latter two are order-of-magnitude proxies,
not validated partially ionized collision coefficients. Bounded stellar
controls vary conduction and residual heat to assess their structural effect.
[Transport assessment](results/convective_transport_integration_sept26_v1.json).

Any region boundary needing microscopic species exchange rejects explicitly;
there is no silent zero-flux substitution in a radiative layer. The source
contains the hot screened transport but does not yet select it for such a
boundary. Initial D exhaustion is not a criterion for losing convection or
requiring finite mixing. Finite convection remains available in the shared
engine for phases where burning and mixing times become comparable.

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
The nonlinear composition tolerance and global inventory budget are separate:
the selected control solves to **1e-15** and checks inventories to **1e-14**.
The tighter solve resolves a marginal budget failure without loosening the
inventory or energy acceptance criteria.

`history.jsonl` retains every accepted model's scalar diagnostics;
`attempts.jsonl` retains acceptance and conservation checks. Only the initial,
latest and final full checkpoints are retained by default. Checkpoints include
CN/D abundances, the luminosity convention and, when applicable, material-heat
rates. Exact restarts require the same configuration, executable and all
referenced EOS/opacity source files. Run duration and compute/step budgets may
change. `report.json` distinguishes the requested age from a planned stop or a
physical/source-domain failure. A bounded run finishing is not a claim that the
full evolutionary track has been completed.
