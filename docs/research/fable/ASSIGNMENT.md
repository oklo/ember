# Fable: literature and simulations of unbound white-dwarf encounters

Updated 2026-09-11. The user explicitly requests simulations IN ADDITION to
literature and physics work. This replaces the literature-only assignment and
its prohibition on SPH surveys. Proceed autonomously. The aim is to advance
understanding and use available compute productively, without allowing this
survey to displace the main Ember trajectory or accumulate large outputs.
There is no arbitrary low particle-count ceiling.

Read [COORDINATION.md](COORDINATION.md) and PRIMARY_NOTES.md first. The protocol
defines numbered messages, acknowledgments, independent tasks, process-aware
resource reservations, pause handoffs and immutable results for asynchronous
review. Follow it when starting, pausing and resuming.

Read HANDOFF.md, docs/COLD_REMNANT.md, docs/ULTRACOLD_PHYSICS.md and
[the approximation remit](../../APPROXIMATION_REMIT.md). The primary agent owns
Ember's stellar evolution, shared code and working paper. Your independent
assignment includes SPH adaptation, controls, simulations, literature and analysis.

The project follows an initially 0.1-solar-mass star through helium-remnant
cooling past 100 K and disappearance under explicit nucleon-decay scenarios.
Environmental heating and stable-baryon alternatives are part of this account.
Neither a fixed proton lifetime nor a known WIMP particle is assumed.

## Scientific question

How much energy, and specifically irreversible retained heat, can a close
WD–WD encounter deposit if the pair REMAINS UNBOUND after passage? Focus on
pericenter just above contact, then ask whether probabilities and the remnant's
cooling response make this relevant around 10^14 yr and beyond. The remote
epoch changes environmental uncertainty; it does not remove the need to
measure the encounter's small energy transfer reliably.

Test the user's suspicion that existing work emphasizes bound systems,
captures, mergers and detonations. Do not claim an unexplored population from
a few examples. Distinguish initial/final orbital energy, WD–WD/WD–black-hole
encounters, fluid/solid objects, and heat/mode energy.

## Local implementation inspected by the primary agent

The curated setup is `/Users/greglaughlin/Projects/mars-earth-collision`.
Read README.md and docs/REPRODUCIBILITY.md, then:

- src/assemble_settled_impact.py: loading relaxed bodies, centering, spin and
  particle IDs. Adapt the orbit to pericenter and speed at infinity.
- src/analyze_body_snapshot.py: basic structure and velocity checks only;
  it does NOT measure a WD energy budget or irreversible heat.
- src/prepare_swift_run.py and configs/mars_earth_grazing_settled_smoke.yml:
  staging, sparse outputs, gravity and hydrodynamics settings.

SWIFT is in `/Users/greglaughlin/Projects/earth-mars-swift`:

- `swift-impact/swift`: energy-evolving planetary SPH with self-gravity.
  config.h selects EOS_PLANETARY; PLANETARY_FIXED_ENTROPY is undefined.
- `swift/swift`: fixed-entropy settling build. Do NOT use it to measure heating.
- `swift-impact/src/equation_of_state/planetary/equation_of_state.h`, near
  line 2680: ten custom SESAME-format table slots. A suitable WD table may
  reuse the impact executable; validate interpolation and thermodynamics.
- An ideal-gas gamma=5/3 implementation exists, requiring a separately
  configured build. Its temperature is NOT a WD temperature.
- The documented builds use CPU threads; MPI and vectorization are disabled.
  No Metal/GPU implementation was identified in this inspection. Measure a
  pilot before contemplating performance changes; do not start a GPU rewrite.

The iron/rock ANEOS bodies are planetary models, not scaled WDs. Earth–Mars
softening lengths, units, box sizes and output cadence cannot be copied blindly.
The source has artificial viscosity and predominantly float thermodynamic
particle fields: measure the error floor for small thermal residuals. The old
visualization series totals about 116 GB. Do not copy it or recreate its dense
movie cadence. No SWIFT job was active at the primary agent's 14:34 UTC process
check; inspect again before launching. The sibling projects remain read-only.

## Physical model sufficient for the question

Begin with nonrotating spherical hydrostatic degenerate fluid bodies. A
nonrelativistic gamma=5/3, n=1.5 polytrope is acceptable for a labelled low-mass
hydrodynamic pilot, deriving its mass-radius scale from electron degeneracy and
composition. Do not interpret its ideal-gas energy or temperature as the real
WD thermal reservoir. A pilot can measure orbital/mode transfer without claiming
a physical temperature or heat capacity.

For physical heating, use consistent cold degenerate-electron pressure and
energy plus an interpretable thermal contribution, or a documented WD EOS
covering encountered states. Check adiabatic response, sound speed and isolated
hydrostatic stability. Use relativistic electron corrections where needed.
Do not obtain a claimed tiny thermal signal from unresolved subtraction of
large degeneracy energies. A thermal-variable or entropy diagnostic may be
needed; choose the smallest adequate implementation.

Start with a 0.1 + 0.1 solar-mass helium pair as a controlled reference. Extend
to a 0.1 + 0.6 solar-mass He/CO pair with an adequate companion structure if
useful and tractable. Reuse equilibria between encounters. The current Ember
endpoint is a hot hydrogen-burning star, NOT a cold-WD initial condition.
Record the source and assumptions of each mass, radius and density profile.

Treat fluid simulations as fluid simulations. Estimate whether solid-core
rigidity or yielding can change the response. Bound the effect or identify the
limitation; detailed elastic/plastic SPH is not required for the initial survey.
Do not call a fluid result a rigorous upper/lower bound without justification.
Radiation and nuclear reactions can be omitted during weak passages if timescale
estimates support that. Flag cases where shocks or destruction invalidate this
approximation and end that branch without following a merger or detonation.

## Adaptive survey

1. Inspect jobs, source revisions and dependencies; record exact inputs and
   executable hashes. Benchmark a small run and choose particle count, thread
   allocation and output cadence from cost and the measurable signal.
2. Prepare reusable isolated bodies. After any settling, run controls using
   the SAME energy-evolving settings and duration as encounters. Measure radius
   drift, residual motion, artificial heating and conservation.
3. Define q = r_peri/(R1+R2) and w = v_infinity/v_escape,contact, where
   v_escape,contact = sqrt(2 G (M1+M2)/(R1+R2)). State how relaxed radii are
   measured and whether an atmosphere is resolved. Distinguish q from impact
   parameter at infinity and inverse penetration factor.
4. Candidate first grid: q = 1.05, 1.2, 1.5, 2.0 and w = 0.02, 0.1, 0.3.
   Use pilots to choose an informative subset or refine a boundary. Roughly
   50000–100000 total particles is a reasonable starting estimate, not a rule;
   200000 or more can be appropriate for selected controls if the benchmark
   and scientific gain justify it. Quote actual counts: the sampler adjusts
   requested counts. Do not sacrifice a measurable signal merely to keep N low.
5. Initialize an inbound hyperbola from exact finite-separation energy and
   angular momentum. Check a point-mass control and starting-separation
   sensitivity. Test particle orientation/sampling on a representative case.
6. Follow departure until outgoing orbital energy and transfer stabilize over
   more than one separation. Identify each self-bound remnant and separately
   account for stripped/transferred material. Initial positive energy does not
   establish final escape. Classify unbound survival, capture, stripping,
   contact, or unresolved outcome. Avoid extending captures into long mergers.
7. Repeat informative cases at higher resolution and with changed viscosity,
   gravity and time accuracy. Choose controls to test the conclusion, not an
   automatic large resolution ladder. Compare weak encounters with analytic
   tidal excitation and a distant passage. Adapt the next batch to the results.

Do not launch the whole grid before controls show what is measurable. If heat
is below numerical background, useful mode-energy bounds or capture boundaries
may still be possible. Prefer those over an unbounded quest for a tiny residual.

## Energy and interpretation

Track mass and particle IDs, total energy, angular/linear momentum, outgoing
orbital energy, spin, residual internal motion, and a separately defined thermal
diagnostic. Distinguish reversible compression from entropy production and
persistent modes. SPH viscosity damping is not by itself a physical damping law.

Compare signals with isolated drift, integration error, particle sampling,
viscosity and resolution changes. Use the largest relevant control variation
as a conservative numerical scale, not an assumed independent Gaussian error.
Where practical keep that scale below 10% of a claimed heat signal. Otherwise
report a range, upper limit or unresolved result. Conservation relative only
to total binding energy is insufficient. Bracket later conversion of measured
mode energy to heat and state the damping and deposition assumptions.

Couple transfer or heat bounds to encounter rates and cooling response. Use
plausible retained/ejected environments and distributions of waiting times and
heat pulses. Do not prescribe one precisely known Galaxy at 10^14 yr. Establish
when a mean heating rate represents intermittent events poorly. A coarse radial
energy-deposition profile is enough unless finer structure changes cooling.

## Compute and data management

Use spare resources autonomously. Benchmark throughput for one larger threaded
run versus concurrent independent runs; select the arrangement that produces
useful results fastest. Keep the main stellar continuation and atmosphere
calculations first in line. Recheck jobs before each batch; reduce new SPH
concurrency when those jobs are active. On this 14-core, 36-GB machine, use up
to about 12 CPU threads across SPH jobs when otherwise idle, leaving room for
the interface. This is not a mandate to saturate threads if throughput falls.
Keep aggregate SPH resident memory below about 8 GB initially.

Work in bounded batches, initially up to 24 aggregate CPU hours of preparation,
controls and encounters. Publish results and resource usage after each batch,
then autonomously continue the next scientifically justified batch within the
same resource discipline. Do not wait for a user nudge or approval at every
boundary. Reassess with the primary agent via files if repeated batches fail
to improve the conclusion; do useful literature/analytic work in the meantime.
Use measured run times to avoid scheduling many long cases without feedback.
Particle count follows measured cost and value; no fixed low ceiling applies.

Cap generated scratch data at 2 GB across active runs, including logs,
snapshots and restarts. Inventory private source/build copies separately and
keep the TOTAL private footprint below 4 GB; do not count existing read-only
shared dependencies as newly generated data. Estimate snapshot bytes first and
enforce limits in the launcher. Reduce output frequency or stop your own jobs
cleanly before breaching the cap. Aim for retained research artifacts below
100 MB, preferably much less. These limits constrain output, not particle count.

Use frequent compact scalar diagnostics and normally no more than eight full
snapshots per encounter, plus at most one rolling restart. Retain initial,
near-passage and separated outgoing states. Add on-the-fly diagnostics when
sparse snapshots cannot measure the physics. Do not create movies or frame
series. Compress small arrays and preserve only necessary particle fields.
After checking analysis and reproduction recipes, delete only your own
redundant intermediate outputs and record that in a manifest. Never delete
the user's planetary data, shared inputs or another task's files.

## Deliverables and continuing work

- STATUS.md at start and milestones: UTC time, actual PID/session, command,
  scratch path, compute/storage totals, finding and next action.
- FINDINGS.md: checked primary-source literature table, simulation results,
  assumptions, numerical controls, limits and recommendations. Literature
  columns: masses/compositions, initial/final energy, pericenter convention,
  solid/fluid model, numerical method and measured energy quantity.
- Runnable survey/analysis recipes, seeds, configs, EOS definition, exact
  source changes/hashes, compact run receipts and CSV/JSON results with errors.
- A few PDF/PNG plots and representative states within the storage cap.
- Encounter probabilities and thermal-response estimates, including solid-body
  uncertainty and retained/ejected environments; specify what further effort
  could change the conclusion. Publish partial findings as they become useful.

Keep the literature work moving alongside simulations and use it to choose
controls and interpret results. Starting references, not an exhaustive search:

- Loren-Aguilar et al. 2010: https://academic.oup.com/mnras/article/406/4/2749/1022010
  (explicitly initially bound systems).
- Aznar-Siguan et al. 2013: https://academic.oup.com/mnras/article/434/3/2539/1042364
  (check actual orbital energies).
- Fuller & Lai 2013: https://arxiv.org/abs/1211.0624 (bound helium-WD tides).
- Schaller et al., SWIFT: https://arxiv.org/abs/2305.13380 .
- Adams & Laughlin 1997: https://arxiv.org/abs/astro-ph/9701131 .

## Coordination

Write only within docs/research/fable/ and scratch paths /tmp/ember-fable-*.
Adapt/build SPH there and retain source changes here. Do not edit the sibling
projects, shared Ember code/tables, main HANDOFF.md or paper. Do not switch
branches, stage, commit or push this shared checkout. The primary agent reviews
and integrates results. Read PRIMARY_NOTES.md at work boundaries; put requests
in QUESTIONS.md. Each agent owns its own communication files; the primary owns
ASSIGNMENT.md and PRIMARY_NOTES.md. Do not overwrite another agent's results.

Shared files do not themselves schedule either agent. Do not claim jobs are
running before launch. Use primary sources and exact URLs, plain writing and at
most four significant figures in prose; preserve numerical precision in data.
