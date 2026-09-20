# Matching the calculation to the question

The objective is the life of an initially 0.1-solar-mass star, including
formation, hydrogen burning, helium-remnant cooling past 100 K, and its fate
under explicit nucleon-decay scenarios. Use the simplest physical model that
can resolve the question being asked. Increase detail when it can change the
answer; carry uncertainty when more detail cannot resolve it.

Elapsed time alone does not decide accuracy. A close encounter at 10^14 yr
still lasts a few stellar dynamical times and must conserve energy accurately
enough to measure the transferred energy. Conversely, the uncertain frequency
of such encounters does not warrant a finely specified future Galaxy. Keep
the accuracy of a calculation conditional on its inputs separate from our
confidence that those inputs describe nature.

## Progression of approximation

| Phase or question | Appropriate calculation | What to report | When to add detail |
| --- | --- | --- | --- |
| Formation and contraction before the initial Ember model | Published low-mass formation and pre-main-sequence models; vary plausible starting entropy and composition | Context and uncertainty in the starting model, separate from Ember's age clock | If starting conditions leave a significant difference during hydrogen burning |
| Hydrogen burning and the turn toward cooling | Resolved 1-D structure, energy and abundance evolution; wavelength-dependent atmospheric boundary; supported thermodynamic and opacity tables | Track, profiles, event ages, and measured numerical and physical sensitivities | Where a missing opacity, equation of state, reaction or transport process changes the track |
| Helium-remnant cooling through 1000, 500 and 100 K | Resolved structure while transport, residual burning, diffusion or phase changes control cooling; then a reduced thermal model calibrated in an overlap interval | Ages conditional on surface composition, stored heat and phase physics; surface and core temperatures separately | If reduction loses an important thermal reservoir or changes a cooling milestone appreciably |
| Very cold remnant in an evolving environment | A few isolated, retained and ejected histories; event rates and distributions coupled to thermal response | Heating probabilities, deposited energy, cooling after events, ranges over environments | If deposition depth, a threshold, or rare-event tails decide the answer |
| Unbound white-dwarf encounters | Analytic orbital and tidal estimates plus an SPH survey; fluid and solid assumptions explicit | Transfer, retained heat or bounds, capture/stripping classification and numerical sensitivity | For cases that decide whether heating matters; particle count follows the signal and available compute |
| Pycnuclear reactions and uncertain dark matter | Competing-rate estimates, selected structure calculations and parameter ranges; a zero-dark-matter-heating case | Whether a process competes with cooling, encounters or decay under stated assumptions | If a plausible parameter range changes the qualitative fate or composition |
| Proton/nucleon decay | Lifetime and branching scenarios, mass and energy balances, simple hydrostatic or mass-radius models, coarse composition groups initially | Conditional phase order, times in units of lifetime, heating, loss of degeneracy and transparency; stable-baryon alternative | If daughter nuclei, energy escape or material phases materially change the order or retained heat |
| Dispersal and last particles | Binding estimates and stochastic particle survival | End of the bound object separately from the distribution of last decay times | Do not continue stellar structure or a blackbody atmosphere after their assumptions fail |

These are defaults, not mandatory reductions on a fixed date. A cold surface
does not imply an equally cold core. Grains and crystallization depend on local
composition and thermodynamic state; evaluate their importance where they occur.

## Decide effort from sensitivity

Before a new calculation, state the question, quantity needed, dominant
uncertain inputs, cheapest useful method, and result that would justify more
work. Start with conservation laws and order-of-magnitude estimates. Survey
parameters before refining one realization.

Compare numerical error with the signal and the decision being made. A
percent-level error in total binding energy can invalidate a much smaller
tidal-heating signal, even if the occurrence rate is uncertain by orders of
magnitude. Report an upper limit or an unresolved signal in that case.
Do not reinterpret numerical dissipation as physical heat.

Where sensitivity can be estimated, aim for numerical uncertainty below about
one tenth of the physical variation relevant to the conclusion. This is an
effort guide, not permission to relax existing Ember acceptance checks. For a
detection or phase boundary, demonstrate that the classification survives the
relevant controls. Do not force a precise answer that the controls cannot
distinguish within the survey budget.

Additional compute should be capable of changing an event's existence, order,
approximate age range, deposited energy or fate classification. Four significant
figures are a maximum in writing, not a target; speculative estimates usually
warrant one or two. Small calculations are not inherently better: use the
resolution the question needs while preventing a secondary survey from
dominating the project's compute or storage.

## Joining models and writing results

When changing approximations, compare descriptions in an overlap region and
transfer mass, composition, stored energy and relevant angular momentum
consistently. Record what was averaged away. A reduced cooling model must
reproduce the resolved model there before extrapolation. For encounter heating,
retain a coarse deposition profile if it affects subsequent cooling; a full
particle history need not be retained.

Label calculated track points, literature comparisons, physical estimates and
conditional scenarios distinctly. Do not draw a continuous computed track into
unintegrated future scenarios. A reference proton lifetime is not a measured
event age. Changing it can alter temperature, chemistry and energy escape as
well as rescale time. Maintain a short account of assumptions and their effects.
Numerical precision cannot reduce uncertainty in an unmeasured interaction or
a future environment. Present robust and conditional outcomes clearly without
recounting the research history.

## Compute and storage

Reuse equilibria and source calculations. Set survey compute and storage
budgets and stop rules; select particle count after a benchmark, not from an
arbitrary low ceiling. Retain input recipes, source changes, hashes, compact
diagnostics, uncertainty estimates and representative states. Save full time
series only to answer a specific physical question. Never delete another task's
data as part of cleanup.

The [Fable assignment](research/fable/ASSIGNMENT.md) applies this remit to the
encounter survey. The [cold-remnant plan](COLD_REMNANT.md) and
[very cold physics plan](ULTRACOLD_PHYSICS.md) describe the stellar work.
