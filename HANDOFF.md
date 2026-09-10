# Completed central-core models and updated PDF — 2026-09-10 17:06 UTC

The user requested an updated PDF. Today's paper was overwritten in place,
now ten pages, with the completed 3.560-trillion-year calculations, LBA97
comparisons, a new convection/hydrogen-profile figure, and numerical grain
results. Four significant figures and faint F77 curves remain in force.

Both stellar jobs previously listed as active are COMPLETE. Process inspection
at 16:43 UTC found no surviving stellar, atmosphere or grain-source calculations.
No new stellar run was launched during this PDF update. The autonomous goal
of actual evolution to a 100 K helium remnant remains unfinished.

- Fresh 512-point v2: `out/evolution-cold-remnant-x015-transition-512-3560gyr-v2.*`.
  Reached 3.560 trillion years, 2354 history states, 198 rejected attempts.
  Central X=0.1518, surface X=0.1747, Teff=3240 K, L=0.002025 Lsun,
  convective mass=0.6277. Native checkpoint intact. Executable SHA256
  `7a4a808e801f1ce523f5c692d54583420a629c3139a64a36843e85f98c8849e0`.
- 1024-point exact continuation v1: same final age, 532 segment states,
  central X=0.1491, Teff=3241 K, L=0.002027 Lsun, convective mass=0.6223.
  No new checkpoint was written, as intended. Older checkpoint/binary preserved.
- BOTH final centers are now stable against convection. Actual Ledoux region
  probe gives stable central mass fractions 37.23% and 37.77%, reaching about
  45% of radius. The old 3.548-trillion-year 512 trial still had a convective
  center; its age, X, He3, radius and luminosity reappear exactly in the new
  fresh history. Central-convection disappearance is bracketed only between
  3.548 and 3.560 trillion years, not yet accurately located.
- At equal final age, 1024/512 luminosity differs by +0.07813%, temperature by
  +0.8654 K, central hydrogen by -1.752%, central density by +1.108%, total
  hydrogen mass by -0.1707%. Further mesh and time-step controls are required.
  Hydrogen burning supplies almost all luminosity; no cooling claim is made.

`scripts/summarize_transition_runs.py --mixing-probe
/tmp/ember-evolution-mixing-probe-v1` checks source outputs and receipts,
finite histories and normalized profiles, then evaluates actual Ledoux regions
and exports the new small data files. Report:
`docs/results/evolution_transition_3560gyr_v1.json`. All captured physical-data
hashes agree between runs; neither receipt reports changes during its run.
All five figures (PDF/PNG) and two LBA97 CSVs reproduce byte for byte from
published data in `/tmp/ember-paper-3560-rebuild-v1`. Hash-checked local archive
recovery reproduces the new 2354-state CSV exactly. Fourteen recovery archives
have verified raw/compressed hashes and remain local. Earlier files are intact.
No C++ changes were made for this paper update; the prior 32-suite pass remains
applicable, with its original test date preserved in validation.json.

Grain atmosphere coupling is still unfinished. The new text quantifies the
existing fixed 2800 K alumina profile only; no additional grain physics was
installed in either stellar run. The previous entry's material/chemistry and
thermodynamic work remains necessary. Continue the scientific work after the
PDF publication; inspect processes again before launching calculations.

---

# Grain optics, LBA97 comparison and transition checks — 2026-09-10 15:35 UTC

The autonomous 0.1-solar-mass evolution to a correctly calculated 100 K helium
remnant remains active and unfinished. The user additionally asked to add grain
opacities and to make LBA97, inferred from its paper where necessary, the main
working-paper comparison. Preserve all earlier writing preferences: at most
four significant figures in prose and graph labels, plain language, one paper
per day, and faint F77 curves with abrupt segments made more transparent.

The September 10 paper now has LBA97 HR and composition comparisons. It records
42 HR positions, 20 hydrogen positions and 15 helium-3 positions from Figure 1,
plus stated values. The original PDF and extracted image were inspected;
selected pixels and calibration are published, not a falsely precise original
history. LBA97 data and both PDF/PNG figures reproduce exactly from the small
published files. F77 remains supplementary and explicitly a recent
reconstruction with changed inputs. Both new figure pages and all nine pages
of extracted text were checked. The completed Ember curve ends at 3.40T; open
orange points identify a later saved trial at 3.548T. No cooling claim is made.
The README and scientific notes are updated with the paper.

The 120-atmosphere gas table is complete and independently checked:
`/tmp/ember-nongrey-exhaustion-warm-x015-v1.dat`, SHA256
`5ae51fad184894503119eb06a751fd8aedafdbcbb531d3683bf055123be8f924`.
Assembly specification and result are under data/atmosphere/sources and
`docs/results/nongrey_exhaustion_warm_x015_v1_audit.json`. No source-atmosphere
jobs remain active. The new hydrogen-poor profiles still need condensation
checks. This is warm gas coverage, not a cold or zero-hydrogen atmosphere.

The first revised 512-point trial stopped near3.548T at a SECOND old counter
check in `apps/evolution_checkpoint.hpp`: lifetime counters were still capped
at10000 accepted and101 rejected despite the evolution-loop correction.
The native checkpoint retains2061 accepted and101 rejected; its age is
3.548T, slightly before the last state printed before the write failed.
The 1024-point trial reached a nonconvective layer near3.545T and stopped at
the same checkpoint check. The top-level error JSON does not contain their
histories; preserve the logs, native checkpoints, receipts and copied binaries.
The corrected checkpoint parser accepts arbitrary representable nonnegative
lifetime counters and rejects negative/overflowing input. Execution limits
remain10000 accepted steps per invocation and100 consecutive rejections.
The expanded exact-restart test passes with15000/102 counters in a test-only
fixture; all32 CTest suites pass in84.76s. No research checkpoint was edited.
Report: `docs/results/evolution_checkpoint_counter_v1.json`.

Active calculations (inspect processes before launching anything):
- Corrected fresh512 run, controller77620/star77621, session82038,
  snapshot `/tmp/ember-cold-remnant-x015-transition-3560gyr-v2`,
  output `out/evolution-cold-remnant-x015-transition-512-3560gyr-v2.*`.
  Target3.56T, checkpoints enabled, same120-atmosphere/EOS/TOPS/plasma inputs.
- Exact1024 continuation, session11867, snapshot
  `/tmp/ember-cold-remnant-x015-transition-1024-3560gyr-continuation-v1`,
  output same descriptive stem under out. It uses the IDENTICAL older copied
  executable and untouched v1 native checkpoint, with no checkpoint output.
  This avoids the old writer limit while preserving all physics and restart
  identity; it cannot produce a new native checkpoint. The corrected fresh
  run provides future checkpoint support. The old1024 source receipt correctly
  records that working-tree sources changed while its copied binary ran.

Offline diagnostics locate the512 stable shell at roughly19–43% of radius,
with central convective mass3.823%, stable mass29.23%, outer convective
mass66.95%. Near enclosed mass21.10%, the required/adiabatic gradient ratio
falls from1.588 at3.30T through1.336 at3.40T to.9274 in the saved trial.
The center remains convective (ratio1.421). Rising T and falling opacity favor
radiation; conduction supplies4.144% of the local diffusive capacity.
`docs/results/evolution_convection_transition_v1.json` records the evidence.
This is not a converged transition age or a central radiative core. The user
asked about partial convection: heat shares vary continuously through MLT;
instability classification is binary and connected convective regions are
mixed instantaneously. Slow mixing exists but semiconvective heat is omitted.

Grain optics implemented and checked, NOT atmosphere-coupled yet:
- `scripts/prepare_grain_sources.py` pins LX-MIE and Optool and builds separate
  offline tools. Normal MIT Ember does not link the GPL LX-MIE adapter.
  Current source receipt `/tmp/ember-grain-sources-v3/prepared.json`.
- `grain_opacity.py`, `grain_mie_probe.cpp`, `generate_grain_opacity.py` provide
  separate absorption/scattering/g, explicit size and density, number-to-mass
  weighting, FastChem condensate mass fractions, no spectral extrapolation.
- Eight independent controls pass, plus Rayleigh, lossless/no-contrast,
  size weighting, zero-grain and invalid-input checks. Optool imposes a
  nonphysical absorption floor for lossless grains; compare scattering/g only
  there and check zero absorption analytically. Refine its g quadrature to1800
  angles without loosening tolerance. Final report
  `docs/results/grain_mie_v1_audit.json`, raw `/tmp/ember-grain-mie-audit-v5`.
- First fixed gas profile at2800K/g5.15/XH.7 contains only alumina. Explicit
  exploratory amorphous optical data and density3.2g/cm3; radii.01/.1/1micron.
  `docs/results/grain_alumina_profile_v1.json`, raw
  `/tmp/ember-grain-alumina-profile-v1`. Full gas/element/pressure closure
  rechecked; maximum condensed mass fraction.0001094. This does not include
  atmospheric feedback or updated gas opacity.
- Source atmosphere extends to.09micron but alumina data begin.2micron.
  Resolve that gap physically before coupling. Phase/temperature/porosity
  require comparisons: amorphous alumina anneals, so it is not automatically
  the equilibrium grain phase. A coupled gas-depletion2800K profile already
  adds CaTiO3 and gehlenite; LX-MIE includes the former, not the latter.
  Do not silently ignore missing grains or replace their optical properties.
- Next: complete optical-material coverage, couple absorption/emission and
  angular scattering with T/P responses, solve the atmosphere again, check
  flux/depth/material support and consistent condensate enthalpy. Existing
  importer guards remain unchanged. See `docs/GRAINS.md`.

Unfinished helium coexistence probes remain untracked development files;
no phase prescription or liquid/solid EOS was installed. The rejected electron
screening cache was fully reverted after a one-ULP output difference. Its
negative result is retained. Continue the scientific work, not just the handoff.

---

# Paper comparisons, writing preferences and transition trial — 2026-09-10 13:10 UTC

The user asks for autonomous evolution of the 0.1 solar-mass star to a correctly
calculated 100 K helium remnant. That work remains active and unfinished.

The September 10 paper has been rewritten for a reader arriving fresh, with
F77 AJR and Ferguson curves against both age and hydrogen abundance. Use thin,
translucent F77 curves; abrupt changes are drawn more faintly, with all saved
data retained. F77 luminosity and radius now use Ember's nominal solar units.
The luminosity difference at the present hydrogen abundance is about 76.4%;
the earlier 74.9% comparison used each code's different solar luminosity unit.
Both conventions are explicit in the matched-state report.

Persistent user preferences: at most FOUR significant figures in writing,
tables and figure labels; retain full precision in machine data and exact
identifiers. Avoid internal workflow jargon. Explain what was calculated and
what it means. Update the current paper in place during the day; keep one dated
paper per day. The prior GitHub publication is commit36819fb; this revision
adds the comparison figures and reader-facing rewrite after checking the PDF
and reproducing both figures exactly from the published data.

The fresh stellar trial with the 108-atmosphere table stopped before its
requested 3.60-trillion-year age. Its last accepted model is near3.544 trillion
years, XH=.1751, Teff3216 K, convective mass fraction.7638. It first lost full
convection near3.543 trillion years. It accepted1988 steps and rejected101.
The actual stop is the cumulative100-rejection guard in apps/evolve.cpp,
not a demonstrated physical-domain failure or a successful completed track.
Its checkpoint preserves the cumulative rejection counter, so blindly restarting
will repeat the stop. Diagnose the numerical transition and execution limits
without changing old checkpoint identities. Outputs and exact copied binary
remain at the paths in the preceding entry. PID72500 is no longer running.
The completed, separately checked reference remains3.40 trillion years.

Both new independent atmospheres at XH=.125/.175 are complete. Against the
108-atmosphere table, matching-temperature differences are1.539%/1.833% and
pressure differences are-.02355%/-1.068%. Full source reassembly, input hashes
and interior checks pass. Report:
`docs/results/nongrey_x010_quarter_points_v1_audit.json`.

The XH=.15 plan still runs in `/tmp/ember-nongrey-x015-refined-warm-v1`.
Two3200 K/logg5.4 initial models stopped at too small an optical depth to
include the matching point; a converged initial guess reached only about88.
Both failed atmospheres now pass after a separate two-model retry with
initial_bottom_tau1200 instead of600:
`data/atmosphere/sources/nongrey_x015_high_gravity_retry_plan_specification.json`,
work `/tmp/ember-nongrey-x015-high-gravity-retry-v1`, session66285. Existing
jobs were inspected before launch. No star inputs or accepted source checks
were changed. Finish and independently assess the refined atmosphere family.

A five-second sample of the owned stellar process is preserved at
`/tmp/ember-stellar-diagnostic-20260910.sample.txt`. Repeated nuclear-screening
electron inversions appear costly; no optimization from that sample has yet
been implemented. Current C++ source is unchanged since the prior32-suite pass.

---

# Autonomous continuation and September 10 publication — 2026-09-10 12:44 UTC

This entry supersedes older live-status paragraphs. The completed reference is
still **3.40T**, fully convective, gas-only. No hydrogen exhaustion or cooling
endpoint has been reached. The user reaffirmed autonomous evolution through
Teff=100 K, then explicitly requested replacing the September 9 working paper
and updating GitHub/README. The calculations continue during publication.

- Initial process inspection found all old atmosphere controls complete and
  no surviving stellar/source/condensate/F77 jobs. No existing job was signalled.
- Finished the XH=.15 interpolation separation audit:
  `docs/results/nongrey_x015_interpolation_separation_v1_audit.json`.
  Composition-only T +1.652%, Pgas -.7842%; temperature/gravity
  T +.4743%, Pgas -.002543%. Full-grid +2.134% T and -.7868% Pgas.
  Log decomposition closes below 9e-16. Hydrogen spacing dominates.
- **Active fresh diagnostic star:** copied same 24a6ec3d... executable,
  snapshot `/tmp/ember-cold-remnant-x010-grid-diagnostic-3600gyr-v1`,
  output `out/evolution-cold-remnant-x010-grid-diagnostic-512-3600gyr-v1.json`,
  checkpoint same stem `.restart`, session59761, PID72500/controller72499.
  Target3.60T, 512 points, same refined EOS/TOPS v4 and plasma losses,
  selects `/tmp/ember-nongrey-exhaustion-warm-x010-v2.dat`. This candidate's
  interpolation is under refinement; the run is diagnostic, not a replacement
  for the completed-reference result or a validated transition. Fresh start
  preserves native restart identity. At12:40 it had reached about2.00T.
- **Active XH=.15 refinement:** new material specification
  `data/atmosphere/sources/nongrey_x015_refined_warm_specification.json`;
  two opacity planes in `/tmp/ember-nongrey-x015-refined-opacity-v1`,
  controller72491/session96163, 2 workers. Atmosphere plan
  `nongrey_x015_refined_warm_plan_specification.json` waits for those tables,
  then computes12 canonical3200/3400K, g4.9/5.15/5.4, X3=0/.12 states.
  Work `/tmp/ember-nongrey-x015-refined-warm-v1`, controller72537/session60968.
- **Independent quarter points:** XH=.125/.175, X3=.005,3300K/g5.1,
  specifications `nongrey_x0125_x0175_heldout{,_plan}_specification.json`.
  Opacity work `/tmp/ember-nongrey-x0125-x0175-heldout-opacity-v1`,
  controller72496/session7003; atmosphere work
  `/tmp/ember-nongrey-x0125-x0175-heldout-atmosphere-v1`,
  controller72540/session51224. Two workers each; source logs have work-name
  `.log` alongside each work directory. Both material families had13/14
  completed isotherms per plane at12:42. No automatic installation is attached.
- New native source diagnostic `scripts/audit_nongrey_abundance_scale.py`
  completed, report `docs/results/nongrey_abundance_scale_v1_audit.json`.
  At4017K, XH=.15/X3=.005, common abundance scales1/.15/3 preserve mass
  fractions to roundoff and electron densities within1.14e-6 relative, but
  sampled absorption coefficients change by up to factor27.3 (scale.15).
  Simple arbitrary abundance rescaling is therefore not verified as a route
  to true zero H. This is one spectral isotherm, not a stellar-boundary error.
  Raw controls `/tmp/ember-nongrey-abundance-scale-v1`; all completed.

Publication files: new8-page paper and joined1861-state0–3.40T history under
`docs/reports/2026-09-10`. September9 published edition removed from the
current tree, with its ignored local `artifacts/` retained, including its old
recovery manifest. New recovery archives contain both current histories and
receipts, the3.40T checkpoint and exact copied executable. Bulk tables and raw
archives remain local. README, roadmap, forward/survey summaries, reproduction
instructions, next-session prompt and cold-EOS status now reflect September10.

Validation: separate Release build passes32/32 CTest suites. Additional
condensate archive tests exposed exact-dictionary mismatches after two new
range diagnostics were added. `source_diagnostics_match` accepts only the
complete legacy or current schema and requires exact equality for every
recorded field; all current physical/source/hash checks remain. After the fix,
13/13 condensate tests and26/26 non-grey import tests pass (targeted CTest
retest1.38s). Paper compiles with no warnings; four representative rendered
pages and all extracted text checked. CSV-only figure rebuilding is byte-
identical in PNG, with all recovery hashes verified. See adjacent validation.

Continue the still-active stellar/source jobs, assemble and independently audit
the refined atmosphere only after canonical completion, then compare evolution.
Do not stop the remnant objective at a documentation handoff. Dense matter,
diffusion, condensates, zero-H boundaries and cooling remain substantial work.

---

# Continued calculations and F77 comparison — 2026-09-10 11:55 UTC

**This checkpoint supersedes the live-status entries below. Furthest completed
Ember evolution is now3.40 trillion years. It is fully convective and gas-only;
no hydrogen-exhaustion or white-dwarf/cooling endpoint has been reached.**

- Inspected processes before new launches; existing productive controllers
  were reused through completion. No parent-owned job was signalled here.
  At11:47 the four old cycling TLUSTY jobs, their controllers and the slow
  direct canonical comparison were no longer listed. Their old work trees
  lack new canonical finals/terminal manifests; do not infer successful
  completion or claim this turn stopped them. Idle condensate waiter27809
  remains. Its missing install receipt was not written.
- Fresh512-zone forward run0->3.30T completed,1822 accepted/3 rejected,
  2351 CPU seconds and4088 UTC seconds. Native output/checkpoint:
  `out/evolution-cold-remnant-forward-512-3300gyr-v1.{json,restart}`.
  Snapshot `/tmp/ember-cold-remnant-forward-3300gyr-v1`, executable SHA
  `24a6ec3de0d06bf1e58e4ff4527fa4b38f4b02774be1703cba4ba371d7f15237`.
  Physical inputs did not change. Receipt records source-tree edits in the
  side conversation; the copied executable was independent of those edits.
- Exact same-executable/physics restart3.30->3.40T also **completed**,
  session44958,97.00 CPU seconds. Output/checkpoint
  `out/evolution-cold-remnant-forward-512-3400gyr-v1.{json,restart}`;
  copied binary and receipt `/tmp/ember-cold-remnant-forward-3400gyr-v1`.
  XH=.2034, X3=.0009703, R=.1459Rsun,
  L=.002016Lsun, Teff3202K, Tc8.031MK,
  rhoc258.3g/cm3. All inputs and incoming restart hashes unchanged.
  Both endpoint reports are in `docs/results/remnant_endpoints_*forward_v1.json`.
  **No stellar process remains active.** The selected96-node atmosphere
  still ends atXH=.2; do not launch blindly beyond that support or bypass
  exact-restart identity to insert a different atmosphere.
- All12 XH=.1 warm source nodes completed. Candidate108-node gas table:
  `/tmp/ember-nongrey-exhaustion-warm-x010-v2.dat`, SHA
  `ca92b939ae386e7d2edc69a5053c807ed4a92fca846963d1c4e717d08609b255`.
  Explicit36-continuation recipe
  `data/atmosphere/sources/nongrey_exhaustion_warm_x010_v2_assembly_specification.json`.
  Source reassembly/EOS/runtime audit completed, all108 nodes supported;
  density difference-.01082..+1.991%. Prior72-node values retained exactly.
  Additional112 queries preserve all96 prior warm nodes plus each old cell
  interior, including derivatives, exactly. Reports
  `nongrey_exhaustion_warm_x010_v2{,_retained_warm_v1}_audit.json`.
  **Candidate remains in/tmp and is not installed or selected by a star.**
- Two new3200/g5.4 depth controls completed under
  `/tmp/ember-nongrey-x010-cool-high-gravity-depth-v1`, session30264 done.
  Actual depths165/339 and170/350; matching-state maxima7.13e-6/1.93e-5 pass.
  Reports `nongrey_x010_cool_high_gravity_y{000,120}_depth_v1_audit.json`.
- Independent XH=.15/X3=.005/3300/g5.1 comparison of the108-node candidate:
  matching T+2.134%, Pgas-.7868%. AtXH=.25/3100/g5.1 the prior+.6904%T
  and-.9341%Pgas differences are unchanged. Local diagnostics, not global
  physical or lifetime errors. The largerXH=.15 discrepancy needs resolution.
- Completed four-sourceXH=.25 interpolation separation audit:
  composition-only+0.3359%T/-.2012%Pgas; T/gravity+0.3533%T/-.7344%Pgas.
  Logarithmic decomposition closes to8.9e-16. New reusable script
  `scripts/audit_nongrey_interpolation_separation.py`, report
  `docs/results/nongrey_warm_interpolation_separation_v1_audit.json`.
- **Active new work:** four3300/g5.1 controls atXH=.1/.2 andX3=0/.12,
  `/tmp/ember-nongrey-x015-interpolation-separation-v1`, controller71615,
  session34141, two workers; same-name.log. Plan
  `data/atmosphere/sources/nongrey_x015_interpolation_separation_specification.json`.
  At11:54 bothXH=.1 canonical finals completed; twoXH=.2 cases remain.
  Run the separation audit with the108-node candidate, these four controls
  and `/tmp/ember-nongrey-x015-heldout-atmosphere-v1/atmosphere` after completion.
  No automatic installation or stellar continuation is attached to this plan.

The user asked about computation cost, reusable atmospheres, F77 core onset
and why F77 differs from Ember. Findings and receipts:

- `docs/COMPUTATIONAL_COST.md`: tables are reusable; source extension and
  verification are separate from stellar integration. Do not attribute all
  prior machine load to the single stellar process. Prior cache benchmarks
  gave~22% combined CPU reduction with byte-identical outputs. No matched
  F77/Ember whole-star performance attribution has been made.
- Local F77 source matches the retained baseline hash exactly. Separate
  `/tmp/ember-f77-core-onset-v1` copies add only a diagnostic print. Unmodified
  controls reproduce original output bytes for both the6000-model Ferguson
  check and full AJR cooling run. AJR center first nonconvective at3.224T,
  XH=.01756; gradient unity crossing one model earlier. Current Ferguson
  option first flags it at2.946T,XH1.064e-5; later exits2 during cooling.
  Do not equate these with the1997 paper's5.742T onset. Diagnostic patch,
  exact decks, hashes and brackets in `f77_core_onset_comparison_v1.json`;
  reusable `scripts/summarize_f77_core_onset.py`. All F77 jobs are done.
- `scripts/evolution_physics_probe.cpp` now accepts explicit EOS/opacity paths.
  Offline probe `/tmp/ember-evolution-physics-exhaustion-probe-v1` at current
  inputs: minimum diffusive/adiabatic gradient1.588->1.336 from3.30->3.40T,
  at massfraction~.211; central3.707->3.090. Approaching loss of full
  convection does not prove imminent central-core onset. Nodal audit, not
  the runtime face criterion. `evolution_convection_margin_3400gyr_v1.json`.
- At matchedXH~.203, AJR F77: age2.507T, L.003525Lsun, Teff3698K,
  R.1452Rsun,Tc8.000MK; Ember:3.40T,L.002016,Teff3202,R.1459,Tc8.031MK.
  F77 He3=.001516 vs Ember.0009703. Swapping only the F77 low-T opacity
  to Ferguson raisesL54.5% at this H abundance. Identical central-state
  pp-heating comparisons differ~9%, not the full75% surface-luminosity gap.
  Atmosphere/opacity and composition evolution are leading explanations;
  EOS, nuclear screening and time-coupling changes remain confounded.
  `f77_ember_matched_hydrogen_v1.json` retains comparisons, rate-probe recipe
  and primary literature context. No claim of uniquely isolated causation
  or of Ember being a strict physics superset/validated truth.

Validation this turn:25 atmosphere import tests pass; new composition-corner
validation checks pass; source/hash/runtime/depth audits above pass. Separate
C++ physics probe builds and evaluates both512-zone profiles. No evolution
executable rebuilt/replaced, no original inputs or restarts modified. New bulk
outputs remain local. Continue the atmosphere interpolation diagnosis, then
identified-family evolution and the outstanding cold-matter work in
`docs/COLD_REMNANT.md`; a100K evolved remnant remains the objective.

---

# Ctrl-C diagnostic inspection — 2026-09-10 11:17 UTC

Read-only process inspection for the user's Ctrl-C question; no jobs were
signalled and no calculation code or input families were changed. The cold-remnant
work checkpoint below remains the implementation handoff.

- Codex PID46513 is in foreground terminal group46512 onttys000. Ember
  controllers occupy separate session/process groups with no controlling TTY.
  A terminal Ctrl-C is therefore not itself a process-tree-wide cancellation.
  Exact user-observed keypress behavior remains unconfirmed; clarification asked.
- Live forward star PID69232 (wrapper69230/group69230) is progressing,
  about1.10T at inspection, progressing toward3.3T. Native checkpoint updated within a second
  of inspection. This is the new fresh run; furthest COMPLETED result remains2.85T.
- Four old source attempts still cycle:64973/64986 under warm-grid controller
  64914, and66795/66841 under CONREF-comparison controller66522. They consume
  about four cores. Repeated corrections persist around80-109 iterations.
  Separate slow direct canonical control66776 under66769 also uses a core;
  it is changing, not proven to be in the same two-cycle.
- warm-grid and cool-initialization cancellation.json records only prevent
  subsequent launches. run_nongrey_plan.py explicitly lets active jobs finish;
  ThreadPoolExecutor also waits for running workers. Source/snapshot wrappers
  wait on children without an explicit interrupt cleanup/final-receipt path.
  No current orphaning or failed Ctrl-C receipt was established from these facts.
- Productive controllers69559 (XH=.1 fine warm extension) and69580
  (interpolation separation) each had two canonical models completed and
  two source children running. Old condensate waiter27809 remains asleep
  awaiting an install receipt that does not exist, using negligible CPU.
- Official CLI documentation distinguishes /ps (inspect background jobs)
  from /stop (stop ALL background terminals in the current session, including
  useful work). No /stop was executed. Read-only ps escalation approved;
  no approval-review rejection in this diagnostic turn.

---

# Cold-remnant work in progress — 2026-09-10

**Latest checkpoint: 10:33 UTC. This section supersedes all older status
paragraphs below. The furthest evolved star remains 2.85 trillion years;
no hydrogen-exhaustion or white-dwarf result has been reached.**

The user authorized autonomous local implementation and calculation through a
very cold white dwarf. The target is cooling-branch Teff=100 K with supported
physics, with intermediate exhaustion and cooling diagnostics described in
[docs/COLD_REMNANT.md](docs/COLD_REMNANT.md). No new publication is requested.
Preserve the original executable, restarts, and input families. New calculations
and source builds use separate paths. Do not interact with parent-owned
controllers. Bulk generated data remain local/ignored; leave `.DS_Store` alone.

## Completed stellar control and implementation

- **Fresh forward512-zone star is running toward3.3T**, session80141,
  `out/evolution-cold-remnant-forward-512-3300gyr-v1.{json,log,restart}`.
  Copied executable and immutable input receipt:
  `/tmp/ember-cold-remnant-forward-3300gyr-v1`; wrapper log same path+.log.
  Uses the96-node gas atmosphere extension, refined64-plane EOS, TOPSv4,
  ledoux-diffusive transport, SFII-SVH rates, WD conduction, and plasma-HRW
  thermal neutrinos. At10:32UTC it reached0.273T, fully convective.
  This is an exploratory gas track, not an exhaustion or cold-remnant result.
- Fresh 512-zone control with refined EOS/TOPS and original gas atmosphere
  **completed 0->2.85T**, 4023 CPU seconds / 6753 elapsed seconds. Output
  `out/evolution-cold-remnant-control-512-2850gyr-v1.json`, native checkpoint
  same prefix `.restart`, immutable binary and receipt in
  `/tmp/ember-cold-remnant-control-2850gyr-v1`. No physical input changed
  during the run; source edits were independent of its copied executable.
  Session16925 done. `compare_evolution_control.py` comparison report
  `docs/results/evolution_cold_remnant_control_2850gyr_v1.json`: relative
  L difference1.369e-7, R3.739e-8, central H-1.413e-7;
  Teff difference0.00004823K. This is consistency, not lifetime convergence.
- Original star `out/evolution-metal-512-2850gyr-gas-restart.json` remains
  untouched: XH=.3029, X3=.004663, R=.1488 Rsun,
  L=.001858 Lsun, Teff3107K, Tc6.879MK,
  rhoc241.2g/cm3, fully convective. Original `build/apps/ember-evolve`
  SHA d04ae649221f0d35c448e438fea620cb72089c6a983498b34c0b2fad7877f974.
  Native restarts require identical binary/physics; no migration bypass.
- Separate Release build `/tmp/ember-cold-remnant-build-v1`: **32/32 CTests
  passed** after plasma-neutrino coupling and driver cache. Full log
  `/tmp/ember-cold-remnant-neutrino-full-tests-v1.log`. Later bare-rate cache
  passed its4 targeted nuclear/mixing tests.25 atmosphere Python tests pass.
- `--thermal-neutrinos none|plasma-hrw` selects an explicit HRW1994 plasma
  sink, with analytic T/rho derivatives, central/zone energy terms and
  checkpoint identity. History adds thermal-neutrino L as column21.
  **Not selected in either2.85T run.** Fully ionized Ye approximation;
  other thermal channels omitted. [docs/NEUTRINOS.md](docs/NEUTRINOS.md).
- Exact caches reuse complete nuclear responses in the single-threaded
  driver and temperature-only SFII bare quadratures in a thread-local
  bounded map. Screening and composition remain state dependent. Two ABBA
  200G benchmarks each reproduced all output bytes (115 accepted steps),
  reducing CPU~12.5% and another10.8%, roughly22% together. Reports
  `nuclear{,_rate}_cache_benchmark_v1.json`; sessions22632/92785 done.
- `--opacity-directory` and snapshot receipts pin all selected input planes.
  Energy history distinguishes deposited nuclear, nuclear-neutrino,
  thermal-neutrino and last-halfstep gravothermal luminosity; initial/restored
  gravothermal value is null. Endpoint summaries retain crossing brackets
  and recrossings; they do not automatically label a model a white dwarf.

## Latest cold-matter implementation

- Cold analytic-electron first thermal derivatives and mixed density response
  now use paired Fermi-surface quadrature. Independent Sommerfeld tests at
  T100/1000/10000K andrho1e3..1e6 pass without loosening tolerances.
  The initial test exposed a second mixed-derivative cancellation and the
  final implementation repairs both. `cold_electron_response_v1_audit.json`.
  **32/32 full CTests passed** again after this change, log
  `/tmp/ember-cold-remnant-cold-electron-full-tests-v1.log`, session3674 done.
  Analytic-component entropy is still unimplemented; no complete remnant EOS.
- Ioffe EOSFI22 forced-phase/free-energy audit isolated two pressure-derivative
  inconsistencies. Optional source build repairs preserve upstream bytes,
  exact patches and separate binaries. The density-derivative discrepancy
  drops from6.77 to1.05e-6 ideal-ion pressure units after Richardson
  differencing; other identities agree to1.42e-8 or better.
  `docs/DENSE_EOS.md` gives the derivation, source caveats and report paths.
  **Offline only**, no source report sent upstream and no installed dense EOS.

## Candidate physical inputs

- EOS `data/eos/exhaustion_refined_v2/freeeos300_gs98_z020.dat`,64planes,
  true XH=0, baryonic isotope mapping.528 fresh heldouts: P/E/cv/cp maxima
  .05213/.1088/.1841/.1922%. First-law defect1.82e-9;2048 original
  profile queries byte-identical. **Composition coverage, not dense/cold EOS.**
- TOPS `data/opacity/hydrogen_poor_refined_v4`, with unchanged AESOPUS
  references.6 fresh X=.225/.275 heldouts at three Z: max active-source
  interpolation difference1.706%, old-profile difference1.105%.
  Source discrepancies are not lifetime error bars. All TOPS work complete.
- Masked atmosphere v2 requires all16 interpolation/derivative corners.
  Exact grid knots may use a complete incident cell; no extrapolation or
  omission of zero-weight corners. Assembler validates canonical source
  receipts, opacities and all archived inputs; old boundary values retained.
  `audit_nongrey_extension.py` needs independent canonical heldouts.
- Canonical source preparation:
  `/tmp/ember-nongrey-hydrogen-poor-warm-v1/prepared.json`;
  initializer `/tmp/ember-nongrey-conref-indexed-v1/source.json`.
  All accepted finals use the original canonical TLUSTY executable. The
  indexed RUSSEL variant is only an initializer. Old unindexed gas control
  is still slow; no new gas initialization-equivalence result yet.
- Four18-isotherm material planes (1075..10103K) at XH=.1/.2 × X3=0/.12:
  `/tmp/ember-nongrey-hydrogen-poor-v1/plane-000..003/opacity/fort.63`.
  Warm14-row subsets1822..10103K in the corresponding `-warm-v1` directory.
  XH=.3 material extension also complete in
  `/tmp/ember-nongrey-x030-material-extension-v1/plane-{000,001}/opacity`.
- `TAULAS` does not truncate a supplied mass grid. Continuations now record
  measured seed trimming and actual final optical depth. Independent depth
  checks at XH=.2,T3200,g5.15 passed max matching-state change3.058e-5
  across actual tau1000/3000/11076; g5.4 comparison passed2.699e-6.
  XH=.1,T3400,g5.15 depths1328/3235 passed1.637e-5. Reports in docs/results.
  These are local depth checks, not global spectral or physical bounds.

## Current source and validation work

All source jobs have same-name `/tmp/WORK.log` files. Inspect
`final/validated.json` and its hashes; initializer diagnostics are not canonical
acceptance. None of these controllers automatically install a table or evolve
an accepted star.

- **Complete24-node warm extension assembled** from explicit unique source
  plan `data/atmosphere/sources/nongrey_exhaustion_warm_v1_assembly.json`.
  Candidate `/tmp/ember-nongrey-exhaustion-warm-v1.dat` has96 source nodes
  (original72 plus24), with24 explicit missing colder/H-poor nodes.
  C++/EOS/reassembly audit completed: all96 knots have supported stencils,
  all retained boundary values reproduce exactly, source/EOS density
  differences+.1112..1.9910%. Six old upper-T knots have changed one-sided
  derivatives because the extension opens an upper cell. Independent
  XH=.25/X3=.005/3100K/g5.1 comparison differs+.6904%T and-.9341%Pgas;
  this is a local interpolation diagnostic, not a global accuracy bound.
  Report `docs/results/nongrey_exhaustion_warm_v1_audit.json`.
  Selected separately as
  `data/atmosphere/nongrey_gs98_z020_exhaustion_warm_v1.dat`, SHA
  4277d9144f7718077b9179a32b6a22a6d1dbded8213207f8b303d76218f087cb,
  with adjacent pinned manifest and source selection receipt. This is the
  family used by fresh forward80141; preserve its bytes throughout that run.
  `assemble_nongrey_grid.py` now accepts
  an explicit `--continuation-plan` and pins that file along with all sources.
- Partial93-node predecessor passed retained-value and EOS checks exactly;
  `nongrey_warm_partial_v3_audit.json`. Density differences+.1114..1.9910%.
  It had only90 supported source nodes and no usable independent heldout.
- `ember-nongrey-warm-extension-shallow-v1` finished16 accepted /8 rejected
  source attempts. Original rejected warm-subset and excessive-depth cases
  remain preserved. Both old deep3200/g5.15 pilots remain the chosen nodes.
- `ember-nongrey-x020-cool-retry-v1` finished4 accepted /2 rejected.
  Rejected3000/g4.9 cells were replaced successfully by same-temperature
  seeds and CONREF10 under `ember-nongrey-x020-cool-short-initialization-v1`
  and `ember-nongrey-x020-cool-short-y120-v1`; both canonical finals passed.
  Old direct canonical comparison2537 still runs; its queued follow-up was
  cancelled between launches and moved to the successful separate controller.
- All4 `ember-nongrey-x020-warm-low-gravity-v1` models completed. Depth
  comparisons passed max6.200e-6 (X3=0) and1.979e-5 (X3=.12).
  Both deeper seed-tau1000 models are chosen in the assembled extension.
- The3000/g5.15 depth control passed max4.225e-6 (actual tau448/1268).
  Low-gravity depth comparison passed9.987e-6 under
  `ember-nongrey-x020-cool-low-gravity-depth-v1`,5095 done.
  High-gravity3000K control under
  `ember-nongrey-x020-cool-high-gravity-depth-v1`,90800 also completed;
  audit10273 passed, report `nongrey_x020_cool_high_gravity_depth_v1_audit.json`.
  All seven local depth audits are recorded in the source selection receipt.
- Independent heldouts **both canonical completed**:
  XH=.25,X3=.005,3100/g5.1 under `ember-nongrey-x025-heldout-shallow-v1`:
  matching T4633,Pgas11400000; actual bottomtau749.3.
  XH=.15,X3=.005,3300/g5.1 under `ember-nongrey-x015-heldout-atmosphere-v1`.
  The latter cannot test the grid until complete XH=.1 stencils exist.
- XH=.1,3200/g4.9 full initializers two-cycle by~.000411 after CONREF
  has turned off. Both CONREF10/20 comparisons still cycle at late iterations;
  shortening CONREF alone did not fix them. The older11-cell plan46166 is
  cancelled between future launches; its first2 active jobs still finish.
  Do not infer these source processes were stopped.
- A finer numerical opacity derivative is now tested in initialization:
  `/tmp/ember-nongrey-conref-indexed-fine-v1/{prepared,source}.json`,
  `prepare_atmosphere_initializer.py --convective-refinement --indexed-russel
  --opacity-step0.0001`. Its canonical final remains the original source.
  The same numerical remedy previously resolved a different condensate
  initialization cycle; this gas case requires its own validation.
  Both jobs under `ember-nongrey-x010-fine-initialization-v1`,16376,
  **canonical accepted**: matchingT5796/5840K and
  Pgas6510000/6463000, actual bottomtau745/754 forX3=0/.12.
  Remaining nine warmXH=.1 cells running97785 with the successful fine
  initializer andCONREF10, work `ember-nongrey-x010-warm-fine-v1`,
  plan `nongrey_x010_warm_fine_plan_specification.json`. Two workers.
  Depth/EOS/independentXH=.15 audits remain before adoption.
  The builder now emits the standard continuation `source.json` receipt.
- Four3100K/g5.1 source models atXH=.2/.3 andX3=0/.12 running41877,
  work `ember-nongrey-warm-interpolation-separation-v1`, two workers.
  Plan `nongrey_warm_interpolation_separation_specification.json`. Existing
  full material opacities reused. These separate composition interpolation
  curvature from T/gravity curvature in the independentXH=.25 comparison.
- Old nonindexed .25/deep gas control29696 and old direct coarse trial38997
  may still run. Previous own SIGTERM and OS stack sampling were denied;
  no profiler sample was obtained. Do not signal parent-owned processes.

## Condensates: corrected interpretation and outstanding physics

- **Correction:** the earlier mutable-bulk adapter audit did not model the
  actual Fortran interface. Both prepared TLUSTY and SYNSPEC retain fixed
  `cbase`, initialized once, separate from mutable gas `ccomp`/`abndd`.
  New v3 audits show **both old and new adapters exactly reproduce fresh
  references on12 fixed-bulk calls**. Only the old adapter fails synthetic
  bulk feedback. Reports `condensate_memory_{before,after}_v3_audit.json`
  supersede v2's interpretation; v2 is annotated. No history dependence of
  actual atmosphere outputs was established. The retained C++ reservoir
  copy is defensive; existing binaries were not changed.
- Full equilibrium Gibbs correction is researched offline using pinned
  FastChem gas-versus-condensed chemical potentials, including gas
  rearrangement.108state survey and finer4state phase stencils completed.
  At1769K/.00406bar Delta cp~117496erg/g/K, versus the incomplete
  grain-only289000. Phase-crossing stencils and very small finite differences
  need care; Delta cp may be negative. No runtime correction installed.
- Existing condensate enthalpy gate remains. Need consistent bulk Gibbs
  baseline, volume/entropy/all responses and opacity mass normalization;
  retain explicit distinction between local equilibrium and rainout.
- Dense EOS EIP pilot remains offline with forced classicalGamma175 phase
  choice; quantum helium melting, density inversion and thermodynamic
  derivatives require validation. No dense remnant EOS installed.

Next: follow fresh forward80141; complete and independently audit theXH=.1
warm atmosphere cells and interpolation separation controls. Preserve active
stellar inputs; any later selected atmosphere is a separately identified family.
Continue
condensate thermodynamics, near-exhaustion network/screening, microscopic
transport, dense quantum EOS and truly cold atmosphere coverage. Do not count
an extrapolated cooling law as the requested evolved cold remnant.

---

# Publication policy update — code and recipes, bulk tables local

The user's latest instruction is to push the code used to generate the tables,
not approximately 2.2 GB of generated tables and source archives. Follow
[DATA_REPRODUCTION.md](docs/DATA_REPRODUCTION.md). Generator code, source patches,
small configurations/manifests, tests and reports are versioned. New bulk data,
raw histories and exact restart executables remain local and ignored. Existing
published data history is retained; no force push or remote-history rewrite.

The prior all-data commit `e095c06` was never published. It is retained only as a
local backup and must not be pushed as an ancestor or through `--all`/`--mirror`.
The replacement publication commit has the prior remote `ff5e8c7` as its parent.
Read git log/status for the final commit and remote state. There are no changes
to the reference stellar model or numerical physics in this publication revision.
The report's `recovery_manifest.json` identifies local files, not bundled archives.

---

# Session-close checkpoint — 2026-09-09, 23:09 UTC physics snapshot

Read [the next-session prompt](docs/NEXT_SESSION_PROMPT.md) and
the September 9 scientific report (retained in Git history; superseded by
[the September 10 edition](docs/reports/2026-09-10/README.md)). This checkpoint
supersedes the status of the historical entries below. The user requested a
report, handoff and GitHub checkpoint before clearing context. No further
stellar evolution or source-physics change was made during this packaging step.

- Furthest accepted stellar model is still **2.85 trillion years**, gas-only,
  original GS98 EOS/opacity inputs, 512 points, fully convective. XH=.3029,
  X3=.004663, R=.1488 Rsun, L=.001858 Lsun,
  Teff=3107 K, Tc=6.879 MK, rhoc=241.2 g/cm3. The atmosphere
  floor remains XH=.3. No hydrogen-exhaustion or lifetime claim.
- Report artifacts archive the original 0->2.5T and both restart segments,
  receipts, final native checkpoint and original macOS arm64 d04ae649... binary.
  Their hashes are in `docs/reports/2026-09-09/artifacts/manifest.json`.
  The stitched 0->2.85T figure checks matching seam states. Preserve the original
  binary: changed executable/physics bytes are incompatible with native restart.
- Low-H EOS is available separately through XH=.1 and passed its source checks.
  Low-H opacity remains **under validation**: independent midpoint discrepancy
  11.90% over all eligible source cells; .9183% in the audit's hot rectangle.
  These are interpolation checks, not physical uncertainty estimates or direct
  errors on the existing star. Current stellar inputs remain unchanged.
- TOPS refinement controller completed all FOUR raw requests: X=.15 at
  Z=.01/.03 and X=.125/.175 at Z=.02. Raw requests/results and completion receipt
  are preserved under `docs/reports/2026-09-09/artifacts/pending_tops_refinement/`.
  They are NOT imported/independently accepted. Read actual request mixtures;
  filenames x012/x018 are rounded labels. Retain the old coarse heldout audit.
- Last collector accepted70/72 condensate cells. Both remaining solar2800K,
  logg4.9 cells now converge but fail the unchanged grain-enthalpy gate:
  X3=.12 gives4.170e-7 and X3=0 gives4.114e-7 convective flux fraction
  in condensing layers, above1e-8. The fine-initializer controller exited on
  that physical audit failure; it is no longer merely a slow source solve.
  Both logg5.0 controls also failed. No production policy relaxation.
- Both artificial extra-cp controls completed. For cp+1e6/+1e7 erg/g/K,
  fractional boundary pressure changes are7.494e-7/7.376e-6. They remain
  diagnostic, not full grain thermodynamics or a rigorous missing-physics bound.
- No condensate runtime install receipt existed at inspection. The old installer
  was stopped earlier; its queued stellar controller may still wait for the
  receipt. No controllers were stopped/restarted for the report. Logs establish
  completion/failure, not current process liveness; inspect before launching work.
- New Release build/regression validation uses a separate directory so it cannot
  replace the checkpoint executable. Final results are in the report's
  `validation.json`. The source/data checkpoint is intended for origin/master;
  inspect git log/status for the actual published commit, not this prose.

Immediate next steps: independently validate/refine the new opacity family;
resolve condensate enthalpy/convection; extend atmosphere H/Teff/gravity support;
add explicit opt-in selectors; then continue with supported inputs under a
consistent initialization/restart policy. Next structural work is microscopic
transport and moving-boundary/shell convergence through H exhaustion. The report
also details remnant atmosphere/ion physics, mass-loss and magnetic/rotation
sensitivities, high-Z composition choices and conditional ultra-long-term physics.

---

# Active forward-evolution work — 2026-09-09

**Latest checkpoint: 21:41 UTC. This supersedes the older checkpoints below.**

- Manual collector89025 completed:70/72 accepted, zero conflicts. Missing
  solar XH=.7,2800K,g4.9 at both He3 values. He3=.12 is the physical grain
  rejection, He3=0 still running fine-derivative initializer6394 (iteration19
  at21:41, smooth slow corrections about2e-4; previous two-cycle removed).
- Actual-corner thermal control12314 still active; cp1e6 complete as below,
  cp1e7 at iteration8 with corrections~3e-7 at21:41. No policy relaxation.
- New independent TOPS heldout X=.15,Z=.02 ARCHIVED under
  data/opacity/sources/hydrogen_poor/heldout_x015_z020. Generic reproducible
  script scripts/audit_tops_heldout.py compares log-opacity interpolation
  between source X=.1/.2 on3496 common original cells, excluding54 substituted
  cells. Maximum discrepancy11.90% atT=.002keV,rho1.585g/cm3; maximum
  hot-rectangle error.9183% atT=.03keV,rho.1585. Report
  docs/results/opacity_hydrogen_poor_heldout.json. Therefore the separate new
  opacity family is UNDER VALIDATION, not ready merely because its runtime
  support/derivative checks pass. Old track inputs remain unchanged.
- ACTIVE64410: /tmp/ember-tops-hydrogen-poor-refinement-controller.{py,log}.
  Sequential remote requests forX=.15 atZ=.01/.03 thenX=.125/.175 atZ=.02,
  work /tmp/ember-tops-hydrogen-poor-refinement. ExistingX=.15,Z=.02 heldout
  may become a refined node ONLY with a new family/report; retain the coarse
  comparison. Quarter-point results provide independent refinement checks.
  Source filenames use rounded100X labels (012/018 for.125/.175); actual
  mixtures are explicit in requests/results and must be checked on import.
  This controller only fetches raw sources, never installs data.
- EOS/opacity data READMEs now document separate families, reproduction and
  interpolation limits. FORWARD includes new opacity heldout and current
  condensate physical issue.2.85T PNG visually checked.13 condensate tests
  PASS in /tmp/ember-plot-env/bin/python (system python lacksnumpy),23nongrey
  testsPASS, Python compile andgit diff --checkPASS. No new C++ driver changes.

**Latest checkpoint: 21:39 UTC. This supersedes the older checkpoints below.**

- Latest user asked whether work continues. Background atmosphere jobs had
  continued; this turn resumed active development and advanced the gas track.
  No goals, subagents, commits, pushes, external messages or pending approvals.
  Preserve the extensive preexisting working tree. AC caffeinate session39008
  expires about22:28UTC.
- FURTHEST COMPLETE: **2.85 trillion years**, gas atmosphere and original
  EOS/opacity selections. Controller39319 finished:37 accepted/0 rejected
  continuation steps, totals1673/1. Output
  out/evolution-metal-512-2850gyr-gas-restart.json SHA256
  b2d67c169fb29d6a87b2fb66ea8722427a2ce35c8f6b9dc475a6b2b6f15ed0b4.
  Native checkpoint out/evolution-metal-512-2850gyr-gas.restart; binary remains
  d04ae649221f0d35c448e438fea620cb72089c6a983498b34c0b2fad7877f974.
  Receipt confirms no data/restart-input changes. Reports
  docs/results/{evolution,transport,nuclear}_metal_m010_2850gyr_gas.json.
  History and PNG/PDF cover ONLY2.75->2.85T. Earlier histories are retained.
- Final XH=.3029, X3=.004663,
  X4=.6724, Z=.02; R=.1488Rsun,
  L=.001858Lsun, Teff3107K,
  Tc6879000K, rhoc241.2g/cm3,
  Pc1.652e+17dyn/cm2, logg5.093.
  Fully convective/homogeneous; minimum diffusive/adiabatic ratio3.199,
  maximum local conductive fraction.03289, reduced ppII/pp.002966.
  Nuclear rest-mass release fraction.002827; gravity remains fixed baryonic.
  XH is now close to the atmosphere floor. Substantial further evolution needs
  atmospheric support below.3. No late-time mesh/time refinements yet; do not
  infer a precise physical radius turning point from the tiny turnover.
- NEW separate EOS family installed:
  data/eos/freeeos300_gs98_hydrogen_poor_z020.dat, SHA256
  93cfb5809dd3548b003fe33aa009f1d572b0303eb58e75fd30c2a961d2cd5924.
  24 planes, XH=.1/.125/.15/.175/.2/.25/.3/.4/.5/.6/.7/.75,
  X3=0/.12. Same material grid, source options, masks and isotope convention.
  All original planes and original family unchanged. Coarse .1 spacing FAILED
  at1.22% cp/cv interpolation error; retained coarse audit and computed extra
  planes before acceptance. Final336 fresh source states PASS: max P.04985%,
  E.1087%, cp.1934%; new XH<.3 only P.03491%,E.09190%,cp.1548%.
  First-law/response errors<7.2e-10. All1536 states from2/2.5/2.75T give
  byte-identical old/new EOS probe output. Reports metal_eos_hydrogen_poor_*.
  Three old source slices (963 states) reproduce exactly with ORIGINAL probe
  /tmp/ember-freeeos-probe (f83fff2498551e7ef02724ede8ea2d594c2487340a6d70c59d456681cdf76114).
  The similarly named /tmp/ember-freeeos-reproduce/probe is different; not used.
  New scripts assemble_metal_eos_family.py and audit_metal_eos_family.py.
  Raw sources, combined manifest and three parent manifests archived under
  data/eos/sources. Source works /tmp/ember-metal-eos-hydrogen-poor{,-refinement}.
  Offline runtime probe /tmp/ember-metal-eos-hydrogen-poor-probe.
- NEW separate opacity directory data/opacity/hydrogen_poor:30 TOPS planes,
  X=.09/.1/.2/.3/.4/.5/.6/.65/.7/.75 at Z=.01/.02/.03. Original21 planes
  and old runtime axes/values retained byte for byte. Nine new verified source
  requests archived in data/opacity/sources/hydrogen_poor. Extra X=.09 is
  necessary because baryonic XH=.1 with He3=.12 maps to source X<.1.
  Original AESOPUS files referenced unchanged. New scripts
  mixture_opacity_probe.cpp and audit_opacity_extension.py. Runtime probe
  /tmp/ember-hydrogen-poor-opacity-probe. Audit PASS:1536 old profile states
  byte-identical;120 new states supported; responses agree to2.66e-9 scaled;
  XH=.08 extrapolation rejected. docs/results/opacity_hydrogen_poor_audit.json.
  This is implementation/source-support validation, not independent opacity
  accuracy. NEW heldout TOPS X=.15,Z=.02 fetch40755 completed in
  /tmp/ember-tops-hydrogen-poor-heldout: archive and compare against.1/.2 next.
  Driver wd still selects ORIGINAL opacity directory; no new selector yet.
- Condensate production rectangle remains original72 cells with g4.9/5.15/5.4.
  Automatic installer30387/PID30687 was STOPPED after a real physical audit
  failure; its65/72 inventory is stale. Queued stellar controller11101/PID27809
  still waits /tmp/ember-condensate-forward-installed.json. No condensate
  runtime table or install receipt exists. Do not restart installer until the
  physical issue is addressed; it otherwise repeatedly revalidates everything.
- Actual production corner XH=.7,X3=.12,2800K,g4.9 converges and conserves
  chemistry, but grain-bearing layers carry4.170e-7 of total flux by
  convection, above unchanged1e-8 guard. grain_enthalpy_supported=false.
  Work /tmp/ember-condensate-forward-retries/focused-plane-005-g490/plane-005/teff-2800-gravity-000.
  Archived diagnostic in sources/condensation_grain_transport_rejections/
  solar_he3120_2800_g490; report condensation_grain_transport_rejection.json.
- Tested g5.0 because it brackets the actual stellar trajectory. BOTH solar
  He3 controls still fail grain guard:6.232e-7 (He3=.12) and7.416e-7
  (He3=0). Moving the boundary does not resolve missing grain thermodynamics.
  Converged canonical outputs archived as diagnostic rejections under
  condensation_grain_transport_rejections/solar_he3{000,120}_2800_g500.
  Do not keep narrowing gravity just to seek a numerically zero convective flux.
- NEW run_condensate_atmosphere.py --checkpoint-initial uses attested,
  unconverged profiles only as starting guesses. Verifies source/input/profile
  hashes and condensate physics; canonical replay supplies acceptance.
  --convective-iterations now allows0 for indexed chemistry with CONREF off.
  Fixed completed.json requirement to apply only when checkpoint is None.
  All13 condensate and23 nongrey tests passed before that last one-line fix;
  actual fine-initializer/canonical replay exercised the fix successfully.
  Source /tmp/ember-condensate-indexed-fine-derivative-source/prepared.json
  combines indexed RUSSEL and opacity derivative step.0001 for INITIALIZATION
  ONLY. It reduced one slow g5.0 solve to4 iterations, then canonical1 iteration.
  No production source physics or grain acceptance gate changed.
- ACTIVE6394: /tmp/ember-condensate-g490-fine-controller.{py,log}, work
  /tmp/ember-condensate-forward-retries/solar-g490-fine/plane-004/{initializer,final}.
  Recovers He3=0,g4.9,2800 from preserved oscillating checkpoint with fine
  derivative initializer, then canonical replay and independent chemistry.
  Old g4.9 trials after CONREF still two-cycled; stopped only owned sources,
  preserving all outputs. Checkpoint attestation explicitly records hashes
  AFTER stopping, without claiming a missing pre-run fingerprint comparison.
- Controller41958 COMPLETE: all four solar g4.9,3000/3200,He3=0/.12 sources
  and chemistry. Work /tmp/ember-condensate-forward-retries/solar-g490-hot.
  Manual collector should now include these and final focused hot cells.
- ACTIVE12314: /tmp/ember-condensate-production-thermal-controls-controller.{py,log}.
  Diagnostic cp+1e6/1e7 erg/g/K controls on the ACTUAL failed production
  XH=.7,X3=.12,2800,g4.9 corner, same canonical21 capacity. Sources
  /tmp/ember-condensate-production-thermal-source-1e{6,7}/prepared.json.
  Work /tmp/ember-condensate-production-thermal-controls/solar_he3120_2800_g490.
  cp1e6 complete: delta P7.494e-7, T-1.555e-7, rho1.151e-6 relative.
  cp1e7 pending at21:35. Archives condensation_grain_transport_controls;
  final report condensation_production_corner_heat_capacity_controls.json.
  These are artificial sensitivity controls, NOT a full grain EOS or rigorous
  missing-physics bound. Default production guard unchanged. Need actual
  condensate formation-enthalpy estimate at flagged face if useful; existing
  estimate_condensate_enthalpy.py does this diagnostically, not full grain EOS.
- Next: finish these checks, archive heldout TOPS and quantify interpolation,
  document new data reproduction, update latest atmosphere inventory, then
  extend atmospheric composition and wire opt-in inputs. No warm-block
  atmosphere format, opacity CLI selector or grain acceptance exception has
  been implemented. Existing2.85T track includes NO condensate feedback.

**Latest checkpoint: 20:32 UTC. Active atmosphere work continues; older checkpoints are historical.**

- NEW FURTHEST COMPLETE:2.75 trillion years. Controller49137 completed
  successfully. The native continuation from2.5T used94 accepted steps and
  zero rejections; total history1636 accepted/one rejected. Output
  out/evolution-metal-512-2750gyr-gas-restart.json SHA
  018d01593338c99562bb662b6d3285402b7d20e64576796adeac62b8824d19a7.
  Receipt confirms no data or restart-input changes. Native output checkpoint
  out/evolution-metal-512-2750gyr-gas.restart. Binary remainsd04... .
- Final: XH=.3198, X3=.005907,
  X4=.6543, Z=.02; R=.1489,
  L=.001826, Teff3092K,
  Tc6720000K, rhoc240.4g/cm3,
  Pc1.645e+17dyn/cm2, logg5.093.
  Fully convective and homogeneous. Min diffusive/adiabatic ratio3.657,
  maxlocal conductive fraction.03222, plasmon/L1.412e-8,
  maxreducedppII/pp.002254, maxactiveppzeta.05597.
- Reports docs/results/{evolution,transport,nuclear}_metal_m010_2750gyr_gas.json;
  evolution PDF/PNG generated and PNG visually checked. This history/plot
  contains ONLY2.5->2.75T. Original0->2.5T remains in2500gyrgas artifacts;
  full-history He3 peak.1002 at.7107T is in that earlier report.
  Radius has a tiny turnover around2.70T, below earlier mesh sensitivity;
  do not claim a physical turning point without late-segment refinements.
  No new2.75T refinement runs launched; previous complete2T refinements remain.
- README, EVOLUTION, FORWARD, LIFETIME_SURVEY, RESTART updated2.75T and
  successful full2.5T checkpoint-build repeat. CONDUCTION and NUCLEAR now
  correctly distinguish sharedGS98 inventory from legacy carriers. DiffcheckPASS.
  Latest user list/status answer should use2.75T, not2.5T. Still no condensate
  feedback, radiative core, hydrogen exhaustion or high-Z lifetime survey.
- Condensate inventory59/72, no pending chemistry/conflicts/rejections.
  Installer30387 remains active; stellar condensate controller11101 queued.
  Solar4/g4.9 initializer is oscillating at18 iterations, but its specified
  CONREF phase lasts20 iterations before ordinary equations resume: allow
  that switch before judging convergence. If it continues oscillating, the
  validated colder seed /tmp/ember-condensate-grid-solar/plane-004/teff-2750-gravity-000
  is available for a fresh2800 continuation using --initial-model. Current
  attempt started from matching gas. Preserve failed attempts and stop only
  verified owned source processes if a replacement becomes necessary.
- AC-only caffeinate renewed as session39008 for7200s, expires~22:28UTC.
  The older75645 expires20:39UTC. No pending approval or rejection, no
  commits/pushes/subagents/goals. Beyond this checkpoint, substantial further
  gas evolution requires source composition coverage belowXH=.3; current
  XH=.3198 is close to that boundary. Do not extrapolate or weaken guards.

**Latest checkpoint: 20:28 UTC. Active work continues; older checkpoints are historical.**

- Furthest completed model remains the 2.5 trillion year gas track described
  below. The new checkpoint-enabled executable has now repeated its ENTIRE
  output byte for byte, including the full trajectory and final profile.
  Report: docs/results/evolution_metal_2500gyr_checkpoint_repeat.json.
  Controller49137 passed this gate and has started the native restart from
  2.5 to 2.75 trillion years. No completed 2.75T result yet.
- Condensate inventory is now55/72, with no conflicting or rejected models.
  Re-audited the completed3200K sources at plane0/g5.15, plane1/g5.15 and
  plane1/g5.4 using the corrected, separately labelled6000K hot-join audit.
  All three PASS: zero condensates at each join pressure, zero convective
  flux in condensing layers, grain_enthalpy_supported true. Source atmosphere
  calculations were retained. Logs /tmp/ember-hot-join-recovery-{0,1,2}.log.
  Installer30387 and condensate stellar controller11101 remain active/queued.
- Solar plane4/g4.9 old canonical source PID27653 was stopped after69
  oscillating complete iterations; original outputs and stopped.json retained.
  Replacement generator session60131, actual PID31246, now runs indexed
  CONREF20 initialization followed by canonical2800/3000/3200 solves under
  /tmp/ember-condensate-forward-retries/solar-g490-warm. Exact17T specification
  /tmp/ember-condensate-solar-g490-specification.json, log
  /tmp/ember-condensate-solar-g490-warm.log. The earlier7-Teff-spec launch
  failed before source execution; its error log was retained and corrected.
- Missing plane2/g5.15/3200 source is already running under old family
  generator PID25415; do not launch a duplicate. Plane5/g5.15 has completed
  its3200 node. Focus38800 continues plane2/3/5g4.9 and plane5g5.4.
- docs/CONDUCTION.md now distinguishes the current shared GS98 inventory
  from the retained legacy metal carrier approximation. git diff --check and
  Python compilation for the four hot-chemistry modules PASS. No further
  C++ changes or stellar-physics selections changed.

**Latest checkpoint: 20:18 UTC. Active work continues; older checkpoints are historical.**

- Latest user requests (1) physics improvements vs F77/LBA97 and (2) furthest
  .1Msun details. Delivered full table/status in commentary. Actual F77
  baseline read from GitHub f77/PLAN.md and henyey77.f, retained as
  /tmp/ember-f77-baseline-{PLAN.md,henyey77.f}. F77 already has He3,
  conduction, CNO/Be7 branching, PMS and adaptive mesh; Ember is NOT a
  strict capability superset. F77 SCVH + KingIVa/AJR83/Ferguson, CF88/GDGC,
  Hubbard-Lampe, grey atmosphere; current Ember inputs summarized FORWARD.
  No email sent. Latest heating/condensates question also answered earlier.
- FURTHEST COMPLETE:2.5T GAS, session73844 completed~19:56.15accepted,
  one rejection. XH=.3600, X3=.01015,
  X4=.6099, R=.1487, L=.001753,
  Teff3062K, Tc6375000K, rhoc239.9,
  Pc1.641e+17, logg5.093. Fully convective.
  Min diffusive/adiabatic ratio4.975; maxlocal conduction3.093%;
  plasmon/L1.350e-8; maxppzeta.05915; maxreducedppII/pp.001203.
  Integrated nuclear rest-mass release .002392 of fixed baryonic mass.
  Output out/evolution-metal-512-2500gyr-gas.json SHA
  3051b3801f1592787184b09fab58033a21a25d438306f600422daa1631a75bfb.
  Same original gas executable b904...; no condensate stellar feedback.
  Reports/plots docs/results/{evolution,transport,nuclear}_metal_m010_2500gyr_gas;
  evolution png/pdf visually checked. No2.5T numerical refinements yet.
- TWO-TRILLION-YEAR refinements40162 now ALL COMPLETE. Full repeat is byte
  identical. Mesh512->1024: L+.06163%, R-.01234%, XH-1.214e-4,
  X3-5.280e-5.4x tighter time: L-.003610%, R-.0005941%,
  XH+1.300e-5, X3+2.216e-6. Summary rerun with both comparisons;
  docs/results/evolution_metal_m010_2tyr_gas_convergence.json. These are
  numerical segment checks, not lifetime/physical error bars.
- NEW NATIVE RESTART SUPPORT IMPLEMENTED AND TESTED in apps/evolve.cpp and
  apps/evolution_checkpoint.hpp. Trailing --checkpoint FILE saves each
  accepted state; --restart FILE restores full internal log variables,
  all8abundances, nextdt, counters; --checkpoint-after N writes once while
  uninterrupted calculation continues. New output path required. Physical
  selections/tolerance/mesh/mass/executable and original data bytes checked;
  FNV1a detects accidental changes, snapshot runner adds SHA256. Additional
  unrelated .dat files do not invalidate old input identities. Output history
  on restart is explicitly only the continuation segment. No oldJSON conversion.
  docs/RESTART.md. Header currently fingerprints runtime*.dat plus selected
  external atmosphere and EOS directory*.dat; future exotic manifests with
  references outside their directory need corresponding fingerprint support.
- New binary SHA d04ae649221f0d35c448e438fea620cb72089c6a983498b34c0b2fad7877f974.
  Tests: native70Myr->1Gyr restart exact all subsequent accepted rows + final
  profile, relocated exe, final checkpoint roundtrip, rejects changed tolerance,
  changed table, truncation. CTest evolution_restart, coupled_evolution,
  secular_mixing allPASS. Native test standalone and CTest both done.
  Snapshot integration29242 also COMPLETE:10Myr->100Myr identical, SHA256
  input/output checkpoint matched, partial-history summary checked.
  docs/results/evolution_restart_validation.json. No more native changes since.
- NEW CONTROLLER49137 ACTIVE:
  /tmp/ember-checkpoint-forward-controller.{py,log}.
  First full2.5T repeat uses new copied binary, stores
  out/evolution-metal-512-2500gyr-gas-checkpoint.json and
  out/evolution-metal-512-2500gyr-gas.restart; snapshot
  /tmp/ember-forward-run-2500gyr-gas-checkpoint. Must reproduce OLD full2.5T
  JSON byte for byte, then writes evolution_metal_2500gyr_checkpoint_repeat.json.
  Only after PASS, continues saved checkpoint to2.75T with same copied exe,
  snapshot /tmp/ember-forward-run-2750gyr-gas-restart,
  out/evolution-metal-512-2750gyr-gas-restart.{json,log,receipt.json}, plus
  out/evolution-metal-512-2750gyr-gas.restart. Automatic audits/partial track
  reports named *_metal_m010_2750gyr_gas. Neither new long result complete yet.
- Condensate collector reached42/72 at20:15, no conflicts/rejections after
  hot-join fix. CURRENT INSTALLER session30387 (95138 and17785 stopped only
  waiting controllers, approved; no source jobs killed). Same recipe/log
  /tmp/ember-condensate-forward-install.{py,log}. Requires full72+EOS+2775/2900
  interpolation+qualified warm material diagnostic before installing.
  Stellar condensate controller11101 still queued, uses OLD b904binary for
  a controlled2T comparison. No condensate runtime table installed yet.
- Warm material diagnostic COMPLETE/PASS with explicit qualification:
  report condensation_warm_material_refinement.json; T+.03271%,P+.2890%.
  Fine33T model REMAINS INELIGIBLE for production (grainenthalpy false).
  Both artificial cp controls complete+archived;1e7 changes matching pressure
  4.388e-6, density6.444e-6;1e6 pressure8.388e-7. Audits check each
  archived source/input/output and same reference, require<1e-5 boundary
  differences. No physical uncertainty bound claimed. Production1e-8 grain
  guard unchanged. Thermal/controller26103,78272,76806 allfinished.
- Heldout2775 COMPLETE source+independentchemistry,grainflux0,deepestcond
  tau.0007968. Archive/heldout comparison will happen under30387 oncefullgrid.
  Old79951 overallmayexit1 dueits obsolete capacity21warm branch, unrelated.
- NEW HOT CHEMISTRY AUDIT: source already uses fully vaporized bulk gas
  above6000K. Old independent audit tried FastChem beyond its tested6000K
  ceiling for3200lowH models, then failed before output. Now actual-state
  chemistry is tested only through6000K; each hotter layer has a separate
  independent6000K query at actual pressure, required to contain no grains.
  Raw join data separately archived/labelled, never relabelled as hot chemistry.
  scripts/condensate_chemistry_profile.py, audit_condensate_atmosphere.py,
  validate_condensate_model.py, archive_condensate_model.py. Failed probe
  stdout/stderr now retained before throwing. All production guards remain.
- First hot control XH.3,X3=0,3200,g4.9 PASS:19hotlayers6020..7157K,
  P23.76..43.99bar;6000K joincond0,closure7.56e-13. Archive
  sources/condensation_validation/vaporized_x300_3200_g490 and report
  condensation_vaporized_atmosphere_join.json. Also auditedplane1g4.9 PASS.
  Tests13condensatePASS;23nongreyPASS before hothelperchange (no importer
  changes). Last git diff --check and Python compilation PASS before final
  doc updates. README/EVOLUTION/FORWARD/LIFETIME updated2.5T and2T convergence.
- Existing generator processes imported OLD validator before hothelper change:
  they may finish source+newchem then markchainfailed at3000/3200 because old
  validator expects300 actual-chem rows. Collector30387 revalidates correctly.
  If failureat3200, allnodesalreadyusable; if at3000, needcompute missing3200
  in a fresh workdir from matchinggas, canonicalsource, thennewchem. Preserve
  failures.0/1g4.9 endedat3200 andarefullyrecovered; no rerunneeded.
- Focus38800 now progressed to2/3/5g4.9;5g5.4 pending/starting. Old family
  2g5.15 at2800 initializer(7iters~20:12),3g5.15 at3200,5g5.15 at3000.
  Solar4g4.9 direct2800 is OSCILLATING (65iters,max~.06 at20:12), running
  under93335. Needs verified stop/preserved outputs and indexed CONREF20
  retry (as successful solar5.4); only that owned source process, not other
  jobs. Check actualiteration_progress.
- Caffeinate75645 expires~20:39UTC; renew if stillactive. No goals/subagents,
  commits or pushes. All narrow process approvals granted; no rejected approval.


**Latest checkpoint: 19:50 UTC. Active work continues; older checkpoints are historical.**

- Latest user asks why condensates matter when the star is warming. Answered
  in commentary: local upper layers are cooler than Teff; main purpose is a
  more accurate earlier/cooler boundary and a controlled rerun. Condensation
  generally retreats with heating. Current approximation is gas depletion
  with zero grain opacity, not full cloud microphysics.
- Gas baseline still COMPLETE at 2 trillion years. Fresh repeat under40162
  finished and gives byte-identical full JSON, SHA17e5c5ebdb7d9b8acfc2f7b17e2d467d2445c0db7ca181678298d9cc119d8d47;
  docs/results/evolution_metal_2tyr_gas_repeat.json. Mesh1024 near1.98T;
  tight512 started after repeat and is still running. No refinement result yet.
- NEW 2.5T gas run73844 started19:09:34UTC, PID28812; about2.09T at19:49.
  Recipe/log /tmp/ember-run-gas-2500gyr-controller.{py,log}; snapshot
  /tmp/ember-forward-run-2500gyr-gas uses copied2T baseline executable.
  Outputs out/evolution-metal-512-2500gyr-gas.{json,log,receipt.json}.
  It starts fresh (native restart support is NOT implemented), and runs
  transport/nuclear audits and summary/plots automatically on success.
- Warm33T model /tmp/ember-condensate-refined-2800 SOURCE CONVERGED, but
  independent chemistry flags Fconv/Ftotal=1.224e-6 in one condensing
  outer layer, above the strict production1e-8 gate. Old controller77176
  exited1 as intended. Tmatch4030K, Pmatch16460000. Numerical
  coarse/fine differences T+.03271%, P+.2890%, rho+.2432%.
  It is NOT a production atmosphere. Archived explicitly as diagnostic at
  sources/condensation_material_refinement/refined-2800-diagnostic.
- Native convection root cancellation was tested and ruled out: report
  condensation_convection_roundoff.json, relative roots agree4.44e-16.
  Diagnostic source/archive retains printing-only changes; no production
  Fortran formula changed. Actual gas cp at flagged face154500891erg/g/K.
- NEW estimate_condensate_enthalpy.py uses pinned FastChem logK temperature
  derivatives and independent fixed-P chemistry. Report
  condensation_warm_latent_enthalpy_estimate.json: at native face1769K,
  P4062dyn/cm2, added ideal-atom-referenced cp~288772erg/g/K
  (~.187% native gas cp), converged finite-difference interval scan.
  This omits atomic excitation/reference EOS differences; not a full grain
  EOS, pseudoadiabat or rigorous missing-flux bound.
- NEW prepare_atmosphere_thermal_control.py adds artificial outer cp1e6/1e7,
  taper2000..2400K, adjusts grad_ad consistently at fixed density derivatives.
  Diagnostic source only, rejected by production validator. Sources
  /tmp/ember-condensate-thermal-control-source-1e{6,7}-v2/prepared.json.
  Both source solves COMPLETE;1e7 required an approved system-sandbox retry
  after dyld failed to map libc++ (no physical iteration in failed launches).
  Session78272 complete.1e6 pressure difference8.388e-7 relative.
  Old controller76806 exited1 due the original1e7 launch failure; NEW replay
  session26103 reuses verified1e6 archive and completed1e7 source, performs
  chemistry diagnostic and archives controls. Recipe
  /tmp/ember-condensate-thermal-controls-controller.py; replay log appended
  name -replay.log. Final report condensation_warm_heat_capacity_controls.json.
- NEW diagnostic_grain_transport keyword in validate_condensate_model and
  archive_condensate_model reports unsupported grain transport without
  accepting it into production. Default production1e-8 guard unchanged.
  audit_condensate_material --thermal-controls independently verifies both
  archived thermal sources/outputs and requires boundary differences<1e-5
  to qualify the fine interpolation diagnostic. Fine model stays explicitly
  unaccepted for production. This is not a physical uncertainty bound.
  Code added; warm material execution/tests and installer update PENDING.
- Forward collector17785 reached22/72, no conflicts/rejected candidates.
  Focused chains38800,96093/90048 and old surviving chains remain active.
  Heldout2775 initializer completed; canonical source active. No condensate
  runtime installation or stellar feedback yet. Installer17785 must be
  updated only after the thermal/material diagnostic checks actually pass;
  current loaded installer still calls strict warm audit and will fail.
- Caffeinate75645 expires~20:39UTC. No new C++ changes or restart implementation.


**Latest checkpoint: 19:02 UTC. Task active; earlier checkpoints historical.**

- Focused forward continuation38800 now replaces unnecessary cold calculations:
  /tmp/ember-focus-forward-chains.{py,log}; status
  /tmp/ember-condensate-forward-retries/focused-status.json. Four workers,
  eight chains2/3g5.4 and0/1/2/3/5g4.9 plus5g5.4. First four running.
  Uses validated indexed CONREF20, minTeff2750, then canonical source.
  New generate_condensate_grid --initial-model accepts a validated colder
  seed only for a single matching composition/gravity chain. Input receipt
  records its path and checksum; all final chemistry/source checks unchanged.
- Verified/stopped six unneeded cold trials26411/27458/27575/27448/27977/28147,
  outputs and stopped.json retained; /tmp/ember-stop-unneeded-cold-chains.py.
  Stopped ONLY old retry controller27005 (session62468), preserving its
  currently running0/1g5.4 child chains. Their static status.json is now STALE;
  watch actual child running/completed/validated/chemistry files. These child
  generators will finish2750/2800/3000/3200 normally. No new queue under62468.
  Old family14019 retains2/3/5g5.15 chains only; other failures include
  deliberate preserved interruptions and previous cold-domain failures.
- Recovery74761 produced accepted2650 controls, then began costly2700
  intermediates below the production range. Verified/stopped28352/28287;
  /tmp/ember-stop-recovery-intermediates.py. New sessions96093/90048 run
  focused0/1g5.15 from those accepted2650 seeds with indexedCONREF directly
  at2750; work /tmp/ember-condensate-forward-retries/focused-plane-00{0,1}-g515,
  logs /tmp/ember-focused-plane-00{0,1}-g515.log. Recovery74761 should exit1;
  its completed controls are preserved. Allnewroots already collected17785.
- Heldout2725 now COMPLETE source and independent chemistry.2775 initializer
  has started under79951;2900 complete. Mainheldout31059 finished successfully.
- Cold material comparison script now actually PASSED:
  docs/results/condensation_cold_material_refinement.json, original17Topacity
  and electron values retained bit for bit in33T; SAME capacity33 executable
  coarse/fine T+.02652%,P+.3605%,rho+.3324%. This is a numerical
  sensitivity, comparable with the cold .8095% condensation pressure effect.
  Warm33T still active;18iterations with nonmonotonic corrections. Installer
  continues to require warm material,2775/2900 and full72(EOSchecked) first.
- Last import tests8condensate+23nongrey PASS after new changes. No C++ edits.
- Caffeinate75645 expires~20:39UTC. Stellar2T repeat/mesh still active;
  no2T refinement completion claimed. Gas baseline2T and its audits/plots done.
- LIFETIME_SURVEY/FORWARD now also record .001810 integrated nuclear rest-mass
  release relative to baryonic mass at2T; gravitating mass is still the
  conserved baryonic coordinate. Match mass convention for future cross-code
  lifetime comparisons; this is a physical approximation, not numerical error.

**Latest checkpoint: 18:52 UTC. Task still active; earlier checkpoints are historical.**

- TWO-TRILLION-YEAR GAS run COMPLETE (78884): 1339 accepted, one rejected.
  XH=.4309, X3=.02575, R=.1478, L=.001632,
  Teff3017 K, Tc5.834 MK, rhoc242.5. Still fully convective;
  minimum diffusive/adiabatic gradient ratio8.243, max conduction flux2.924%.
  Full output out/evolution-metal-512-2tyr-gas.json SHA
  17e5c5ebdb7d9b8acfc2f7b17e2d467d2445c0db7ca181678298d9cc119d8d47.
  Snapshot binary SHA b904902f4f86843dc957f3b004f1db78acfa983f93d370a083eecb581f00500c.
  UTC2365s, child user1629s/system2.144s, concurrent atmosphere jobs.
  Transport and nuclear audits COMPLETE. New summarize_forward_evolution.py
  produced docs/results/evolution_metal_m010_2tyr_gas.{json,png,pdf}; figure
  visually checked. Generic utility verifies receipts, age, same binary/data/
  physical selections, and mesh/time arguments; comparisons await refinements.
- Refinements40162 ACTIVE since18:36, two workers: mesh1024 and fresh512repeat;
  tight512 follows repeat. At18:51 repeat~.52T, mesh~.31T. Same copied binary.
  The repeat must be byte-identical; no completed2T refinement claim yet.
- Indexed CONREF control15531 COMPLETE and archived. Canonical final relative
  T-2.208e-9, P+1.209e-8, rho+1.738e-8 versus original2750 control.
  Initializer/final output NOT byte-identical; roundoff-level differences.
  NEW validated indexed initializer explicitly selected in new plan
  /tmp/ember-condensate-grid-plan-indexed.json and current grid-plan.json.
  Old plan preserved as /tmp/ember-condensate-grid-plan-original-conref.json.
  Running generators retain their previously loaded source; retries newly
  launched by62468 will read the indexed plan. Use exact old plan if resuming
  old directories. Canonical source unchanged.
- Solar center3750 COMPLETE through3200. Forward inventory still6/72 at18:51,
  installation17785 waits. It collects family-chains,grid-solar,grid-solar-central,
  family-warm-recovery-v2,forward-retries. All complete canonical+chemistry
  accepted nodes; no conflicting alternatives or source gaps filled.
- Solar g5.4 2800 direct model REJECTED despite small temperature correction:
  deepest-layer total flux residual .003044 exceeds .002. Retained source
  inputs/outputs. NEW retry88017 starts from matching gas2800 with validated
  indexed CONREF20, then canonical; follows3000/3200. Work
  /tmp/ember-condensate-forward-retries/solar-g540-warm, log
  /tmp/ember-condensate-solar-g540-warm.log. This root is already collected.
  Solar g4.9 direct2800 still running under93335.
- Old family14019 now has6cold failures (0/1g5.15;0/1/2/3g5.4).
  Retry62468 currently0/1g5.4, later2/3g5.4 will use indexed initializer.
  Recovery74761: both0/1g5.15 CONREF2650 succeeded; canonical replay active.
  Old2/3g5.15 CONREF2700 also nearlydone/done; not stalled.
- Fine33T cold2600 COMPLETE source+chemistry,10iterations. Relative against
  original capacity21 coarse: T+.02650%,P+.3606%,rho+.3325%.
  Archive condensation_material_refinement/refined-2600; report
  condensation_refined_2600.json. New material audit compares SAME capacity33
  coarse/fine and verifies original opacity/electron values retained exactly.
  Warm33T2800 still active under77176,11iterations, corrections nonmonotonic.
  Its acceptance remains prerequisite for production grid installation.
- Heldout2725 initializer COMPLETE; canonical active (one iteration at18:51).
  Heldout2775 under79951 still waits,2900 already complete. Warm-capacity21
  obsolete branch in79951 will fail guard; separate77176 is correct fine run.
- Installer17785 (/tmp/ember-condensate-forward-install.{py,log}) requires
  complete72, EOS audit,2775/2900 heldouts and warm33T material comparison;
  then installs nongrey_gs98_z020_condensate_tau100.dat and writes receipt.
  NEW stellar controller11101 (/tmp/ember-run-condensate-2tyr-controller.{py,log})
  waits for that receipt, uses copied gas baseline binary and otherwise same
 512/tol10 settings, then runs transport/nuclear audits. Condensate numerical
  refinements have not yet been queued. No condensate stellar result yet.
- Caffeinate renewed under75645 at18:39 for7200s (expires~20:39UTC), AC-only.
  Older50864 has expired. Stop/renew as appropriate while working.
- No C++ changes since prior checkpoint. Python compilation, new2T summary,
  and git diff --check PASS. README/EVOLUTION/FORWARD/LIFETIME docs now
  distinguish completed2T gas baseline from pending refinements/condensates.
  Latest Fred question already answered in commentary with draft and scope;
  final answer must be self-contained. No email sent; no high-Z lifetime claim.

**Latest checkpoint: 18:24 UTC. Task still active. Earlier checkpoints are historical.**

- Extended GAS family is COMPLETE, archived, audited and installed separately:
  `data/atmosphere/nongrey_gs98_z020_extended_tau100.dat`, 72 cells,
  SHA eec2e0553fb17e17fe16b82260a04977de6e0058d03b272992e638a1150d9532.
  Archive `sources/nongrey_gs98_z020_extended`; manifest SHA
  c6d0c7a5d2bdb5a5689f8692ad2b7e678321baa515b10b2b25bb0f2b0484aa2a.
  New EOS density comparison +.6785..1.9911%; old gas heldouts unchanged.
  Gas29152 ended with the expected 3 interrupted-initializer failures;
  replay71136 independently assembled all72 accepted canonical cells.
  Original installer89258 caught the expected manifest/table race and ended;
  rerun91772 completed installation. Log /tmp/ember-install-extended-gas-replay.log.
- Two-trillion-year GAS stellar run ACTIVE, controller78884, PID26838,
  started17:56:50 UTC. Snapshot /tmp/ember-forward-run-2tyr-gas; output
  out/evolution-metal-512-2tyr-gas.{json,log}. Around1.06T at18:22, fully convective.
  NEW controller40162 (/tmp/ember-gas-2tyr-refinements.{py,log}) waits for
  successful2T, then runs1024mesh and independent512repeat concurrently,
  followed by512tight. Same copied binary; repeat requires byte-identical
  stellar JSON and writes docs/results/evolution_metal_2tyr_gas_repeat.json.
- Condensate FORWARD production target is now explicitly 72 cells with
  Teff2750/2800/3000/3200, same3H/2He3/3g axes. Versioned
  nongrey_condensate_forward_specification.json. This brackets the star's
  initial2768K. The attempted126-cell2600..3200 rectangle remains rejected
  at cold XH=.3,g5.15/5.4 (bothHe3) because trialT crosses source1000K floor.
  No cold gaps are filled. Retain solved colder nodes as controls.
- Full oldcond family14019 continues; failures now0/1g5.15 and0/1g5.4.
  NEW retry controller62468 (/tmp/ember-condensate-forward-retries.{py,log})
  restarts failed chains at2750 with verified CONREF20 first guesses, two
  workers. Work /tmp/ember-condensate-forward-retries. It skips0/1g5.15,
  which are handled by separate recovery74761 (see below). At18:16 it had
  started0/1g5.4 retries. Other original chains retain their loaded policies.
- Direct2650 recovery45169 was stopped and preserved (PIDs26157/26163),
  recipe /tmp/ember-stop-cold-recovery.py. NEW recovery74761 uses verified
  CONREF20 for the first shifted-temperature gas seed, then canonical replay,
  at2650 through3200 for0/1g5.15. Work/log
  /tmp/ember-condensate-family-warm-recovery-v2[.log]. New generator options
  --initialize-first, --minimum-teff, --convective-iterations(default20).
- Solar warm2800 COMPLETE:28canonical iterations plus independent chemistry.
  Matched17T gas comparison P+.1246%,T-.02347%,rho+.1828%.
  Archive sources/condensation_atmosphere_x700_2800_g515 and report
  docs/results/condensation_atmosphere_warm.json. New
  compare_condensate_control.py verifies source/input/opacity physics and
  revalidates compressed gas/condensed originals. Figure now3columns,
  docs/results/condensation_atmosphere_controls.{png,pdf}, visually inspected.
- Heldout2900 COMPLETE, no condensates in its actual profile, independent
  element closure1.77e-13.27 initializer still runs under31059. Independent
  heldout opacity archived at sources/condensation_heldout_validation/opacity.
  NEW extra controls79951 (/tmp/ember-condensate-extra-controls.{py,log})
  queues2775 after2725 chemistry, using CONREF20 then canonical. Its other
  warm-refinement branch still points to capacity21 and will FAIL the new
  capacity guard; its2775 thread can continue. Correct refined models use77176.
- All33 solar material isotherms COMPLETE (original17 reused exactly).
  Initial96067 failed only at Python merge's old21-row guard, retaining all
  isotherms. Increased reader/merger allowance to the native reader's100-row
  temporary limit; source launch separately checks actual compiled MTABT.
  Rerun87076 reused all outputs and merged successfully; log appended
  /tmp/ember-condensate-opacity-refinement-controller.log. Table SHA
  4778a22aa790fb751919b8b49b23dbab3784e691f0ecd8d10a419598788312d4.
- Native TLUSTY was compiled MTABT21. NEW prepare_atmosphere_capacity.py
  builds isolated capacity33 source (only BASICS.FOR MTABT changes; same
  equations/object/library/flags). /tmp/ember-condensate-source-capacity33,
  compile13916 DONE. Correct controller77176 waits for merged33 table and
  runs cold/warm refined models after same17T controls. Both capacity-only
  controls PASSED source+chemistry: maxrelativeT/P/rho1.56e-6. Reports
  condensation_capacity_control_{2600,2800}.json. Refined work dirs
  /tmp/ember-condensate-refined-{2600,2800}; controls/archive under
  sources/condensation_material_refinement. Recipe/log
  /tmp/ember-condensate-capacity-controls.{py,log}. Refined results pending.
- Profile of owned CONREF initializer26907:1536stacks in2seconds,1341 leaf
  stacks in native RUSSEL (not FastChem); repeated element-list scans dominate
  that sampled phase. docs/results/condensation_initializer_profile.json,
  rawgzip in condensation_continuation_validation. No whole-run speed claim.
  NEW prepare_atmosphere_initializer --indexed-russel replaces linear element
  lookup with counts, preserving floating update order and duplicate entries.
  It is only an auxiliary initializer. Source
  /tmp/ember-condensate-conref-indexed-source/prepared.json built22885 DONE.
  Control15531 /tmp/ember-condensate-indexed-control.{py,log} compares2750
  initializer and canonical replay against originalCONREF2750, including
  byte comparisons and separate childCPU timing. At18:23 initializer3iters;
  NOT accepted or selected in production plans yet. Wait for validation.
- collect_condensate_grid.py inventories actual source models, chemistry,
  plane hashes, missing nodes and alternate-solution differences. Forward
  inventory /tmp/ember-condensate-forward-inventory.json; initially4/72,
  no conflicts. Roots: family-chains,grid-solar,grid-solar-central,
  family-warm-recovery-v2,forward-retries, plus known model paths in plan.
  Archive only after complete; heldout2775/2900 and33T checks must precede
  installation/stellar condensate feedback. Full72cond family still pending.
- Python import tests23PASS and condensed archive/control tests8PASS. Tests
  now include33-row binary assembly; no C++ changes. Update test_cases in
  condensation_import_validation.json from6 to8 after latest rerun32209.
- Caffeinate50864 still expires18:47UTC. Renew before then while working.

**Latest checkpoint: 17:49 UTC. Earlier checkpoints below are historical.**

- User asked how to respond to Fred Adams about converged low-mass lifetimes
  versus mass/metallicity and high-Z EOS/opacity. Answered in commentary with
  a draft and explicit limitations: current 1.3T segment is not a lifetime;
  detailed EOS/atmospheres fix Z=.02; He-rich is not high-Z. Added
  docs/LIFETIME_SURVEY.md and updated README's current checkpoint. Continue
  the active evolution/atmosphere task; no email sent, no survey yet.
- Gas70/72 accepted. Three CONREF200 starting models plateaued at corrections
  2.43e-6..3.68e-6. Verified and stopped owned PIDs24190/24246/24259; original
  inputs/output snapshots and stopped.json preserved in their v3 directories.
  NEW scripts/replay_nongrey_initializer.py ran canonical source from these
  guesses: all3 PASSED2iterations, corrections1.58e-7..3.15e-7. No relaxed
  final physics or tests. Stop recipe /tmp/ember-stop-gas-plateaus.py.
  Controller71136 (/tmp/ember-replay-gas-plateaus.{py,log}) now waits for
  final2 gascells (v3 plane0/1,3200,g5.4), then independently reimports all72.
  Those last2 initializers were still converging at12iterations, not stalled.
- Original gas29152 remains running the final2; its eventual exit will report
  the3 intentionally interrupted auxiliary jobs. New71136 preserves that
  manifest and reimports all72 actual accepted cells. Installer89258 may see
  the original manifest before the replacement table and exit; if so rerun
  /tmp/ember-install-extended-gas.py after atmosphere.dat exists. Do not alter
  final acceptance to compensate.2T gas78884 still waits for installed receipt.
- Condensate2600,XH=.3,bothHe3,g5.15 STOPPED when outer trial T fell below
  the source partition-function floor1000K. These are rejected, retained
  models (not a solver tolerance issue). The full126-cell2600..3200 rectangle
  is not established. Warmer production lower bound may be necessary; actual
  star track starts2768K. NEW warm recovery45169 runs bothfailedplanes at
  g5.15 starting2650, then2700/2750/2800/3000/3200. Work/log
  /tmp/ember-condensate-family-warm-recovery[.log]. Need inspect2650 support.
  --minimum-teff explicitly supports partial chains; eventual archive still
  requires a complete rectangle. Never extrapolate partition functions.
- Fullcond controller14019 runningplanes2/3/5g5.15 and0g4.9; initial2failed
  chains recorded.2/3cold2600 sourceconverged10iterations, audits in progress.
  Otherpendingchains stillstart2600; mayrequire warmer retries after evidence.
- CONREF2750 control fully PASSED:11initializer+3canonical. Agreement with
  direct2750 relativeT-1.90e-7,P+1.04e-6,rho+1.50e-6. Originalcompressed
  archives in condensation_continuation_validation/{direct-2750,conref-2750};
  docs/results/condensation_continuation_validation.json. Newstandalone
  archive_condensate_model.py revalidates compressed originals and initializer.
- New generator processes use CONREF only2750/2800 and now20 iterations before
  auxiliary equations resume ordinary solves (separate canonical replay still
  required). Existing first4 chains and45169 loaded earlier200iteration policy.
  Final source inputs/provenance record actual choices.20 is an initializer
  strategy, not an acceptance tolerance. Nativegasdefault remains3; v3used200.
- Warmcanonical2800 approachingconvergence (23iterations,correction7.51e-4).
  Heldout opacity X=.55,Y3=.1 COMPLETE SHA
  c5cc9d3b63da3c71fb5853a18c72716d81f01d844c65da0188c8f9c85ae3ac2e.
  NEW controller31059 runs2900 frommatchinggas, then2725 withCONREF200 guess;
  /tmp/ember-condensate-heldout-atmospheres-controller.{py,log}.
- source_failure now detects explicit fatal Fortran STOP messages, including
  cached failed outputs, even withzeroexitstatus/apparentearlierfinalprofile.
  tests/test_nongrey_import.py23PASS; condensate6PASS. NoC++change.
- Caffeinate50864 expires18:47UTC; renew whileworkcontinues.

**Latest checkpoint: 17:25 UTC. The 17:11 and 16:48 details remain below.**

- Gas grid66/72 accepted. Remaining3200 cells still run; two newest3200g5.4
  initializers started after3000g5.4 passed. All5 remaining condensate opacity
  planes COMPLETE; controller14019 is finishing their archival and launching
  first4 composition/gravity chains. Source archives use names
  `condensation_opacity_x{300,450,700}_he3{000,120}_extended` (solarx700/y0
  retains `condensation_opacity_x700_extended`).
- Solar2650 atg4.9 and5.4 now also PASSED. Warm central2750 PASSED16
  canonical iterations and independent chemistry;2800 is running, with
  large deep corrections still declining after8 iterations. Do not claim
  the2800 result or full condensate family accepted yet.
- CONREF2750 initializer PASSED11 iterations; canonical replay currently
  underway at `/tmp/ember-condensate-conref-2750/final`, controller72845.
  Validate/audit and compare against independent direct2750 once complete.
- Plan changed at17:20 to explicitly select the verified CONREF initializer
  for2650/2700/2750/2800 guesses in NEW generator jobs. Every cell still runs
  separate canonical source and full independent chemistry checks. Original
  plan saved at `/tmp/ember-condensate-grid-plan-direct.json`.
  Existing solar93335 loaded the original plan and continues direct; use
  the direct-plan backup if restarting that existing work directory.
- generate_condensate_grid now accepts --convective-initializer or explicit
  plan key. run_condensate_atmosphere records initializer source/output
  provenance. Validator checks every initializer input/output hash and
  reproduces canonical fort.8 exactly; family archive preserves originals.
- Chemical acceptance now also rejects nonfinite/nonpositive conserved
  nuclei densities and uses absolute convective heat flux in the grain
  enthalpy restriction.6 archived-model tests PASS; numerical outputs of
  accepted models unchanged (their condensing heat flux is zero).


- Independent new gas heldout COMPLETE: XH=.375,X3=.06,3000K,g5.15.
  Canonical2 iterations afterCONREF14. Interpolation errors T+.1958%,
  P+.03401%,rho+.1295%. Original opacity and atmosphere/initializer artifacts
  archived under `nongrey_extended_validation`; report
  `docs/results/nongrey_extended_heldout.json`. All4 required corners revalidated.
- NEW2T gas controller78884 waits for installed extended72 table, then uses
  `scripts/run_evolution_snapshot.py` with current build,512points,tolerance10,
  GS98 EOS,ledoux-diffusive. Work `/tmp/ember-forward-run-2tyr-gas`, output
  `out/evolution-metal-512-2tyr-gas.json`. Recipe/log
  `/tmp/ember-run-gas-2tyr-controller.{py,log}`. It runs transport/nuclear audits.
  Queued only: not launched yet. Snapshot records UTC,awake,child CPU separately.
- Condensate solar2600K at g4.9 and5.4 both source/chemistry PASSED in
  `/tmp/ember-condensate-grid-solar`;2650 next. Warm2750 canonical has recovered
  from initial deep corrections and is converging. CONREF2750 control continues.
- NEW full condensate family controller14019 waits for all5 remaining opacity
  planes, archives them and starts15 non-solar composition/gravity chains.
  Capacity4 chains before gas-grid installation,10 afterward; one source worker
  per chain. Work `/tmp/ember-condensate-family-chains`; recipe/log
  `/tmp/ember-condensate-family-controller.{py,log}`. No chains launched yet.
  It writes `/tmp/ember-condensate-grid-archive-plan.json` with opacity archives.
- NEW solar-center completion controller3750 waits for validated/audited2800,
  then solves gas-seeded3000/3200 with generate_condensate_grid (known lowerT
  cells reused). Work `/tmp/ember-condensate-grid-solar-central`; recipe/log
  `/tmp/ember-condensate-solar-completion-controller.{py,log}`.
- Retain all solved2650/2700/2750 nodes in final condensate table:126cells,
  dense spec `nongrey_condensate_dense_specification.json`. Generator now saves
  `continuation-chain-*.json` for all7T; existing solar93335 predates this, so
  gather its intermediate validated files explicitly at assembly.
- NEW independent cold condensate heldout opacity controller5516 now running
  XH=.55,X3=.1,17T; `/tmp/ember-condensate-heldout-opacity`, logs/recipe
  `/tmp/ember-condensate-heldout-opacity-controller.{py,log}`. Planned target
  Teff2725,g5.1 lies between the retained grid nodes; atmosphere NOT launched.
  Good same-composition gas initializer: `/tmp/ember-nongrey-final-heldouts/composition`
  (2800K,g5.15); another at `.../four-axis` (2900K,g5.1).
- NEW material refinement controller96067 waits for gas-grid installation,
  then computes16 NEW midpoint opacity isotherms with2workers and runs2600K
  condensed control plus independent chemistry audit.33T,19rho,30000freq;
  original17 isotherms already copied/revalidated exactly, no recomputation.
  Spec `nongrey_condensate_refined_material_specification.json`; work
  `/tmp/ember-condensate-opacity-solar-refined`, atmosphere
  `/tmp/ember-condensate-atmosphere-cold-refined`; recipe/log
  `/tmp/ember-condensate-opacity-refinement-controller.{py,log}`.
  generate_condensate_opacity now has --reuse-opacity; original receipts and
  input paths retained with reused.json provenance. Full refined plane pending.
- New validator/importer/archive/generator scripts and5 actual archived-model
  rejection tests pass;23 existing import tests also pass. Python syntax and
  diff--check pass. Build current/no work needed. No new C++ changes here.
- `plot_condensate_controls.py` produces PNG/PDF in docs/results, visually
  inspected. Helium gas-control originals now archived alongside its model.
  At fixed shallow optical depth, depletion changes gas pressure by factors;
  deep matching differences stay below1% in these two controls.
- Caffeinate50864 expires18:47UTC. Old17865 has expired. Continue work.


- GAS grid session29152: `/tmp/ember-nongrey-extended-v3`,6workers,
  log `/tmp/ember-nongrey-extended-v3.log`.64/72 cells accepted; eight remaining
  use pressure-preserving continuation and native CONREF200, then canonical
  replay. Originalv2 remains stopped/preserved. Do not resume old strategy.
- Install/archive controller89258 waits for completed72 grid, archives to
  `data/atmosphere/sources/nongrey_gs98_z020_extended`, runs GS98 EOS family
  audit, installs separate `nongrey_gs98_z020_extended_tau100.dat`.
  Recipe/log `/tmp/ember-install-extended-gas.{py,log}`; completion receipt
  `/tmp/ember-extended-gas-installed.json`. Old48 runtime remains intact.
- Heldout XH=.375,X3=.06,3000K,g5.15: opacity COMPLETE SHA
  e1e7916f6cfabaaa41acddd45dbab4c7323428a79c8dfac20b6bd2212a8103a8.
  Controller38972 runs `/tmp/ember-nongrey-heldout-atmosphere/{initializer,final}`;
  initializer9 iterations, correction.00333 at16:47. Canonical replay and
  interpolation comparison still pending. Recipe/log `/tmp/ember-heldout-continuation.{py,log}`.
- Condensed2600K solar model and SAME17T gas control both validated/archived.
  Isolated depletion effect: Pgas+.8095%,T-.1349%,rho+1.104%.
  `docs/results/condensation_atmosphere_cold.json` updated; gas control lives
  inside `condensation_atmosphere_x700_2600_g515/gas-control/`.
  Matching gas2800 also completed at `/tmp/ember-condensate-atmosphere-gas-warm-v3`;
  not archived yet. Gas17T opacity archived in `condensation_gas_opacity_x700_extended`.
- Helium-rich condensed2800K XH=.3,X3=.12 PASSED54 iterations, chemistry
  element closure3.56e-9; condensing layers have zero convective flux.
  Same16T gas comparison P+.4833%,T-.0974%. Archived in
  `condensation_atmosphere_x300_he3120_2800_g515`; compact report
  `docs/results/condensation_atmosphere_helium.json`. No stellar feedback yet.
- Warm solar temperature chain99884:2650 and2700 PASSED source/chemistry.
  Active2750 has large deep corrections; canonical2800 queued after it.
  Work `/tmp/ember-condensate-atmosphere-continuation-{T}`; recipe/log
  `/tmp/ember-condensate-warm-chain.{py,log}`.
- NEW native CONREF condensed initializer source `/tmp/ember-condensate-conref-source/prepared.json`
  built with original bounds patch and original Newton/opacity derivative.
  Controller72845 tests2750 from accepted2700, then canonical replay plus
  chemistry audit: `/tmp/ember-condensate-conref-2750/{initializer,final}`.
  Recipe/log `/tmp/ember-condensate-conref-control.{py,log}`. Not validated yet.
- Fine opacity derivative control8716/PID23541 was explicitly stopped after
  persistent unconverged corrections; stopped.json preserves identity/output.
  Previous60487 and originalwarm60078 are also stopped. Do not claim solved.
- Full condensate family: opacity controller84722 generates5 remaining17T
  composition planes,1worker each, `/tmp/ember-condensate-opacity-x{300,450,700}-y{000,120}-extended`;
  at16:46 each3/17 isotherms accepted. Solarx700/y000 uses completedv3 plane.
  Recipe/log `/tmp/ember-condensate-remaining-opacity.{py,log}`.
- Family generator93335 runs solarplane4 atg4.9/5.4,2workers, from2600 through
  50K continuations to2800 then gas-seeded3000/3200. Work/log
  `/tmp/ember-condensate-grid-solar[.log]`, plan `/tmp/ember-condensate-grid-plan.json`.
  Other plane chains not launched yet. Plan gas source is extendedv3.
- New scripts: generate_condensate_grid,archive_condensate_grid,
  validate_condensate_model. Condensate importer uses explicit calculation
  label, full input fingerprint (including gz archives), original chemistry
  conservation/profile T/P and no heat-carrying condensing layers.
  Cold archive revalidation and altered-input rejection checks PASS; report
  `docs/results/condensation_import_validation.json`.23 nongrey import tests PASS.
- Star remains at completed1.3T with mesh/timestep refinements.2T not launched.
  After gas72 installation and heldout audit, snapshot current binary/inputs
  and launch2T, measuring UTC, awake elapsed and child CPU separately.
- NEW caffeinate session50864 expires18:47UTC, AC-only; stop when finished.
  Older17865 expires16:54UTC. No display-sleep change.

Task remains in progress: extend atmospheres/composition, condensates and
radiative-core transport, then evolve beyond 1T. No commit/push this task.
Preserve the historical handoff sections below.

Implemented and installed:
- GS98 baryonic metal EOS, 12 H/He3 potential planes, XH=.3..75, X3=0/.12.
  `data/eos/freeeos300_gs98_z020.dat`; original raw source and manifest
  archived alongside. Same representative element numbers as atmospheres;
  FreeEOS omits K (baryonic mass 4.56e-6); helium isotope entropy added once.
  144 independent source checks: max P .06164%, E .1675%, cp .2594%.
- Shared GS98 ion inventory for screening/conduction. Ioffe source tables
  expanded through Zn to bracket Ni, without changing original values.
  Opacity conversion preserves aggregate metal atomic mass, but its fixed
  native GS98 atomic-weight pattern retains ~1% individual-metal number
  differences from the representative-isotope pattern. This is approximate.
- Ledoux MLT with full-EOS face buoyancy and analytic Jacobian. Schwarzschild
  remains an explicit control. Conservative implicit burn/diffusion solver
  includes Langer mixing-only semiconvection and Kippenhahn thermohaline
  mixing; convective regions mix instantaneously. `ledoux-diffusive` selects
  alpha_sc=.1, alpha_th=1, explicit uncalibrated efficiencies. No microscopic
  diffusion/remeshing yet. Synthetic transport, stiff conservation, nonlinear
  burning and analytic Jacobian tests pass; stellar radiative core not reached.
- EOS interpolation optimized by combining four node stencils before common
  Hermite evaluation. A local F/T offset is removed before differentiating
  to prevent third-derivative cancellation. New metal response error4.14e-8;
  first-law error3.75e-10. All27 CTest suites pass (13:13 UTC).

Completed new consistent-mixture 512-point1T run:
`out/evolution-metal-512-1tyr.{json,log}`; 925 accepted/1 rejected.
R=.1442, L=.001369, Teff2923,
XH=.5331, X3=.09106, Tc5154000,
rhoc257.7. Fully convective. Uses old48 non-grey gas grid,
GS98 EOS/screening/conduction/opacity conversion, Ledoux. Run predates the
EOS optimization and slow diffusion (zero coefficients, fully convective).
Exact executable hash/wall time were NOT recorded at launch: do not invent
those. Further runs must snapshot executable and provenance before launch.

Extended source work:
- `/tmp/ember-nongrey-extended-v2`: all6 opacity planes complete,16T x19rho
  x30000 frequencies, H=.3/.45/.7 and X3=0/.12. Spec in data/atmosphere/sources.
  T1076..7762.30 retains old low-T isotherms exactly. No extrapolation.
  Original 15000K material rectangle exposed negative native HeI4471 profile
  extrapolation above its ne<=1e16 tabulation; exact probe retained in
  `/tmp/ember-helium-profile`, no ad hoc line correction.
- Active generator session97824, `/tmp/ember-nongrey-extended-atmospheres-tail.log`.
  8 workers,72 cells. Old48 revalidated. First3 newH=.3 canonical cells passed; several2800K initializers now finished and in canonical replay.
  First cold initializers take tens of iterations. Do not overwrite their
  active files or executables. Old failure.json files are historical until
  each restarted cell completes; inspect running receipts/logs.
- Repaired optional CONREF initializer source `/tmp/ember-nongrey-conref-source`;
  delta patch in data/atmosphere/sources/tlusty208-conref.patch. Canonical
  final source unchanged. Bottom20 trial T(m) gradients are starting guesses
  only; final canonical replay required. x=.3,y0,2800K,g5.15 initializer
  passed14 iterations, correction3.16e-7, flux1.04e-5, then canonical replay
  validated. Plain divergent control explicitly stopped/recorded.
- Prepared canonical source receipt
  `/tmp/ember-nongrey-reproduce/prepared-opacity-derivative.json`.
  Active TLUSTY `/tmp/ember-nongrey-opactr-source/tlusty/tlusty.exe`.

Condensates:
- FastChem4 pinned ae67cbd559bc64a3233a1cee6030b8e6b50520de. Reproducible
  offline source builder `scripts/prepare_fastchem_sources.py`; fresh receipt
  `/tmp/ember-fastchem/prepared.json`, verified-probe/static library built.
- Full48 equilibrium+rainout fixed-profile audit archived under
  `data/atmosphere/sources/condensation_gs98/` (~102MB, original source outputs,
  chemistry data, license/provenance). Report docs/results/condensation_nongrey_family.json.
  Max element closure3.66e-9; condensation reaches tau.01356 at2600K,
  .001784 at2800K, absent3000/3200K, none tau>=.1. This does NOT by itself
  show negligible radiative feedback.
- Active development: `scripts/fastchem_depletion.cpp` offline GPL bridge
  for local equilibrium removal with zero grain opacity (settled-grain
  endmember), and `scripts/prepare_condensate_sources.py`. New isolated source
  build `/tmp/ember-condensate-source-v2` session from current tool context.
  Gas opacity control matches original2371K isotherm byte for byte; gas
  atmosphere2800K,g5.15 converges1iteration, matchingT/P within1.1e-9.
  Bridge independently matches300-layerchemistry2.24e-11; 65pressure
  scans at3000/4000/5000/6000K show no condensation. Exact coupled
  equilibrium-depletion atmosphere still pending.
  No grain-size/Mie cloud model or rainout transport coupling claimed.

Remaining: finish/reimport/archive72 canonical cells; independent interpolation
and atmosphere/EOS checks; complete condensate radiative control; audit and
save compact new1T result; evolve >1T (2T initial target) with snapshotted binary;
update physics docs and final tests. Keep all source approximations explicit.

---

# ember — handoff

The helium-rich non-grey grid is **installed and validated** (2026-09-08
local evening). Its 512-point 0.1 Msun evolution reaches **one trillion
years**. Read **`docs/NONGREY.md`** for source physics, repairs, reproduction,
independent numerical checks and the new stellar comparison. The established
`cond-corrected` atmosphere remains the default; select the new grid with
`nongrey:data/atmosphere/nongrey_gs98_z020_tau100.dat`.

At 1e12 yr: R=.1436 Rsun, L=.001330 Lsun, Teff=2909 K,
Tc=5.123e+6 K, rhoc=261.6 g/cm³, XH=.5339, X3=.09509 and
X4=.3510. The star remains fully convective. He3 peaks at .1026
near 736.8 Gyr. There are 906 accepted macrosteps and one rejected attempt;
stellar computation takes 612.3 wall seconds on the M4 Max. The age excludes
formation and pre-main-sequence evolution. Hydrogen exhaustion, the blueward
turn and remnant cooling have not been reached.

Against the previous 512-point corrected-COND track with identical interior
physics and timestep tolerances, final radius changes +.0201%, luminosity
-1.174% and Teff -8.90 K. The largest local luminosity and nuclear mass-energy
imbalances are 1.12e-8 and 1.51e-8. The older COND mesh/time refinements are
not convergence studies of the new non-grey track. New non-grey mesh/time
refinements are a useful next numerical step.

**Saved results:** `docs/results/evolution_nongrey_m010_1tyr.json`, its PNG
and PDF, `nongrey_boundary_sensitivity.json`, `nongrey_family_audit.json`,
`transport_nongrey_m010_1tyr.json` and `nuclear_nongrey_m010_1tyr.json`.
Full JSON/profile and progress log are retained as
`out/evolution-nongrey-512-1tyr.{json,log}`. The final transport audit finds
zero species spread, minimum diffusive/adiabatic gradient ratio 14.70,
maximum estimated local conductive flux fraction 3.07%, and opacity densities
at least a factor 30 below their ceiling. Plasmon-only luminosity is 1.49e-8
of the surface luminosity. Maximum active-layer pp zeta=.0771 is below the
.2 screening guard. Missing nuclear channels and total thermal-neutrino
losses remain documented limitations.

**Installed family:** `data/atmosphere/nongrey_gs98_z020_tau100.dat`;
original atmosphere inputs/outputs, source receipts and hashes are in
`data/atmosphere/sources/nongrey_gs98_z020/`. The 48 cells cover
XH=.45/.7, X3=0/.12, Teff=2600/2800/3000/3200 K, logg=4.9/5.15/5.4,
fixed Z=.02, alpha=1.9 and tau=100 matching. Each uses 300 depths and
20,000 transfer frequencies. Four opacity planes each use 21 temperatures,
19 densities and 30,000 frequencies. Material temperature support starts
at 1076 K. No atmosphere or interior source extrapolation is enabled.
The filename/source label `nongrey_candidate_specification.json` is retained
as part of the archived generation provenance; the complete family is installed.

Every cell passes the unchanged flux, correction, hydrostatic, independent
chemical-density and full-profile opacity-support checks. Maximum local flux
error is .0391% and maximum temperature correction is 9.45e-7. Offline
reimport reproduces the installed grid. Runtime interpolates log T/Pgas in
four dimensions with analytic derivatives and obtains density from Ember's
EOS. The complete-family atmosphere/FreeEOS density mismatch spans
-.266% to +1.266%. Independent four-axis interpolation changes matching
T/Pgas by +.327%/-.255%; these are local checks, not global uncertainty bounds.

The source patch repairs array bounds, baryonic element masses, Rayleigh
frequency indexing, molecular-density rollback, convection Jacobian entries,
and opacity derivative temperature/state restoration. Molecular chemistry
uses tolerance 1e-8 and convection temperature derivatives use DERT=1e-5.
Native analytic probes and original control outputs are archived under
`data/atmosphere/sources/nongrey_validation/`; see NONGREY.md for details.
Final TLUSTY patch SHA256:
`0b3214a4b18b6a927f2d512475f94c5895b2790e29185e82108db35724218be1`.
It reproduces the source from the pinned tarball and a fresh build compiles.

The final high-gravity 3200 K/XH=.7/X3=.12 model needs a convective starting
profile. The archived `solar-he3-hot-convective-initialization` control uses
native refinement only for its first three iterations and converges in 11.
A separate canonical solve starts from its completed target-composition
structure, converges in one iteration and reproduces T/Pgas within 8e-9
relative. No flux law or final acceptance tolerance is relaxed.
`--initial-models` now reads completed compressed models directly from the
installed archive, verifies the original receipts against decompressed
bytes, and supplies starting structures for future source rebuilds.
Attested interrupted checkpoints remain initial guesses only.

**Retained temporary work:** complete `/tmp/ember-nongrey-family-v10`;
prepared source `/tmp/ember-nongrey-reproduce/prepared-opacity-derivative.json`;
executable `/tmp/ember-nongrey-opactr-source/tlusty/tlusty.exe`.
The V10 completion controller reports `completed` in
`/tmp/ember-nongrey-completion-v10.json`. Older V9 and derivative-control
runs retain historical/stopped outputs. An obsolete waiting controller's
`/tmp/ember-nongrey-completion.json` is not the current workflow; do not
confuse its eventual timeout with a failed current calculation.

Performance: C++ is native ARM64, targets `apple-m4` and links Accelerate;
its small Henyey blocks currently use the portable pivoted kernel. Fortran
source generation is native ARM64, generic -O2 with bounds checks, up to
eight independent atmosphere jobs. The archived source replay attributes
54% of phase CPU time to chemistry; the dense 300-by-300 solve takes .0033 s.
The radiative block repeats tridiagonal factorizations, a future optimization
opportunity. No GPU or internal Fortran parallelization is implemented.
The installed stellar-boundary table is only 2.5 KB.

All 24 CTest suites pass. The final source-reader change separately passes
its CTest suite and all 22 Python cases. Static comparisons and three
helium-rich coupled steps pass. The final plot has been visually checked;
`git diff --check` passes. No commit or push was made. Preserve the broader
preexisting extended-evolution changes and historical references below.

Further physics work includes non-grey composition coverage below XH=.45,
a matched metal-bearing interior EOS, and transport/remeshing validation
before a radiative core develops. Gas-only atmospheres omit condensates,
rainout, irradiation, magnetic support and spots. The present grid and
trillion-year experiment are intermediate main-sequence capabilities.

---

Updated 2026-09-08 after extending atmosphere/opacity composition coverage,
adding electron conduction, and reaching **one trillion years** in the
0.1 Msun experiment. Read this, **`docs/EXTENDED_EVOLUTION.md`**, then
`docs/EVOLUTION.md` for the conserved abundance/energy conventions and
coupled solver, and `docs/ATMOSPHERE.md`, `docs/CONDUCTION.md`,
`docs/NUCLEAR.md` and the data READMEs. Use git status/log for actual history.

The extended model remains fully convective, with H1≈.533 and He3≈.093 at
one trillion years. He3 peaks near .102 around 720 billion years before
its destruction catches up with production. The elapsed time starts from
a specified static X=.7, He3=0 composition; it excludes formation and
pre-main-sequence evolution. This is an intermediate main-sequence model,
not the blueward turn or a helium-white-dwarf cooling track.

At **2048 points and 1e12 yr**: R=.1435 Rsun, L=.001350 Lsun,
Teff=2920 K, Tc=5.134e+6 K, rhoc=261.9 g/cm³, X=.5325,
He3=.09344. The run accepts 926 macrosteps and rejects three attempts.
Doubling 1024→2048 points changes R/L by .00683%/.0731%; fourfold tighter
time tolerances at 512 points change them by .000917%/.00190%. Versioned
records, transport/nuclear audits and a standalone figure are under
`docs/results/*1tyr*`. These numerical differences do not calibrate the
larger atmosphere and mixture uncertainties.

Eight FreeEOS source planes now cover X=.3 through .75. AESOPUS and TOPS
have Z=.01/.02/.03 families with broader H coverage. The new opacity wrapper
preserves elemental number densities and extinction under the baryonic/atomic
and He3/He4 mapping. All original source responses and offline import scripts
are versioned; every selected input reproduces byte for byte. Normal builds
require neither Fortran nor network. Metals remain a declared EOS/opacity
proxy, and isotope-dependent cross sections are approximated.

The default atmosphere applies a grey **convective composition correction**
to the original solar COND T/Pgas boundary. This supplies an evolving mixture
response while preserving the reference boundary exactly; it is not a new
non-grey helium grid. Static alternatives change L by 7–14% when the solar
non-grey anchor is removed, identifying atmosphere physics as a leading
uncertainty. Source Teff/g and EOS/opacity bounds remain strict. The historical
frozen atmosphere and its small-abundance caps survive under `early`.

Electron conduction is now implemented from the Ioffe source tables, with
weakly damped, classic and undamped controls. Mixture resistivities and the
neutral-envelope join are documented approximations. It carries a modest
fraction of local flux on this fully convective branch and changes the
static structure below relaxation accuracy. It must remain available when
a radiative core develops. A direct metal-EOS sensitivity and a plasmon-only
loss audit are also recorded; neither is a replacement for a complete matched
mixture or thermal-neutrino module.

Twenty-four CTest suites pass, including independent source queries,
composition/thermal derivatives, atmosphere sensitivities, helium-rich
coupled evolution, conservation, rollback and timestep checks. The nuclear
prescription remains SFII/SVH, with the same documented limits: approximate
screening separate from FreeEOS, no pep/hep/ppIII/CNO, and omitted pp curvature.
Do not extrapolate SFII's 10–16 MK pep fit into this 5 MK partly degenerate
core. Conserved baryonic mass still enters Newtonian gravity; atomic nuclear
rest-mass defects supply heating without abundance renormalization.

---

## 1. What this is, and why it exists

`ember` is a ground-up C++23 stellar evolution code for the **lowest-mass
stars** — objects that burn hydrogen for trillions of years, turn *blue*
rather than red when their fuel runs out, and end as helium white dwarfs that
outlive everything else in the galaxy.

It is the successor to a FORTRAN line living at
`/Users/greglaughlin/Projects/low_mass_stars` (GitHub: `oklo/Henyey`), which
reconstructed the code of **Laughlin, Bodenheimer & Adams (1997)**. That
reconstruction succeeded and is finished; it stays as a historical artifact
with its period-styled website. `ember` is where the physics moves forward.

**Motivation from the user (Greg Laughlin):** he is discussing with **Fred
Adams** a full-scale update of their **1997 Rev Mod Phys** paper. That update
will want the ultra-cold evolution of *massive* white dwarfs — a solar
remnant, possibly objects near the Chandrasekhar mass. Those are **not to be
built yet**, but the structure must accept them smoothly. Section 6 covers
how.

**The milestone that defines "working":** a 0.1 M☉ star evolved end to end —
Hayashi track, trillions of years of hydrogen burning, the blueward turn, and
down the helium-white-dwarf cooling track below 10⁻⁶ L☉ — in one run.
**GitHub push is now authorized.** On 2026-09-07 the user explicitly asked
to push the current checkpoint and continue development, superseding the
earlier requirement to wait for the end-to-end run. The scientific milestone
is unchanged. The user then explicitly confirmed creation of the private
repository **https://github.com/oklo/ember**. `master` tracks `origin/master`;
use `git log` and the configured remote for current history.

---

## 2. Build and test

```
cd /Users/greglaughlin/Projects/ember
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

Host is an **Apple M4 Max**, 10 performance + 4 efficiency cores. `clang`
accepts `-mcpu=apple-m4`. Accelerate is linked and available, while the
current small Henyey blocks use a portable pivoted kernel. CMake and Ninja
were installed via Homebrew this session.

`-ffast-math` is **deliberately not used** — it licenses the compiler to
assume no NaN/Inf, and a stellar model legitimately probes states where a
table returns one. Those must surface.

---

## 3. Current state

Twenty-two suites pass: eos, opacity, dense_opacity, nuclear, structure,
convection, atmosphere, henyey, relaxation, stellar_equilibrium, cms19,
tops, low_mass_equilibrium, helmholtz, nongrey_equilibrium,
composition_physics, burning_mixing, nuclear_prescriptions, coupled_evolution, extended_physics, convective_atmosphere and extended_evolution. The independent n=3 radiative-polytrope
benchmark and the older ionized-EOS 0.5 Msun benchmark are retained.

| Module | Implementation | State |
|---|---|---|
| Constants/composition | `constants.hpp`, `composition.{hpp,cpp}` | CODATA 2018, nominal solar units, 8 species, AAG21 metal helper |
| Analytic EOS | `eos_components.cpp`, `eos_composite.hpp`, `fermi.{hpp,cpp}` | ions, radiation, relativistic FD electrons; component Hessians |
| CMS19 EOS | `eos_cms19.{hpp,cpp}`, `src/jet2.hpp` | original pure-H/He TP density/entropy, additive-volume mixture, actual interpolation Hessians, strict fluid support; **static only** |
| Helmholtz EOS | `eos_helmholtz`, `eos_composition` | C2 FreeEOS potentials at eight X values, He3 number-density mapping, analytic thermal/composition derivatives; explicit metal/isotope approximations |
| EOS interface | `eos.hpp` | transport response derivatives, optional density bounds, virtual PT inversion, explicit internal-energy availability |
| Opacity | `opacity_table.{hpp,cpp}`, source wrappers | strict native log R/rho tables; X/Z families and elemental He-isotope number-density mapping |
| Low-T opacity | `opacity_aesopus.hpp`, `opacity_ferguson.hpp` | AESOPUS 2.1 gas log R to 6; Ferguson with grains also available |
| Hot opacity | `opacity_opal.hpp`, `opacity_tops.hpp` | OPAL log R to 1; TOPS two un-clamped rectangles, seven X planes from .3 to .75 and Z=.01/.02/.03 |
| Blending | `opacity_blend.{hpp,cpp}` | smooth ln-kappa blends, analytic derivatives, strict overlap intersection |
| Nuclear | `nuclear_pp.cpp` | pp chains, explicit He3, atomic/baryonic bases; SFII quadrature, finite-degeneracy SVH screening; legacy control; analytic derivatives |
| Structure | `structure.{hpp,cpp}` | four residuals and analytic Jacobian; independent numerical checks; unsupported positive-dt energy rejected |
| Convection | `convection.{hpp,cpp}` | BV58 MLT, Schwarzschild criterion, bounded cubic and analytic response |
| Atmosphere | `atmosphere_grey.cpp` | radiative Eddington T(tau), variable opacity, adaptive integration and analytic sensitivities; EOS/opacity bounds intersect |
| Atmosphere table | `atmosphere_table.{hpp,cpp}` | AMES-COND tau=100 physical states; untouched MESA cells; strict bounds and explicit solar-mixture proxy; synthetic tests retained |
| Convective atmosphere | `atmosphere_composition.{hpp,cpp}` | Henyey element cooling, analytic column sensitivities, differential composition correction anchored to solar COND |
| Electron conduction | `conduction_table.{hpp,cpp}` | Ioffe source tables, three prescriptions, approximate mixture collision sum, explicit neutral-envelope join |
| Boundaries | `boundary.{hpp,cpp}` | regular unresolved central sphere and interchangeable atmosphere surface, analytic derivatives |
| Henyey/relaxation | `henyey.cpp`, `relaxation.cpp` | pivoted blocks, iterative refinement, damped Newton with residual AND undamped correction checks |
| Seeds/apps | `examples/stellar_seed.hpp`, `apps/equilibrium.cpp` | n=3 or n=1.5 Lane–Emden seed; rebuild mass after envelope adjustment; structured diagnostics |
| Evolution | `evolution.{hpp,cpp}`, `apps/evolve.cpp` | coupled backward-Euler burning/mixing/thermal structure; adaptive step-doubling on fixed mass mesh; trillion-year main-sequence experiment |
| Future cooling inputs | `losses.hpp` | full thermal neutrinos and dense-ion/remnant physics remain pending |

Data travel with the code. Source hashes, strict coverage and reproduction
commands are in `data/eos/README.md`, `data/atmosphere/README.md` and
`data/opacity/README.md`. The FreeEOS source builder/generator/importer and
COND importer reproduce the new data; all 180,865 direct source evaluations
were reproduced byte for byte from a clean FreeEOS build. Normal builds
need no Fortran or network. The earlier six CMS19/opacity imports were reproduced byte for byte with stdlib Python scripts:

- `scripts/import_opal.py`: 9,880 original GS98 Z=.02 cells, log T=4..7.1,
  log R=-8..1, 10 X planes. The original source is ragged outside this box.
- `scripts/import_aesopus.py`: 149,810 original AESOPUS 2.1 **gas** cells,
  log T=2..4.5, log R=-8..6, 10 X planes. Archive filenames retain 2.0 but
  headers identify 2.1; selected tables do not include grains.
- `scripts/import_cms19.py`: 121×441 original (log rho, log S) pairs per
  pure component on the source TP grid. Runtime masks reject unphysical
  corners. No source energy values are imported or repaired.
- `scripts/import_tops.py`: 3150 low-rectangle and 2556 high-rectangle cells
  (overlapping), native rho grid, only X=.7, Z=.02. Every server-substituted
  density is excluded. Exact returned text and request are versioned under
  `data/opacity/sources/`; a fresh service request is not required to build.

The generic `TabulatedOpacity::Range` now names its density fields
`logD_min/max` and includes the `DensityAxis`. log R requires the -3 term
in the temperature derivative; native log rho does not. All strict source
bounds remain enforced. TOPS itself blends rectangles over log T=5.6..5.7;
the stellar driver blends AESOPUS to its selected hot opacity over 4.4..4.5.

---

## 4. Reproduce and continue

### Current coupled evolution reference

```sh
build/apps/ember-evolve 2048 1e12 1e8 10 sfii-svh wd > out/evolution-extended-2048-1tyr.json
build/apps/ember-evolve 512 1e12 1e8 2.5 sfii-svh wd > out/evolution-extended-512-1tyr-tight.json
python3 scripts/verify_composition_data.py --extended
```

Positions are points, duration years, initial step years, tolerance multiplier,
nuclear model, transport and atmosphere. Default transport is `wd`; default
atmosphere is `cond-corrected`. `--help` lists the sensitivity choices.
Use `sfii-svh early` for the historical 20-Gyr physics or `legacy early`
for the 10-Gyr checkpoint. The latter still reproduces the original 1024-point
radius/luminosity within 2e-11 relative after the EOS performance changes.

JSON stdout includes accepted history and the final profile; stderr reports
progress. The longer references use local timestep tolerances ten times the
short-run defaults, with a separate fourfold tightening comparison. Each
macrostep compares one full versus two half steps and retains the two-half-step
solution. Failed steps roll back. Convection uses connected Schwarzschild
regions and instantaneous mass-conserving mixing, iterated with thermal
structure. New history columns include convective mass, Tc/rhoc and surface
abundances. There is no checkpoint/restart interface or adaptive mass mesh yet.

Source hashes, grid coverage and all reproduction/audit commands are in
EXTENDED_EVOLUTION.md and the data READMEs. The original external FreeEOS
probe is `/tmp/ember-freeeos-reproduce/probe` on the development host;
rebuild instructions are in FREEEOS.md. Independent reference generation:

```sh
python3 scripts/generate_freeeos_composition_reference.py /tmp/ember-freeeos-reproduce/probe /tmp/extended-reference.dat --extended
```

### Historical FreeEOS/COND static reference at `6cc4fb4`

```
build/apps/ember-equilibrium 4096 .1 .15 --eos freeeos --hot-opacity tops --seed-index 1.5 --atmosphere cond-solar-proxy > out/equilibrium-m010-nongrey-reference.json
```

The following values belong to `6cc4fb4` before screening correction;
the same command now gives R=.1294, L=.0009404, Teff=2810 K.
The old checkpoint had nine accepted updates, residual 6.60e-13, correction 6.76e-12, nuclear
balance about 2e-14 relative and independent virial error -6.396e-7.
The 1024/2048/4096 models agree to <.014% in R and <.078% in L between
successive meshes; virial error falls by four each doubling, while R/L
changes are nonmonotonic. Versioned metadata and the separate EOS/atmosphere
comparison are in `docs/results/equilibrium_m010_nongrey.json` and
`docs/results/eos_atmosphere_comparison.json`.

FreeEOS source options are EOS1 `(3,1,-2)`, H=.7/He=.3, all metals zero.
The C2 potential grid is .0125 dex in log T/log Q with strict stability
masks. Read `docs/FREEEOS.md` before modifying it: the H2 partition-function
join at 9000 K means local source identities alone are insufficient.
Energy consistency is exact for the implemented potential; the source
response agreement has a separate, sampled 1–2% precision qualification.
A fresh external FreeEOS process per isotherm avoids a source cache failure
at dense-to-dilute row resets; never accept nonzero source `info`.

COND uses only untouched Teff=1800..3300 K, log g=3.5..6 MESA cells, matched
at Rosseland tau=100. `--tau-top` is grey-only. The seed must use the
atmosphere's local T, not Teff. Do not relabel its native GN93 mixture as an
exact match to the interior, or silently extrapolate outside this rectangle.

### Retained CMS19/grey reference

```
mkdir -p out
build/apps/ember-equilibrium 4096 .1 .15 --eos cms19 --hot-opacity tops --seed-index 1.5 > out/equilibrium-m010-reference.json
```

Use 2048 points for a smaller working model. Positional arguments are
points, mass/Msun, seed radius/Rsun. Radius is free to relax. `--tau-top`
can be set explicitly; default .001 for CMS19, 1e-6 for the ionized EOS.
All composition is fixed at X=.7, Y=.28, Z=.02, He3=0. CMS19 explicitly
represents metals as helium (effective Y=.3). MLT alpha=1.9 is uncalibrated.
The surface is a radiative grey atmosphere matched at tau=2/3. No grain
opacity, conduction, composition mixing or age advancement is included.

At 4096 points: 12 updates, residual 4.103e-10, undamped correction
3.063e-10, R=.1286 Rsun, L=.0009756 Lsun, Teff=2845 K,
Tc=4.566e+6 K, rhoc=366.9 g/cm³. Nuclear luminosity agrees with
surface L to ~2e-15 relative; independent virial error is -6.585e-7.
The unresolved center is 1.23e-10 of the mass. All profiles are ordered.

| Points | R/Rsun | L/Lsun | Teff (K) | Absolute virial error |
|---:|---:|---:|---:|---:|
| 128 | .1291 | .0009629 | 2830 | 6.856e-4 |
| 256 | .1294 | .0009420 | 2812 | 1.699e-4 |
| 512 | .1289 | .0009637 | 2833 | 4.229e-5 |
| 1024 | .1286 | .0009747 | 2844 | 1.055e-5 |
| 2048 | .1286 | .0009751 | 2844 | 2.635e-6 |
| 4096 | .1286 | .0009756 | 2845 | 6.585e-7 |

Coarse thermal results are nonmonotonic. Fine meshes agree to ~.1% in L
and .02% in R, but R/L do not yet show clean second-order convergence;
virial error does. Do not call a tiny nonlinear residual spatial accuracy.
Doubling the atmosphere starting column to .002 changes the independent
256-point R/L by .000402%/.00259%; a same-mass-mesh 2048-point perturbation
also passes in `test_low_mass_equilibrium`. This does not validate grey
atmospheric physics. Full current JSON/log is in ignored `out/`; compact
reference metadata are versioned in `docs/results/equilibrium_m010.json`.

### Retained CMS19 source audit — its energy guard still applies

`docs/CMS19.md` gives the implementation and independent source audit.
CMS19 response derivatives follow the **actual interpolated entropy**.
The Maxwell identity D=P*delta/(rho*T*cp*grad_ad)-1 is not artificially
forced to zero. Max |D|=.2418 in the 4096-point star at T=27304 K,
rho=.1239 g/cm³, m/M=1.000; mass RMS=.004349. Small outer mass is
not evidence of small influence on R or L.

`python3 scripts/audit_cms19_energy.py /tmp/ember-eos2019.tar.gz` reproduces
both defects directly from original source columns:

- Pure H at log T=4.45, log P[GPa]=1.35 has D=-.2650, before ember
  interpolation. This is a source-version consistency issue as well as an
  interpolation accuracy question.
- Pure He at rho=1 g/cm³ drops in internal energy from 8.886e+13 to
  8.227e+13 erg/g between T=891251 and 1e6 K. Secant dU/dT=-6.059e+7
  erg/g/K while listed entropy cv=+9.725e7.

Do not fabricate internal energy, silently repair table cells, loosen the
Maxwell check, or replace an entropy-derived adiabat just to make an
identity pass. Exploratory U-TS/free-energy and pressure-integrated-potential
reconstructions developed negative heat capacities at source joins and
were rejected. They exist only as ignored scratch, not library physics.
`Cms19Eos` has no caloric state or He3 support; energy and unavailable
abundances are NaN and dt>0 zones/center throw explicitly.

The runtime conservatively masks whole 4×4 stencils by source density,
fluid/quantum limits, and log T=3.2..7.3. It is a selected computational
subset, not proof of accuracy everywhere inside. Grey tau_top=1e-6 falls
below its density support; defaulting CMS19 to .001 is explicit and tested.
Never extrapolate a surface integration below the EOS density floor.

### Next actions, in order

1. Obtain a non-grey atmosphere family with H/He composition coverage and a
   consistent interior mixture. The differential convective correction is an
   explicit intermediate approximation; the static sensitivity is not an error
   calibration. Keep COND Teff/g and all source support guards intact.
2. Validate later nuclear evolution: screening from a consistent free energy,
   cold/degenerate pep capture, pp curvature and other channels as needed.
   Audit ion quantum parameters and do not extend formulas outside their
   stated temperature domain just because a run continues to converge.
3. Add checkpoint/restart and adaptive mass meshes for economical longer
   tracks. Validate Ledoux/diffusive mixing and moving convective boundaries
   before following a radiative core or shell. Repeat age/composition-aware
   comparisons with independent models and observations.
4. Broaden composition/thermal sources again when actual mapped source bounds
   approach; the present EOS/TOPS minimum source X=.3 and COND maximum
   Teff=3300 K still preclude the complete lifetime. Improve physical
   metal/isotope EOS and opacity and the FreeEOS source-fit joins.
5. Add full thermal neutrinos, diffusion and dense-ion/cooling physics for
   the end-to-end 0.1 Msun milestone. Conduction is implemented but its mixture
   and partial-degeneracy uncertainties need continued checking. Massive WD
   physics remains future work. GitHub pushes remain authorized in §1.

### Retained old runs and practical source notes

```
build/apps/ember-polytrope 256 > out/polytrope-256.json
build/apps/ember-equilibrium 512 .5 .6 > out/equilibrium-m050.json
build/apps/ember-equilibrium 128 .1 .2 > out/old-ionized-m010-failure.json
```

The first two converge. The last uses the old ionized EOS/OPAL/n=3 defaults
and still fails at the OPAL density edge; it is retained as an explicit
failure diagnostic. The .5 Msun case has R=.9422 Rsun, L=.01988
Lsun, Teff=2233 K; its inflated cool envelope is not realistic. The n=3
stellar seed also proved poor for the new CMS19 .1 Msun model; n=1.5
converges. This changes only initialization, never the solved equations.

Original EOS archive URL:
`https://perso.ens-lyon.fr/gilles.chabrier/DirEOS/DirEOS2019.tar.gz`.
The source README and 2021 README recommend the 2019 IVL tables for stars,
2021 interactions for brown dwarfs; never use 2021 effective H as pure H.
Temporary originals are `/tmp/ember-eos2019.tar.gz`,
`/tmp/ember-eos2021.tar.gz`, `/tmp/ember-aesopus21-gs98.zip`, and
`/tmp/ember-GS98hz.gz`. Hashes and commands are documented with the data.
They are not needed for normal compilation or tests.

TOPS: `/results` can return a stale prepared calculation despite new
parameters. Submit through `/submit`, then verify returned composition and
grid dimensions before import. The final archived request uses the public
form identifiers, a blank mixture name, and the exact 21-element GS98 mass
mixture. Its returned warnings identify clamped densities; every such
pair is excluded. Request-specific curl `-k` was needed for the local
certificate-chain failure. All data needed for TOPS reimport are versioned.

---

## 5. Hard-won lessons — these cost real days in the FORTRAN line

Do not rediscover these.

1. **Fixed-format input is a trap.** `NPRIN=100000` is six digits in an `I5`
   field; it shifted every later field and silently ate the Hayashi-start
   flag. A whole session was spent misdiagnosing the resulting physics as an
   "MLT regression". *`ember` must never column-count input.*
2. **Fixed-width output is the same trap.** A model counter in `I4` printed
   `****` past 9999 and corrupted downstream analysis — I misreported where
   runs ended because of it. *Structured output only.*
3. **An equation form can be valid only in a limit.** Scaling the radiative
   equation by the convective efficiency φ = ∇/∇_rad is sound *only where φ is
   order unity*. Deep in a giant envelope φ ≈ 1e-6 — perfectly efficient
   convection — and multiplying the matrix row by a millionth decoupled those
   zones' temperature from the luminosity. **This was the cool-giant wall that
   stopped every giant for weeks.** Fix: superadiabatic zones use the plain
   gradient form, well conditioned at any efficiency. *Where a formulation has
   a domain, assert it.*
4. **Constants tuned for the Sun break elsewhere.** A fitting point at a fixed
   mass fraction drifts *above the photosphere* when a star swells to 300 R☉
   (measured: T(N)/Teff = 0.86, i.e. cooler than the effective temperature).
   Anything that is "a good value for the Sun" is suspect.
5. **Silent table extrapolation is lethal.** Below its floor the 1983 opacity
   table returned κ ~ 1e-9; the envelope went transparent and cool giants lost
   their Hayashi limit entirely, expanding to 314 R☉ and breaking. *In `ember`,
   `FergusonOpacity` throws outside its table. Keep it that way.*
6. **Write tests as physics statements, not numbers.** Every test caught
   something real this session:
   - electron normalisation short by √π (caught by `n = 2e^η/λ³`)
   - NaN under strong degeneracy: at β ~ 1e-4 the FD derivative integrand is a
     shell of width ~β that a single quadrature panel never samples
   - a wrong ppII branch ratio (caught by the mass defect, **not** by
     "mass fractions sum to zero" in the atomic-mass convention, where
     Σ dX/dt = −(ε+ε_ν)/c². In the new baryonic convention their sum is
     zero and atomic/integer mass ratios give the released energy)
   - **and twice the *test* was wrong, not the code**: a point labelled "the
     ideal limit" with 13% of its pressure in radiation, and a shell built by
     one-point integration checked against a centred-difference equation.
     Suspect the test too.

Two more lessons from connecting the atmosphere in ember:

- **Differentiate the actual interpolant.** Interpolating opacity derivatives
  separately missed the dependence of the temperature slope limiter on the
  input opacities. The full Ferguson atmosphere test saw a ~0.6% mismatch in
  its Teff pressure response; a direct opacity regression also failed. Both
  pass after propagating derivatives through the limiter. No opacity values
  were intentionally changed.
- **An inversion must report failure.** A small pressure residual can conceal
  a bad density when radiation dominates, and an exhausted iteration must not
  return a density as a success. `rho_from_PT` now checks both and throws.

Two more from assembling the analytic zone Jacobian:

- **Accurate values do not guarantee accurate second derivatives.** The
  original 64-node Fermi shell produced plausible pressures but roughly 0.2%
  errors in degenerate response derivatives, unchanged when the finite
  difference step shrank. Splitting the shell at its center resolved them.
  Use centered number-weighted kernels for the density-constrained electron
  Hessians, and preserve `f(1-f)` even when `f` rounds to one. The derivation
  is in `docs/JACOBIAN.md`; do not replace it with large cancelling terms.
- **Differentiate the implemented heating rate.** Approximate reaction
  energy weights and omitted screening derivatives did not differentiate the
  pp mass-defect heating value. Its derivatives now use the same masses,
  neutrino losses, and screening cap as its value. The cap is still only the
  existing weak-screening approximation, not new dense-matter physics.

From the first complete boundary-value solve:

- **Newton convergence is not spatial accuracy.** The coarse polytrope can
  converge to a tiny residual while its radius differs by percent levels
  from the continuum solution. Keep the independent Lane–Emden mesh-refinement
  test and the check that the unresolved central sphere shrinks correctly.
- **Damping is not convergence.** A small accepted step may only reflect a
  line-search restriction. Require a small undamped remaining correction
  as well as a small residual. Physics-domain errors may reject a trial;
  they never authorize table extrapolation or an implicit fallback.

From connecting physical opacity to a stellar trial:

- **Changing a seed's density changes its enclosed mass.** Keeping the old
  Lane–Emden mass mesh after adjusting envelope density produced folded
  surface layers. Rebuild the mesh consistently and inspect monotonicity.
- **Stored points do not establish atmospheric coverage.** The earlier
  Ferguson trial had every stored point below log R=1, yet its atmosphere
  reached the edge. AESOPUS has now supplied that missing cool density data;
  OPAL's hot density limit remains when that source is selected. The new
  TOPS import supplies the .1 Msun hot dense coverage. Report the failing module explicitly.
- **A checked correction can need refinement.** Finer stellar meshes exposed
  a linear backward-error failure absent in the small benchmarks. Correcting
  the original linear residual fixes this without weakening acceptance.
  Independent virial integration then verifies spatial hydrostatic accuracy.

---

## 6. Hooks for massive white dwarfs (the Adams & Laughlin update)

Deliberately not built. Nothing about them should require rearranging what
exists. **The electrons are already exact at Chandrasekhar-mass densities** —
what such a star additionally needs is *ion* physics, which is exactly why the
EOS is a sum of `EosComponent`s:

| Need | Where it goes |
|---|---|
| Coulomb energy of the ion lattice | new `EosComponent` |
| Crystallisation, latent heat, Debye solid | new `EosComponent` |
| C/O phase separation on freezing | new `EosComponent` |
| Plasmon/pair neutrino losses | `NeutrinoLosses` (declared, null) |
| Electron conduction | `Conduction` + `CombinedOpacity` (declared) |
| C/O interiors | `Composition` already carries C12, O16 |

Noted but **not** designed: general relativity in the structure equations,
which the last percent of mass before the Chandrasekhar limit would want.

---

## 7. State of the FORTRAN line (context, not work)

`/Users/greglaughlin/Projects/low_mass_stars`, branch `main`, pushed through
`6689b6d`. This session fixed three walls there:

- `cd26677` — the φ-scaling wall (see §5 item 3)
- `1807d06` — adaptive fitting point (§5 item 4)
- `6689b6d` — Ferguson 2005 opacities replacing AJR83 (§5 item 5)

Solar calibration under Ferguson: **X = 0.737, α = 1.93 → L = 1.004,
R = 1.000, Teff = 5791**. (α rose because Ferguson is ~3× more opaque near
10,000 K.) Note these tracks are **no longer LBA97's** — `IOPC=2` preserves
the 1997 physics exactly for reproducing the paper.

Results: 0.30 and 0.35 M☉ complete their lives and **do not** ignite helium
(peak Tc 4.5e7 and 6.0e7). No star has reached 1e8 K, so on these models the
**flash mass is above 0.55 M☉** — still unmeasured.

**Known regression:** with Ferguson opacities the 0.10 M☉ flagship no longer
reaches the 1e-6 L☉ stop; it stalls at L = 5.2e-6, Teff = 1704 K. Open
question worth the user's judgement: *should* Ferguson's grain opacities apply
to a high-gravity helium-white-dwarf atmosphere at all?

There may still be a background 0.55 M☉ run (`fg055`) in the session
scratchpad; it was holding its Hayashi line correctly at Teff ≈ 3400 K.

---

## 8. Practical gotchas in this environment

- **The shell's cwd resets** between tool calls. Use absolute paths or `cd`
  inside every command.
- **`pkill -f "henyey77 <"` matches nothing** — the `<` is shell syntax, not
  part of the command line. Two processes once wrote the same output file and
  produced an interleaved, backwards-running track. Use `pkill henyey77`.
- The session scratchpad is **cleared between sessions**; anything worth
  keeping goes in the repo.
- Long runs: launch with `run_in_background`, watch with a `Monitor` whose
  filter matches **every** terminal state, not just success.


Latest active checkpoint (14:18 UTC):
- Source build `/tmp/ember-condensate-source-v2/prepared.json` complete;
  new scripts fastchem_depletion.cpp, prepare_condensate_sources.py,
  generate_condensate_opacity.py, run_condensate_atmosphere.py.
  Explicit hybrid settled-grain endmember, not a full cloud model.
- Condensate opacity solar x=.7,y0 active session38437, log
  `/tmp/ember-condensate-opacity-x700-resumed.log`, 3workers, first10/16
  isotherms complete. Work `/tmp/ember-condensate-opacity-x700-v1`.
  First wrapper mistakenly expectedopacity.bin instead ofnativefort.63;
  4completed sourceoutputs recovered only after fulltable/source validation,
  with explicit recovery.json and unknownwalltime, no numericalchanges.
  New jobs usecorrectoutputreceipt. Original firstfailurelogretained.
- Gas atmo control session64711 done/passed, files
  `/tmp/ember-condensate-atmosphere-gas-v1`, useasreferenceforcoupledcond.
  Solar2800K,g5.15 usesoriginal16Ttable andsamegasEOSsource. Eqtable pending.
- Independent chemistry and gas controls archived43artifacts in
  `data/atmosphere/sources/condensation_validation/`; report
  `docs/results/condensation_source_validation.json`.
- Timed512point1.3T run active session84690, output
  `out/evolution-metal-512-1300gyr.{json,log}`. Snapshot executableandinput
  hashes `/tmp/ember-forward-run-1300gyr/receipt.json`, which recordswalltime
  oncompletion. Usesold48gasgrid,metalEOS,ledoux-diffusive. ~570Gyr last.
- New1Tcompactreport,plot(PNG/PDF),transportandnuclearaudits saved under
  docs/results/*metal*m010_1tyr*. Plotvisuallychecked. ExactGS98selectedin
  probes (metadata is descriptive stringbeginningGS98, notliterallygs98).
  Full48atmoEOSaudit newdensitydifference+.738..1.991%; report
  docs/results/nongrey_gs98_eos_family_audit.json.
- Newdocs/FORWARD_EVOLUTION.md describescurrentworkandlimits;README,
  docs/EVOLUTION.md,data/eos/README.md linkit. AddedmetalEOSLedouxfinite
  difference regression passes. All27suites previouslypassed, latesttwo
  impactedtestspass. Do notclaimnewextendedgrid orcondensateradiative
  feedback installedyet. Needfinishsourcegrid,extendstar>1T,sourceaudits.

14:23 UTC continuation: solar depleted-opacity plane COMPLETE, hash
0927d66f247971b1d4b0c45f6a549cfc07eb765b5e07336d3f892cb4a6b8e04d.
Archived original isotherms/table/provenance in
`data/atmosphere/sources/condensation_opacity_x700/` (completed archive job88816).
Absorption Rosseland diagnostics docs/results/condensation_opacity_comparison.json;
hot isotherms>=3000K identical, coldest grid opacity reductions up to99%.
Coupled clear-limit atmo2800K,g5.15 active session25191, work/log
`/tmp/ember-condensate-atmosphere-equilibrium-v1[.log]`;2600K samegravity
active session51617, `/tmp/ember-condensate-atmosphere-cold-equilibrium-v1[.log]`.
Both are only a fewiterationsin, outertemperaturecorrectionslarge; no acceptance
orboundarycomparison yet. Oncecomplete verifycondensation remains in radiative
layers (grainenthalpyomitted), archiveoriginals andcompare matchingboundary.
Timed1.3T run84690 near950Gyr, fullyconvective. Extendedgrid97824 now6newcells
validated; warmerinitializers finishing, cold2600highgravity correctionsstill
steadilydecreasing after~40iterations. Sourcegridnotinstalledyet.

14:29 UTC: 512point1.3T run84690 COMPLETE,1418wallsec underconcurrent
sourcework. FinalR.1462,L.001465,Teff2953,
X.5057,He3.071,Tc5306000,rhoc248.3.
Compactreport/completeinputhashes+transport/nuclearaudits in
`docs/results/*metal_m010_1300gyr.json`; fullout andoriginalreceipt preserved.
Mintransportgrad/ad12.99, maxconductiveflux2.837%,ppzeta.07284.
Two reference refinements active controller39959:1024points tolerances10,
512points tolerances2.5, both1.3T with SAME snapshotted executable asbaseline;
`out/evolution-metal-{1024,512}-1300gyr[-tight].{json,log,receipt.json}`.
Independent helium-rich condensed opacity plane now active3056,2workers,
`/tmp/ember-condensate-opacity-x300-y120-v1[.log]`, H.3,He3.12.
Sameverifiedsource-v2; oncecomplete run2800K/g5.15 atmosphere fromfinished
extendedgasgridplane001/model001001, andcomparelate-mixtureboundary.
Solar2800Kcondcontrol25191 stilloscillates outercorrections(.05,.576,.0464
iterations5..7),2600Kcontrol51617 decreasing(.0392,.00704,.00107).
Do notacceptorstopjustbecauseoflargeinitialcorrections. Addedpostsolver
`scripts/audit_condensate_atmosphere.py` (usevenvpython/numpy) checksindependent
FCelementpressureclosure andrejectscondensates inheat-carryingconvectivelayers
becausegrainenthalpyomitted. Needrunoncompletedcontrols,thenarchive.

## Active continuation — 2026-09-09 15:36 UTC

The 512-point 1.3T run and both refinements COMPLETE. New numerical comparison
`docs/results/evolution_metal_1300gyr_convergence.json`: 512→1024 L +.01917%,
R -.01022%, He3 -7.51e-5; fourfold tighter timestep L +.000620%, He3 +8.59e-6.
All still fully convective. Same snapshotted executable, old 48-cell gas grid.
Correction to previous timing wording: perf_counter excludes host suspension.
1418 s was awake elapsed; UTC interval ~2973 s. See compact report/docs.

Host sleep interrupted source throughput. AC-only `caffeinate -s -t 7200`
PID 20431, session 17865, started 14:54:37 UTC; expires 16:54 UTC. Stop with
Ctrl-C on session 17865 when finished, renew only if still needed.

Extended gas-grid session 97824 remains active with 8 workers, work
`/tmp/ember-nongrey-extended-v2`, log
`/tmp/ember-nongrey-extended-atmospheres-tail.log`. 48 old + 8 new cells pass.
Cold 2600 K high-gravity initializers are cooling their deep layers slowly;
3000 K initializers approaching physical temperatures; 3200 K still queued.
Only complete 300-row fort.9 iterations give valid correction maxima.
Canonical replay and strict source support remain required.

Solar condensate opacity v3 COMPLETE: 17T from1001..7762 K, 19rho, 30000freq,
SHA07875e83e55e912bd5e527b0df2cf30058292b30c24e904413c825ea0f8f4a05,
`/tmp/ember-condensate-opacity-x700-v3/opacity.bin`. New spec in repo.
Earlier 2600 K condensed model converged but was REJECTED for Tmin1054 K
below old1076 K rectangle (21 optically thin layers); no validated.json.
Do not claim its result accepted. Attempts below1000 K fail original source
partition functions. Exactly1000 K also rounds below source limit; v3 starts1001.
Helium-rich16T table complete at `/tmp/ember-condensate-opacity-x300-y120-v1`,
SHA454f72038e0f3bd38bda4efcea00e5e2bc4e83f4b8fdc26bc5e7fe01af3e215d.
Still needs archive. Solar16T archive is retained separately.

New supported solar condensed atmospheres: session11722 at2600 K,
`/tmp/ember-condensate-atmosphere-cold-v3[.log]`; session60078 at2800 K,
`/tmp/ember-condensate-atmosphere-warm-v3[.log]`. Helium-rich2800 K
condensed control session92684 `/tmp/ember-condensate-atmosphere-helium-v1`,
small per-layer step1.003, corrections slowly decreasing ~.08 at15iterations.
Old solar warm default control25191 was stopped after a persistent2-cycle;
old warm smallstep71956 and coldseed1694 still run but are not convergent yet.
Neither condensed source family nor stellar feedback installed. Need full
source validation, independent chemistry/convective-enthalpy audit, then
boundary comparison. Warm controls may require better nonlinear numerics.

New reproducible optional initializer `scripts/prepare_atmosphere_initializer.py`
scales the whole Newton temperature step with one factor, preserving direction,
instead of independently clipping each layer. Convergence still measures the
UNDAMPED correction; original-source replay required. Prepared condensed
initializer `/tmp/ember-condensate-global-step/prepared.json`; gas initializer
building `/tmp/ember-nongrey-global-step/prepared.json`. New warm2800 condensed
control `/tmp/ember-condensate-atmosphere-warm-global-v3[.log]` just launched.
`run_condensate_atmosphere.py` writes initializer.json, never validated.json,
when prepared source has initialization_only.

Other gas controls: pressure-preserving temperature continuation from2800 K
to2600 K atH.3,Y3=0,g5.15, session16765 work
`/tmp/ember-nongrey-pressure-continuation-cold`, still deep2-cycle.
Same seed plus CONREF everyiteration, session65884,
`/tmp/ember-nongrey-conref-continuous-cold`, only2iterations yet.
T-only hot continuation9895 was stopped for worsening deep corrections;
Pgas was not preserved in that earlier guess. Stops recorded in stopped.json,
SIGTERM only after verifying identified PID command and cwd.

Do not finalize yet: extended72 source grid and coupled condensate controls
remain pending; target next star run2T after grid acceptance. No actual
radiative core reached. All27C++ tests passed, only scripts/docs changed since.

15:47 UTC continuation:
- New mixture guard in src/mixing.cpp rejects different metal inventories
  across separate diffusion regions; no effect on uniform stellar references.
  Ledoux/secular_mixing/coupled_evolution and nongrey_import tests PASS.
- Gas source continuation now reproducible: generator --continuation-models
  selects nearby accepted same-composition profiles, preserving n*k*T and g*m
  while scaling Teff/gravity. New --convective-iterations (default3, up to200)
  supports continuing native CONREF throughout initializer; canonical replay
  remains mandatory. Active generator97824 uses earlier in-memory code.
- Continuous CONREF control65884 is improving: max5.39,2.52,.109,.0742,
  .0526,.0258,.00559 through7 complete iterations. Its dir
  /tmp/ember-nongrey-conref-continuous-cold. Need canonical replay when done.
- Same pressure continuation/continuous CONREF coldg5.4 control55143 started
  /tmp/ember-nongrey-conref-cold-g540; recipe /tmp/ember-conref-cold-g540.py.
- Gas common-step control83593 /tmp/ember-nongrey-global-cold active,
  max5.39→1.11 in2iterations. Condensed common-step control72146 warm2800
  /tmp/ember-condensate-atmosphere-warm-global-v3 active but oscillating.
- Found original OPACTR opacity derivative uses a 1% forward T step. New
  optional initializer --opacity-step .0001 tests a finer numerical Jacobian
  around condensation fronts, original-source final replay still required.
  Builder /tmp/ember-condensate-global-derivative active; no atmosphere yet.
- Old unconverged trials71956,1694,16765 STOPPED after verified PID/cwd checks;
  stopped.json preserves reasons. Continuous CONREF and new controls continue.
- Both new condensed opacity archives COMPLETE: condensation_opacity_x700_extended
  (223artifacts,76MB) and condensation_opacity_x300_he3120 (210,72MB).
  New scripts/archive_condensate_opacity.py makes validation/archive reproducible.
- Independent helium-rich hot join:260states,T3000/4000/5000/6000,
  P1e-13..1000bar, zero condensates; element6.51e-11,pressure2.66e-10,
  nuclei4.32e-10 closure errors. Addedoriginals+report tocondensation_validation.
- Sourcegrid still56/72 accepted, no installationyet. Coldcondv3 correction
  .00696 atiteration6; warmv3 stillcycling. Needcompletedcontrols, independent
  chemistry/convective-enthalpy audit, boundarycomparison, next2T star.

## 16:01 UTC checkpoint

The first coupled condensed atmosphere PASSES source and independent chemistry:
`/tmp/ember-condensate-atmosphere-cold-v3`, solar 2600 K, g5.15, new 17T plane.
11 iterations, correction5.12e-7, flux4.72e-5, hydro1.52e-7, chemical rho2.92e-13.
Full T1054..4528.09 K is supported. Matching T3864, P20770000.
Independent FC element error2.97e-9, all condensing layers Fconv0;
condensation reaches tau.00631. Archived20 artifacts in
`data/atmosphere/sources/condensation_atmosphere_x700_2600_g515`, report
`docs/results/condensation_atmosphere_cold.json`. Original16T gas comparison
T-.1283%, P+.8841% includes material-grid differences: do not attribute all
to condensation. Gas17T solar opacity control81926 now4/17 complete at
`/tmp/ember-condensate-opacity-gas-x700-v3[.log]`,2workers. When complete,
solve gas2600 and2800 on this same grid before isolating the boundary effect.

Gas continuous CONREF cold control65884 COMPLETE14 iterations, correction
2.45e-7, flux9.85e-6, full profile supported. Its misleading validated.json
was RENAMED initializer.json. Independent canonical replay19903 active at
`/tmp/ember-nongrey-conref-continuous-cold-canonical`, log
`/tmp/ember-cold-canonical-replay.log`, recipe same stem.py. Once canonical
passes, can switch remaining grid cells to continuation plus CONREF200.
New reproducible `scripts/run_nongrey_continuation.py` automatically performs
initializer then canonical replay, archives both. Hot3000 K,g5.15,H.3,y0
control35686 active at `/tmp/ember-nongrey-continuation-hot[.log]`.
Coldg5.4 control55143 at4 iterations max.0722 after initial10.2,6.08.

Warm condensed finer-derivative control60487 active at
`/tmp/ember-condensate-atmosphere-warm-derivative-v3[.log]`, prepared
`/tmp/ember-condensate-global-derivative/prepared.json`. OPACTR step1e-4 plus
common Newton limiter; initializer only, canonical replay required.
Through5 iterations max2.13,.366,.18,.217,.0629. Independent analytic
opacity/state probe passes, report under prepared directory
`analytic-opacity-check.json`. The probe script now reads actual source
DELT, tests its secant and state restoration, and reports analytic truncation
separately. Common-step-alone controls83593 gas and72146 condensed STOPPED;
failed trials retained with stopped.json. Original warm condensed v3 session
60078 still cycling. Helium condensed92684 remains active31 iterations
max.0506 under step1.003; no final solution yet.

Independent heldout XH.375,X3.06,T3000,g5.15 specification:
`/tmp/ember-nongrey-heldout-specification.json`. Gas opacity COLD HALF ONLY,
indices0..7 active18675,1worker, `/tmp/ember-nongrey-heldout-opacity`, log
`/tmp/ember-nongrey-heldout-opacity-cold.log`. Need launch hot indices8..15
when a core is available, with same prepared/spec/work, x.375,y.06,modegas.
After both halves finish, rerun without indices to validate caches and merge
opacity.bin. No heldout atmosphere yet; intended to check interpolation in
the new H<.45 interval. Gas-mode adapter identity was validated earlier.

`scripts/nongrey_status.py` now includes correction maxima for only COMPLETE
depth sweeps (last4), avoiding misleading partial fort.9 progress.
All27 tests passed earlier; latest C++ mixture guard passes3 impacted tests;
new continuation hydrostatic-invariant and import tests pass. Python compile
and diff checks pass. Caffeinate17865 expires16:54 UTC; stop on completion.
No2T stellar run launched yet. Extended source grid remains in progress.
