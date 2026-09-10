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
[the scientific report](docs/reports/2026-09-09/ember_status_and_future.pdf)
([LaTeX](docs/reports/2026-09-09/ember_status_and_future.tex),
[recovery and validation](docs/reports/2026-09-09/README.md)). This checkpoint
supersedes the status of the historical entries below. The user requested a
report, handoff and GitHub checkpoint before clearing context. No further
stellar evolution or source-physics change was made during this packaging step.

- Furthest accepted stellar model is still **2.85 trillion years**, gas-only,
  original GS98 EOS/opacity inputs, 512 points, fully convective. XH=.3028875983,
  X3=.0046633380, R=.14876088175 Rsun, L=.001857570598 Lsun,
  Teff=3106.8375 K, Tc=6.878544 MK, rhoc=241.23194 g/cm3. The atmosphere
  floor remains XH=.3. No hydrogen-exhaustion or lifetime claim.
- Report artifacts archive the original 0->2.5T and both restart segments,
  receipts, final native checkpoint and original macOS arm64 d04ae649... binary.
  Their hashes are in `docs/reports/2026-09-09/artifacts/manifest.json`.
  The stitched 0->2.85T figure checks matching seam states. Preserve the original
  binary: changed executable/physics bytes are incompatible with native restart.
- Low-H EOS is available separately through XH=.1 and passed its source checks.
  Low-H opacity remains **under validation**: independent midpoint discrepancy
  11.9011% over all eligible source cells; .918319% in the audit's hot rectangle.
  These are interpolation checks, not physical uncertainty estimates or direct
  errors on the existing star. Current stellar inputs remain unchanged.
- TOPS refinement controller completed all FOUR raw requests: X=.15 at
  Z=.01/.03 and X=.125/.175 at Z=.02. Raw requests/results and completion receipt
  are preserved under `docs/reports/2026-09-09/artifacts/pending_tops_refinement/`.
  They are NOT imported/independently accepted. Read actual request mixtures;
  filenames x012/x018 are rounded labels. Retain the old coarse heldout audit.
- Last collector accepted70/72 condensate cells. Both remaining solar2800K,
  logg4.9 cells now converge but fail the unchanged grain-enthalpy gate:
  X3=.12 gives4.169866e-7 and X3=0 gives4.113889e-7 convective flux fraction
  in condensing layers, above1e-8. The fine-initializer controller exited on
  that physical audit failure; it is no longer merely a slow source solve.
  Both logg5.0 controls also failed. No production policy relaxation.
- Both artificial extra-cp controls completed. For cp+1e6/+1e7 erg/g/K,
  fractional boundary pressure changes are7.4938e-7/7.3755e-6. They remain
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
  cells. Maximum discrepancy11.9011% atT=.002keV,rho1.5849g/cm3; maximum
  hot-rectangle error.918319% atT=.03keV,rho.15849. Report
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
- Final XH=.3028875982546768, X3=.004663338019452628,
  X4=.6724490637258706, Z=.02; R=.14876088174924135Rsun,
  L=.0018575705980597992Lsun, Teff3106.83751480519K,
  Tc6878544.13316771K, rhoc241.23194367687486g/cm3,
  Pc1.6524980327519686e17dyn/cm2, logg5.093090139533256.
  Fully convective/homogeneous; minimum diffusive/adiabatic ratio3.19869346,
  maximum local conductive fraction.03288828, reduced ppII/pp.00296634.
  Nuclear rest-mass release fraction.00282709; gravity remains fixed baryonic.
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
  E.10866%, cp.19344%; new XH<.3 only P.03491%,E.09190%,cp.15484%.
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
  chemistry, but grain-bearing layers carry4.169866e-7 of total flux by
  convection, above unchanged1e-8 guard. grain_enthalpy_supported=false.
  Work /tmp/ember-condensate-forward-retries/focused-plane-005-g490/plane-005/teff-2800-gravity-000.
  Archived diagnostic in sources/condensation_grain_transport_rejections/
  solar_he3120_2800_g490; report condensation_grain_transport_rejection.json.
- Tested g5.0 because it brackets the actual stellar trajectory. BOTH solar
  He3 controls still fail grain guard:6.231848e-7 (He3=.12) and7.415962e-7
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
  cp1e6 complete: delta P7.4938e-7, T-1.5554e-7, rho1.15105e-6 relative.
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
- Final: XH=.31976956298654213, X3=.005907354423384522,
  X4=.6543230825900733, Z=.02; R=.14885436723893694,
  L=.001825613477418555, Teff3092.416617007074K,
  Tc6720471.802976825K, rhoc240.35825208978363g/cm3,
  Pc1.6446842720858246e17dyn/cm2, logg5.092544465428368.
  Fully convective and homogeneous. Min diffusive/adiabatic ratio3.65651,
  maxlocal conductive fraction.0322201, plasmon/L1.412159e-8,
  maxreducedppII/pp.002253555, maxactiveppzeta.05596646.
- Reports docs/results/{evolution,transport,nuclear}_metal_m010_2750gyr_gas.json;
  evolution PDF/PNG generated and PNG visually checked. This history/plot
  contains ONLY2.5->2.75T. Original0->2.5T remains in2500gyrgas artifacts;
  full-history He3 peak.1002124 at.710749T is in that earlier report.
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
  XH=.31977 is close to that boundary. Do not extrapolate or weaken guards.

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
- FURTHEST COMPLETE:2.5T GAS, session73844 completed~19:56.1542accepted,
  one rejection. XH=.35995867683748073, X3=.010152583279596424,
  X4=.6098887398829226, R=.1487487109232278, L=.0017530076547748674,
  Teff3062.287264314398K, Tc6374821.889243306K, rhoc239.8557623365566,
  Pc1.6405161792728554e17, logg5.093161205780678. Fully convective.
  Min diffusive/adiabatic ratio4.97484; maxlocal conduction3.092686%;
  plasmon/L1.349805e-8; maxppzeta.0591508; maxreducedppII/pp.00120310.
  Integrated nuclear rest-mass release .0023918915 of fixed baryonic mass.
  Output out/evolution-metal-512-2500gyr-gas.json SHA
  3051b3801f1592787184b09fab58033a21a25d438306f600422daa1631a75bfb.
  Same original gas executable b904...; no condensate stellar feedback.
  Reports/plots docs/results/{evolution,transport,nuclear}_metal_m010_2500gyr_gas;
  evolution png/pdf visually checked. No2.5T numerical refinements yet.
- TWO-TRILLION-YEAR refinements40162 now ALL COMPLETE. Full repeat is byte
  identical. Mesh512->1024: L+.0616321%, R-.0123403%, XH-1.214223e-4,
  X3-5.280259e-5.4x tighter time: L-.00361005%, R-.000594084%,
  XH+1.300106e-5, X3+2.216246e-6. Summary rerun with both comparisons;
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
  report condensation_warm_material_refinement.json; T+.0327057%,P+.288987%.
  Fine33T model REMAINS INELIGIBLE for production (grainenthalpy false).
  Both artificial cp controls complete+archived;1e7 changes matching pressure
  4.387998e-6, density6.444117e-6;1e6 pressure8.388176e-7. Audits check each
  archived source/input/output and same reference, require<1e-5 boundary
  differences. No physical uncertainty bound claimed. Production1e-8 grain
  guard unchanged. Thermal/controller26103,78272,76806 allfinished.
- Heldout2775 COMPLETE source+independentchemistry,grainflux0,deepestcond
  tau.00079676. Archive/heldout comparison will happen under30387 oncefullgrid.
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
  independent chemistry flags Fconv/Ftotal=1.224272682e-6 in one condensing
  outer layer, above the strict production1e-8 gate. Old controller77176
  exited1 as intended. Tmatch4030.3551233K, Pmatch16463729.2576. Numerical
  coarse/fine differences T+.03270566%, P+.28898665%, rho+.24322585%.
  It is NOT a production atmosphere. Archived explicitly as diagnostic at
  sources/condensation_material_refinement/refined-2800-diagnostic.
- Native convection root cancellation was tested and ruled out: report
  condensation_convection_roundoff.json, relative roots agree4.44e-16.
  Diagnostic source/archive retains printing-only changes; no production
  Fortran formula changed. Actual gas cp at flagged face154500891erg/g/K.
- NEW estimate_condensate_enthalpy.py uses pinned FastChem logK temperature
  derivatives and independent fixed-P chemistry. Report
  condensation_warm_latent_enthalpy_estimate.json: at native face1769.41749K,
  P4061.88979dyn/cm2, added ideal-atom-referenced cp~288772erg/g/K
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
  coarse/fine T+.0265223%,P+.3605097%,rho+.3323842%. This is a numerical
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
  XH=.4308797380, X3=.02574632445, R=.1478253064, L=.001631953818,
  Teff3017.374 K, Tc5.834241 MK, rhoc242.4931. Still fully convective;
  minimum diffusive/adiabatic gradient ratio8.2433, max conduction flux2.9244%.
  Full output out/evolution-metal-512-2tyr-gas.json SHA
  17e5c5ebdb7d9b8acfc2f7b17e2d467d2445c0db7ca181678298d9cc119d8d47.
  Snapshot binary SHA b904902f4f86843dc957f3b004f1db78acfa983f93d370a083eecb581f00500c.
  UTC2365.262s, child user1629.476s/system2.144s, concurrent atmosphere jobs.
  Transport and nuclear audits COMPLETE. New summarize_forward_evolution.py
  produced docs/results/evolution_metal_m010_2tyr_gas.{json,png,pdf}; figure
  visually checked. Generic utility verifies receipts, age, same binary/data/
  physical selections, and mesh/time arguments; comparisons await refinements.
- Refinements40162 ACTIVE since18:36, two workers: mesh1024 and fresh512repeat;
  tight512 follows repeat. At18:51 repeat~.52T, mesh~.31T. Same copied binary.
  The repeat must be byte-identical; no completed2T refinement claim yet.
- Indexed CONREF control15531 COMPLETE and archived. Canonical final relative
  T-2.20825e-9, P+1.20904e-8, rho+1.73761e-8 versus original2750 control.
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
  deepest-layer total flux residual .0030436 exceeds .002. Retained source
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
  original capacity21 coarse: T+.026503%,P+.360625%,rho+.332540%.
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
  initial2767.52K. The attempted126-cell2600..3200 rectangle remains rejected
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
  Matched17T gas comparison P+.124609%,T-.0234705%,rho+.182795%.
  Archive sources/condensation_atmosphere_x700_2800_g515 and report
  docs/results/condensation_atmosphere_warm.json. New
  compare_condensate_control.py verifies source/input/opacity physics and
  revalidates compressed gas/condensed originals. Figure now3columns,
  docs/results/condensation_atmosphere_controls.{png,pdf}, visually inspected.
- Heldout2900 COMPLETE, no condensates in its actual profile, independent
  element closure1.77e-13.2725 initializer still runs under31059. Independent
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
  star track starts2767.52K. NEW warm recovery45169 runs bothfailedplanes at
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
  Canonical2 iterations afterCONREF14. Interpolation errors T+.19584%,
  P+.03401%,rho+.12946%. Original opacity and atmosphere/initializer artifacts
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
  Isolated depletion effect: Pgas+.8095%,T-.1349%,rho+1.1038%.
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
  144 independent source checks: max P .061635%, E .16747%, cp .25943%.
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
R=.14422884460648513, L=.0013686915161790307, Teff2923.32373045,
XH=.5331385347948866, X3=.09105648824134456, Tc5153914.71795,
rhoc257.7343772393. Fully convective. Uses old48 non-grey gas grid,
GS98 EOS/screening/conduction/opacity conversion, Ledoux. Run predates the
EOS optimization and slow diffusion (zero coefficients, fully convective).
Exact executable hash/wall time were NOT recorded at launch: do not invent
those. Further runs must snapshot executable and provenance before launch.

Extended source work:
- `/tmp/ember-nongrey-extended-v2`: all6 opacity planes complete,16T x19rho
  x30000 frequencies, H=.3/.45/.7 and X3=0/.12. Spec in data/atmosphere/sources.
  T1075.69..7762.30 retains old low-T isotherms exactly. No extrapolation.
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

At 1e12 yr: R=.14357836 Rsun, L=.0013302286 Lsun, Teff=2909.13 K,
Tc=5.12268e6 K, rhoc=261.626 g/cm³, XH=.53389752, X3=.09509109 and
X4=.35101139. The star remains fully convective. He3 peaks at .10256641
near 736.81 Gyr. There are 906 accepted macrosteps and one rejected attempt;
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
at 1075.69 K. No atmosphere or interior source extrapolation is enabled.
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

At **2048 points and 1e12 yr**: R=.14351457 Rsun, L=.0013495294 Lsun,
Teff=2920.28 K, Tc=5.13439e6 K, rhoc=261.915 g/cm³, X=.53254997,
He3=.09343622. The run accepts 926 macrosteps and rejects three attempts.
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
the same command now gives R=.12935124, L=.00094042956, Teff=2810.43 K.
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

At 4096 points: 12 updates, residual 4.10338e-10, undamped correction
3.06322e-10, R=.12860617 Rsun, L=.00097562953 Lsun, Teff=2844.5719 K,
Tc=4.5663102e6 K, rhoc=366.87522 g/cm³. Nuclear luminosity agrees with
surface L to ~2e-15 relative; independent virial error is -6.5848e-7.
The unresolved center is 1.23e-10 of the mass. All profiles are ordered.

| Points | R/Rsun | L/Lsun | Teff (K) | Absolute virial error |
|---:|---:|---:|---:|---:|
| 128 | .12911938 | .0009629408 | 2829.64 | 6.856e-4 |
| 256 | .12935265 | .0009420455 | 2811.62 | 1.699e-4 |
| 512 | .12885921 | .0009637381 | 2833.08 | 4.229e-5 |
| 1024 | .12862818 | .0009746770 | 2843.63 | 1.055e-5 |
| 2048 | .12861667 | .0009751461 | 2844.10 | 2.635e-6 |
| 4096 | .12860617 | .0009756295 | 2844.57 | 6.585e-7 |

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
forced to zero. Max |D|=.24178 in the 4096-point star at T=27304 K,
rho=.12394 g/cm³, m/M=.99999310; mass RMS=.0043487. Small outer mass is
not evidence of small influence on R or L.

`python3 scripts/audit_cms19_energy.py /tmp/ember-eos2019.tar.gz` reproduces
both defects directly from original source columns:

- Pure H at log T=4.45, log P[GPa]=1.35 has D=-.26498098, before ember
  interpolation. This is a source-version consistency issue as well as an
  interpolation accuracy question.
- Pure He at rho=1 g/cm³ drops in internal energy from 8.88587e13 to
  8.22697e13 erg/g between T=891251 and 1e6 K. Secant dU/dT=-6.05890e7
  erg/g/K while listed entropy cv=+9.72478e7.

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
failure diagnostic. The .5 Msun case has R=.94219963 Rsun, L=.019877124
Lsun, Teff=2232.77 K; its inflated cool envelope is not realistic. The n=3
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

14:29 UTC: 512point1.3T run84690 COMPLETE,1417.86085wallsec underconcurrent
sourcework. FinalR.14620742152,L.0014647156256,Teff2953.1139155,
X.505654541594,He3.070718566033,Tc5306449.997,rhoc248.3491230.
Compactreport/completeinputhashes+transport/nuclearaudits in
`docs/results/*metal_m010_1300gyr.json`; fullout andoriginalreceipt preserved.
Mintransportgrad/ad12.989, maxconductiveflux2.837%,ppzeta.07284.
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
1417.9 s was awake elapsed; UTC interval ~2973 s. See compact report/docs.

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
Earlier 2600 K condensed model converged but was REJECTED for Tmin1053.993 K
below old1075.69 K rectangle (21 optically thin layers); no validated.json.
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
Full T1054.09..4528.09 K is supported. Matching T3863.670177, P20772434.3287.
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
