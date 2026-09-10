# Prompt for the next Ember session

Continue the Ember low-mass stellar-evolution project in this repository.
First read `HANDOFF.md` (newest checkpoint first), `docs/EVOLUTION.md`,
`docs/FORWARD_EVOLUTION.md`, `docs/NUCLEAR.md`, `docs/RESTART.md`,
`docs/LIFETIME_SURVEY.md`, and the report
`docs/reports/2026-09-09/ember_status_and_future.tex` or its PDF. The report
README explains its archived histories, exact executable and native checkpoint.
The earlier portions of HANDOFF and ROADMAP are historical, not current job status.

Publication preference: upload code/generators, small configurations/provenance,
and reports, not the bulk generated tables or raw archives. Read
`docs/DATA_REPRODUCTION.md`. Large inputs and report `artifacts/` remain local and
ignored; keep them intact. The small committed `recovery_manifest.json` describes
local recovery files, not files bundled in a fresh clone. Never push the local
full-data backup branch or an ancestor containing the rejected large checkpoint.

The scientific objective is reproducible, numerically converged lifetimes of
low-mass stars as a function of mass and initial composition, followed eventually
through helium-remnant cooling. Do not assume the LBA97 lifetime, a metallicity
maximum, or its location. Distinguish central hydrogen exhaustion, the end of
sustained hydrogen burning, and later cooling endpoints. Specify Y(Z) and metal
pattern when presenting a lifetime(M,Z) result.

Verified baseline: the fixed-baryonic-mass 0.1 Msun, initial XH=.7, X3=0,
Z=.02 GS98 model reaches 2.85 trillion years from its static initial MS model.
It is still fully convective, with XH=.3028875983, X3=.0046633380,
R=.14876088175 Rsun, L=.001857570598 Lsun, Teff=3106.8375 K,
Tc=6.878544 MK, and rhoc=241.23194 g/cm3. The production trajectory uses
the original GS98 EOS/opacity selections and the extended 72-cell GAS atmosphere
family. It includes NO condensate feedback and has not reached core hydrogen
exhaustion. Numerical refinements cover the track through 2T, not a lifetime
or the subtle late radius turnover. The 2.75->2.85T output is only a continuation
segment; the report archives and plots the stitched 0->2.85T history.

Immediate work, in this order unless new evidence changes the priorities:

1. Inspect current files, completed receipts and any surviving controllers before
   starting work. Do not restart old installers or duplicate source calculations.
   Last log inspection (2026-09-09, 23:09 UTC): the TOPS refinement fetch and both
   artificial heat-capacity controls completed; the final solar X3=0, 2800 K,
   logg4.9 atmosphere converged but failed the grain-enthalpy audit. The old
   condensate installer was stopped; its queued stellar controller may still
   wait for `/tmp/ember-condensate-forward-installed.json`, which did not exist.
   Process liveness was not established by that log inspection. The report's
   artifacts preserve the logs and raw refinement sources if /tmp disappears.
2. Finish the hydrogen-poor opacity interpolation validation before selecting it.
   `data/opacity/hydrogen_poor` has support and derivative checks, but the independent
   X=.15,Z=.02 TOPS comparison found 11.9011% worst source-cell discrepancy and
   .918319% in the defined hot rectangle. Do not call this family validated on
   the strength of derivative checks. Completed raw X=.15 at Z=.01/.03 and
   X=.125/.175 at Z=.02 are archived in the report's
   `artifacts/pending_tops_refinement/`, also originally in
   `/tmp/ember-tops-hydrogen-poor-refinement`. Read actual request compositions:
   rounded filenames 012/018 mean .125/.175, not .12/.18. Preserve the coarse
   heldout audit, assemble a separate refined family, and use genuinely independent
   comparisons on the appropriate active opacity/blend domain. Those raw fetches
   have not yet been imported or independently accepted by this handoff.
3. Resolve the physical condensate limitation. Last collection accepted 70/72
   cells; the two solar 2800 K, logg4.9 corners fail the unchanged 1e-8
   convective-flux gate in condensing layers (about 4.1e-7 each). Moving to logg5.0
   did not solve it. Artificial added-cp controls are diagnostics, not a full
   grain EOS or a rigorous error bound. Do not relax the guard merely to obtain
   an installed table. Develop thermodynamically consistent condensate enthalpy
   and convection, or explicitly justify and document a physically controlled
   alternative before changing the production acceptance policy. Grain opacity,
   sizes, settling and equilibrium/rainout assumptions must remain explicit.
4. Extend non-grey atmospheric composition below XH=.3 and extend temperature
   and gravity along the actual stellar path. The separately accepted EOS family
   already reaches XH=.1, but no atmospheric counterpart is installed. Avoid
   extrapolation, clipping or unsupported He I broadening. Choose useful source
   cells based on coverage and physics, not on which corners happen to converge.
   Add explicit opt-in input selection after independent validation. The current
   driver has no new hydrogen-poor opacity selector or warm-block grid format.
5. Continue forward with supported inputs, preserving the gas baseline and
   comparing the effect of any physics changes from a consistent initial model.
   Native restarts require identical executable and physical-input bytes. An
   updated EOS, opacity, atmosphere, or executable cannot silently reuse the old
   checkpoint. Preserve the original d04ae649... executable; the report includes
   its compressed macOS arm64 bytes. Use a separate build for development and
   fresh output paths. See RESTART.md before deciding how to start a changed model.
6. Before radiative-core formation and exhaustion, assess microscopic diffusion,
   gravitational settling, moving-boundary/shell resolution and adaptive mesh,
   semiconvective heat transport, nuclear branch completeness (including pep/CNO
   where relevant), screening/EOS consistency, and full thermal-neutrino losses.
   Add inputs when their relevance is demonstrated; quantify safe omissions.
   Converge endpoint ages, not just luminosities at an early fixed age.

The wider physical agenda from the side discussion is in the report: surface
composition evolution and diffusion-fed residual burning; helium-remnant dense
atmospheres, ion quantum thermodynamics and crystallization; integrated wind and
rotation/magnetic sensitivities; high-Z and helium-enrichment uncertainty;
minimum-mass and higher-mass evolutionary branches. Do not assign C/O-white-dwarf
phase separation, binary-stripped flash thresholds, or helium ignition to the
0.1 Msun star without demonstrating applicability. Environmental heating,
encounters, pycnonuclear reactions and hypothetical proton decay are later,
explicitly conditional extensions.

Work autonomously through justified, reversible development and validation.
Keep progress updates concise. Preserve provenance, failed controls, source-domain
guards, and the distinction between implementation correctness, interpolation
accuracy and physical uncertainty. Do not spend substantial compute repeatedly
checking an already-established fact. Do not modify the historical Henyey/F77
repository. If you publish further work, describe the actual final state and its
limitations rather than claiming a converged stellar lifetime prematurely.
