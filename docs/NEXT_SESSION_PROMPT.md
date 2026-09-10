# Continue the 0.1 solar-mass remnant calculation

Read the newest entry in `HANDOFF.md`, `docs/COLD_REMNANT.md`, and the
[September 10 working paper](reports/2026-09-10/ember_status_and_future.pdf).
Inspect processes and completed receipts before launching anything. Reuse
productive work; old running markers alone do not establish liveness.

The user has authorized autonomous work through an evolved helium-remnant
cooling branch to Teff=100 K. Preserve that objective across status questions
and documentation updates. A cooling-law extrapolation or unsupported material
table does not meet it. Distinguish core hydrogen exhaustion, extinction of
sustained burning, and the 1000/500/100 K cooling milestones and recrossings.

The completed reference reaches **3.40 trillion years**, still fully convective:
XH=.2034, He3=.0009703, R=.1459 Rsun,
L=.002016 Lsun, Teff=3202 K, Tc=8.031 MK,
rhoc=258.3 g/cm3. It uses the 64-plane `exhaustion_refined_v2` EOS,
`hydrogen_poor_refined_v4` opacity, 96-node `exhaustion_warm_v1` gas atmosphere,
Ledoux/slow mixing and HRW plasma-neutrino losses. Its fresh 0–3.30T run
and exact 3.30–3.40T restart are complete and archived locally with receipts.

The 108-node candidate adds twelve XH=.1 warm atmosphere nodes. Source,
depth, EOS and runtime audits pass. The independent XH=.15 heldout differs
by +2.134% in matching T; completed controls attribute most of it to
hydrogen interpolation. New XH=.15 nodes and independent XH=.125/.175
heldouts have explicit plans and controllers. A separate diagnostic fresh
512-point stellar run targets 3.60T with the candidate. Inspect its actual
status; do not promote it to a physically converged result merely because
the integrator completes. Follow the latest handoff for exact paths.

Next, finish that refinement and compare the resulting stellar boundary and
transition. The selected reference atmosphere ends at XH=.2, so blindly
continuing the old checkpoint would soon leave coverage. Native restarts
require identical executable, tables, tolerances and physical selections.
Use separate builds, source work directories and output paths for changes;
never bypass restart identity or replace the original executable/input files.

The remaining programme includes zero-H atmosphere normalization, condensate
thermodynamic consistency, microscopic diffusion/settling and energetics,
moving-boundary/shell resolution, nuclear branch and screening audits, full
thermal-neutrino losses, and dense helium EOS/phase/transport/boundary physics.
Read the explicit acceptance order in `docs/COLD_REMNANT.md`. The electron
response/entropy and isolated Ioffe derivative work are preparation, not an
installed remnant EOS. A native opacity normalization test also rejects simply
rescaling source abundances without checking absorption coefficients.

The separate F77 comparison locates the AJR central transition at 3.224T,
XH=.01756. Its age cannot be equated to Ember's age at XH=.203. Consult
`docs/COMPUTATIONAL_COST.md` and the retained comparison reports; do not
reinterpret the F77 result as an original-paper or modern-physics benchmark.
Do not edit the historical Henyey/F77 repository.

Publish source, generators, small configurations/provenance and reports;
new bulk tables, raw archives, native checkpoints and executables stay local.
See `docs/DATA_REPRODUCTION.md`. September 9 local recovery artifacts remain
intact even though September 10 replaces that paper in the current tree.
Never push a full-data backup branch or force-add ignored bulk data.

Keep concise progress updates, provenance and failed controls. Run checks
appropriate to each change and measure missing physics along actual profiles.
Do not repeatedly spend compute re-establishing a result already verified.
