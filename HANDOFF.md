# Current handoff — September 27, 2026

## September 27 update — 2026-09-27T03:46:28.379892+00:00

The live job is still `out/continuous-lifetime-sept26-v3` (PID 92428):
**1.755 Tyr**, **2993 K**, **5363** accepted intervals,
**0** rejected; still fully convective. Its frozen executable/configuration
are unchanged. The separately prepared `out/continuous-lifetime-sept27-v1`
failed its FIRST 1000-year interval; it accepted no evolved model and is terminal.
Do not report it as a new trajectory or restart it unchanged.

The common driver now selects the restored metal-dependent atmosphere chain
from Z=.02 to .0005, with a smooth overlap onto the PMS/MS boundary. All **5932**
saved early surface queries are exactly unchanged; all **5942** atmosphere
checks pass. The first reference still ends at H=.85. Fable is recovering the
remaining H-rich reference/interval files from retained numerical reports,
privately under `out/fable-atmosphere-recovery-sept27-v1` (P886). Primary owns
selection. Existing high-metal, low-metal and hydrogen-rich interior opacity
approximations are consolidated in `RadiativeOpacity`, including the dense-H
source-slope join and pure-H endpoint. All **34304** layers in the **67** saved
metal-settling checkpoints have numerical support; derivative checks pass.
These source approximations retain their stated physical limits.

An actual radiative-core integration (3.551-Tyr saved star) and exact restart
pass unchanged inventory/energy audits. The useful local composition tolerance
is **1e-12** for stratified diffusion; **1e-15** fails a line search there. The
fresh PMS negative control shows that **1e-12** throughout is too loose during
initial D burning. The driver now uses **1e-15** when the actual instantaneous
mixing partition has one region and the configured tolerance otherwise. It
checks the partition inside every coupling iteration; no phase-dependent
program or relaxed physical audit is introduced. **All 76 tests pass in both
the working and isolated checkout** after this change and opacity integration.

Prepared input set: `out/lifetime-materials-sept27-v1/portable/data/`, **35**
verified files, no external runtime table paths. Its fresh PMS control accepts
five intervals through **7647 yr** with unchanged audits. It is a short control,
not the replacement production sequence. The reference/response sources are
tracked under `data/atmosphere/lifetime` and `data/opacity/lifetime`.

Fable has directly solved pure-H tau=100 columns at log g **6.0** and **6.1**;
independent seeds/depths agree within **0.013% T / 0.054% Pgas**. New temperature
and gravity points are running in `out/fable-highg-pureh-sept27-v1` with at most
four one-thread workers, four aggregate CPU hours, 1 GB output and an 8 GiB
free-disk floor. Check its reservation/receipts before launching. Fable is also
preparing a private **31–32-page** draft candidate under
`docs/research/fable/short_draft_20260927`; the current paper remains 63 pages.
Its corrected conservation report is reviewed but NOT adopted.

Source work is copied into `/private/tmp/ember-reviewed-source-sept26-v1`,
branch `review/end-to-end-sept26`; the new changes are not committed/pushed yet.
Its `origin` is a LOCAL publication checkout, not GitHub. Public target remains
`github.com/oklo/ember`, branch `codex/atmosphere5000-density-checks`; the user
authorizes publication. Next: finish H-rich atmosphere joins, freeze a better
prepared fresh Hayashi sequence, then retire v3 only after its successor passes
initial checks. Preserve all accepted histories and terminal checkpoints.

Detailed evidence: `docs/results/atmosphere_lifetime_atlas_sept27_v1.json`,
`screened_radiative_core_sept27_v1.json`, and
`opacity_lifetime_extension_sept27_v1.json`. The original detailed progress notes are archived locally under
`out/document-history-sept27-v1`; older source history is also retained in Git.

## Objective and protected work

One modern modular program must evolve the same initially 0.1-solar-mass star
from Hayashi contraction through hydrogen burning and WD cooling, then use
appropriate conditional models through proton decay. Also support neighboring
masses and the stellar/brown-dwarf boundary. Do not substitute separate PMS and
post-main-sequence programs. The immediate flash test needs a continuous
atmosphere history; the old switched-boundary tracks are diagnostic comparisons.

The furthest completed continuous Hayashi calculation remains
`out/continuous-lifetime-sept26-v2`: **2.102 Tyr**, **3020 K**, **X=.4170**,
**5931** accepted intervals, **66.16 CPU min / 35.41 min wall**. All accepted
audits pass. Its overly small atmosphere metal allowance caused retries and
termination; see `docs/results/continuous_lifetime_sept26_v2.json`.
The live v3 target is **4.5 Tyr**, four threads, **14400 CPU s** and **18000**
accepted intervals. Read its atomic checkpoint and history before making any
run decision. Never edit a live executable, configuration, or selected table.

Protect the current draft's figure PDFs and `figure_inputs`, active runs and
source workers, `data/production`, `out/binary-eos-sept26-v1`,
`out/atmosphere-lifetime-atlas-sept27-v1/portable`, the new lifetime material
package, both science environments/build checkouts, and all source manifests.
Also preserve `out/convective-transport-sept26-v1`,
`out/pms-mainsequence-v4/seg-008`, the old metal-sequence v48–v67 checkpoints,
selected pulse-front v33–v56 and onset checkpoints/accepted histories,
`out/material-input-recovery-v1`, `out/preserved-material-inputs`,
`out/completed-metal-segments`, and Fable's current review/recovery directories.
Do not resume the obsolete homogeneous pp/D experiment stopped at **415.9 Gyr**.

Cleanup has removed duplicate MESA caches and **230000** obsolete PMS trial
checkpoint copies, preserving all retained segment histories and endpoints.
See `docs/results/obsolete_pms_cleanup_sept27_v1.json`; two unmatched/incomplete
segments remain held. Measured free space after this cleanup is **47.38 GiB**.
The broader cleanup still requires checking figure and physical-input dependencies;
the current pulse comparison scripts read some `out/pulse-*/queries` files.
Do not delete those caches as a group. The trace-metal endpoint array has not
been recovered: protect its existing figure until the paper no longer uses it.

## Reproducibility and coordination

The 248 original runtime inputs are packaged in `data/production` with
`data/MANIFEST.json`; a separate checkout reproduced the full **1020-state**
contraction record and checkpoint payloads exactly. The **960.7 MB** binary EOS
preserves the source doubles, masks and interpolation and reduces loading cost.
The complete test-data bundle is
`out/test-fixture-release-sept26-v1/ember-test-data-complete.tar.gz`, with a
checked installer and manifest. Publication of the large input bundles remains
separate from source publication; their external source terms are documented.

Shared build: `out/unified-evolution-sept25-v1/build`.
Science Python: `out/science-env-sept18-v1/bin/python`.
Writable publication checkout: `/private/tmp/ember-reviewed-source-sept26-v1`.
The original `.git` is read-only and the original worktree has historical changes;
never mass-stage, reset, or remove them. Stage explicit reviewed files in the
publication checkout. The user authorizes GitHub publication; no repeat approval
is needed. Last verified public commit is **3b0d5a4** before this update.

Fable coordination: `docs/research/fable/PRIMARY_NOTES.md` (latest P887),
`HANDOFF_FABLE.md` (read for new deliveries), and
`coordination/primary_review_state.json`. Its CLI watches new P headings.
Check at least hourly while working and resume useful independent assignments
when it recovers from usage limits. All eleven initial reviewers were interrupted;
the main report and focused follow-ups are now complete. Review recommendations
before adopting them. Fable currently owns the bounded atmosphere sources,
private source recovery, and private shortened-paper candidate. Primary owns
all shared-code changes and stellar-run decisions.

The paper is **63 pages**, with four demonstrated errors corrected. The user
requests about **31–32 pages** without losing scientific qualifications, clear
writing without internal research-history jargon, and at most four significant
figures for reported physical quantities. Current draft:
`docs/reports/2026-09-26/ember_status_and_future.pdf`. Keep only the figures used
by the eventual shortened draft, but preserve the present figures and their
inputs until the replacement is verified.
