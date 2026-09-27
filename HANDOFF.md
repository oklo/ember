# Current handoff — September 27, 2026

## Active work at 2026-09-27T04:22:51.937480+00:00

**Primary run:** `out/continuous-lifetime-sept27-v2`, PID **6843**, caffeinate
**6844**, fresh Hayashi origin and exact five-interval restart prefix. It is
at **52.16 Gyr**, **1097** accepted / **0** rejected. Target
4.5 Tyr, 25000 steps, 21600 CPU seconds, four threads. Frozen binary/input
hashes and the exact command are in its `PLAN.json`. Its self-contained
`data/MANIFEST.json` verifies **50** files. Do not edit these inputs.

The selected common code has initial D, pp+closed CN with SFIII rates,
plasma losses, hot screened H/He/metal transport, partition-dependent local
composition convergence and unchanged physical audits. Restored interior
opacity and metal/H-dependent atmosphere sources are selected from the start.
The atmosphere queries cover all 5932 earlier PMS states exactly and old
metal models v1–v45 through 3.749 Tyr; v46–v67 remain unsupported pending
trace-He/pure-H integration. Receipt:
`docs/results/atmosphere_hydrogen_recovery_sept27_v1.json`.

**Parallel fixed-Z run:** `out/continuous-lifetime-sept26-v3`, PID92428,
now **2.967 Tyr**, **3114 K**, **7099** accepted / **0** rejected.
Its 14400-CPU-second cap will stop it normally soon; preserve the completed
history. It uses the old local1e-15 setting and fixed-Z atmosphere allowance.
Do not mistake its larger age for a result with the restored metal atlas.

**Fable:** session limit until 03:50 America/New_York (07:50 UTC). Primary
owns atmosphere recovery and shorter-paper corrections (P890/P891). Source
stages3/4 ended unaccepted; no active source column remains. The 4800K
logg6.0/6.1 source columns pass; cooler cases oscillate and 6.2 remains
unsupported. Do not duplicate exhausted attempts without diagnosing them.

**Paper:** private candidate under `docs/research/fable/short_draft_20260927`
has been corrected and visually checked at **30 pages** (from63). All paper
figures are still unchanged. Peak luminosity values now match the actual
pulse figure: 666/848/1695 grids give5.537e36/4.835e36/3.281e35erg/s,
not a converged peak. Publication of the shorter September27 paper remains
pending. Code through d67b4b6 is public; H-interval/recovery work is copied
into the isolated publication checkout and **all76 isolated tests pass**.
Next: commit/publish those reviewed changes, publish the shorter paper, then
finish trace-He and hydrogen-dominated atmosphere integration while the star
runs. Do not stop active jobs simply because this handoff is being updated.

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
