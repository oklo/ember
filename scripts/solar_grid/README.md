# Solar evolution grid

The solar-v1 campaign follows 30 masses from 0.01 to 0.30 solar masses at initial Z = 0.016482259798651158. The release board also includes the retained 0.07145 and 0.07150 solar-mass comparison histories. Their source-limited endpoints and temperature reversals are preserved separately from the 30 requested 1 K crossings.

[Release assumptions](solar-release-v1.json) specify the birth clock, surface boundary, nominal cold source, conditions for a reduced continuation and evidence for unusual phases. Freeze this file when issuing a release. Record the implemented nuclear rate, material data, transport and surface closures in each reviewed history. Sensitivities qualify the nominal result. Failed domain, continuity or energy checks require repair.

## Running an existing campaign

Run `python3 scripts/solar_grid/coordinate.py`. Set `EMBER_GRID_ROOT` to select another campaign. The controller needs the campaign manifest, local solver, physics inputs, a bounded recovery plan and its referenced scripts. The public repository contains orchestration and retained comparisons. Large campaign inputs and run receipts live on the operator machine.

A campaign opts into release routing with `release-config.json` containing `release_sha256`, the SHA-256 of the exact release specification. The coordinator reconciles running workers before dispatch. A runnable native segment reserves a lane. Ready recovery jobs inherit the priority of the work orders that need them, with the existing plan order breaking ties. An independent assessment can gate the pilot that uses a table. Interrupted atmosphere columns can resume from a changed finite profile. `PAUSE_NATIVE`, `STOP`, the operator stop and the disk reserve hold new launches.

## Decisions and progress

[Work orders](work-orders.json) name the affected masses, the question, the concrete jobs and the criterion for adoption. To extend the campaign, place the complete revised list in its `work-orders.json`. Each target must exist in the recovery plan. Its dependencies identify required table work, controls and evolution. Use a boolean evidence gate for a scoped pilot. Set `restart_on_adoption` to existing downstream job IDs when a production restart requires a reviewed decision. Review targets must precede that restart. A decision binds its evidence hashes, the work-order definition and the release.

Record an adjudication with `python3 scripts/solar_grid/adjudicate.py CAMPAIGN ORDER adopted --reviewer NAME --rationale TEXT EVIDENCE.json`. Use `sensitivity` or `rejected` for alternatives and failed proposals. Specify evidence paths relative to the campaign. Changed evidence reopens review. The decision history is retained in `decisions/`.

The scheduler heartbeat and workers are in `status.json`. Read `release-status.json` for each mass's birth-clock frontier, separate dispatch state, endpoint review, and next work order. Scientific completion requires the independent endpoint verifier, reviewed initial contraction and a matching `release_sha256` in the accepted endpoint record. The comparison archive is checked against its pinned hash. To audit the underlying history and control evidence, run `python3 scripts/solar_grid/verify_endpoints.py CAMPAIGN`.

## Low-cost shepherd

Run `python3 scripts/solar_grid/shepherd.py CAMPAIGN` for a board refresh. For model-assisted prioritization, create `shepherd-config.json` with `enabled` set to true and `binary` set to the absolute Codex CLI path. The coordinator normally refreshes every five minutes. An independent timer can run the same command and set `independent_shepherd` in the release configuration. Its lock serializes observations.

The agent sees a compact queue summary and can nominate an existing work order. The dispatcher applies the preference for up to six hours while its exact job remains ready. The scheduler’s resource eligibility and lane limits take precedence. The shepherd uses gpt-6-luna at low reasoning effort in an ephemeral, read-only Codex session. It receives the queue summary and returns a schema-constrained nomination. A 120-second timeout and two attempts per rolling day limit usage. The CLI reports token usage when available. The call ledger reserves each attempt before invoking the CLI. Evidence and cooldown suppress repeated requests. Read `shepherd/luna-calls.json` for recommendations, flaws and token usage, and `shepherd/health.json` for observation failures. Set `enabled` to false to continue deterministic operation without model calls.

For the upper grid the controller expands the existing composition-edge recipe in memory. A saved same-mass trial must pass its ancestry, convergence and inventory checks. A new independent column is measured at the actual crossing and compared with a replayed table prediction. A passing review makes the interval eligible for the birth-clock history. The first new 0.11 crossing review is a prerequisite for launching additional new columns in that family. Existing named jobs retain their receipts and limits.

Inspect `dispatch-review.json` for the scheduler’s eligible jobs and `release-dispatch.jsonl` for actual launches. A model nomination can wait behind running work or an eligible job of a preferred resource class.
