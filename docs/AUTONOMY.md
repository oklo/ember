# Autonomous calculations

`scripts/ember_supervise.py` keeps an authorized scientific assignment running
across agent turns. The initial backend is the Claude CLI with a Max subscription;
it does not fall back to API billing. The supervisor uses Python's standard
library and runs on macOS or Linux.

The stellar solver decides whether a step is acceptable. The agent diagnoses
failures, implements changes and tests them. The supervisor preserves the
assignment, supplies recovery context and checks completion receipts. These are
separate responsibilities: an agent finishing a response is not a completed
stellar calculation.

## Start and monitor

Create a study directory under `out/`, copy
[`examples/agent-control.json`](../examples/agent-control.json) to `control.json`,
and replace its session ID with a UUID. Set the objective, reference directory,
resource limits and completion checks for the calculation. Write `brief.txt`
with the physical assumptions, initial inputs, permitted approximations and
required comparisons. The agent owns the study directory; shared inputs and
accepted reference calculations should remain unchanged.

```sh
python3 scripts/ember_supervise.py run out/my-study
python3 scripts/ember_supervise.py status out/my-study
python3 scripts/ember_supervise.py stop out/my-study --reason 'review accepted models'
```

The run command stays in the foreground and can be managed by a terminal or
service manager. `--resume` continues an interrupted assignment within its
original wall-time allowance. It does not silently extend the deadline or clear
an operator stop. To resume after an intentional stop, archive
`.supervision/operator.json`, then run with `--resume`. A process lock prevents
duplicate supervisors; `--wait-for-handover` waits for the previous owner.

## Pausing and stopping

| Signal | Owner | Effect |
| --- | --- | --- |
| `.supervision/operator.json`, written by `stop` | Operator | Stops the whole assignment |
| `PAUSE_NATIVE` | Agent | Pauses stellar jobs so the agent can diagnose them |
| `STOP` | Legacy native runner | Treated as a native pause; does not stop this supervisor |
| `SESSION_ENDING` | Supervisor | Requests a saved checkpoint and handoff before the next agent turn |

Native runners must honor the native-pause signals and save accepted state.
The agent removes its own pause when ready to continue. Separating these controls
prevents a diagnostic pause from terminating the agent that should repair the run.
Older study-specific supervisors may interpret `STOP` differently; migrate at a
checkpoint and update the study instructions before changing controllers.

## Recovery

Native jobs register their exact argv and PID in `cases/<name>/launch.json`,
write `completion.json` on exit, and retain the usual `evolution/report.json`,
`history.jsonl` and checkpoint/group pair. Set `progress.json`'s `current_case`
to the primary case name, updating it when the lineage advances. The supervisor
selects that case only if it has a launch receipt, then reads its native artifacts.
Without a matching declaration it observes the latest launch. The newest case
may be a control: its age alone does not establish primary-lineage membership.

`ember_recovery.py` distinguishes input-domain failures, convergence problems,
conservation failures, segment budgets and process failures. Its recovery order
is to inspect the failing state, search existing compatible inputs and nearby
successful calculations, test a small intervention, and resume the accepted
lineage. Table extensions require overlap and composition checks. Numerical
changes retain the final scientific acceptance checks. The controller never
clips a table query, changes composition or turns off a source automatically.

Each resumed turn receives the native failure, accepted state, checkpoint path
and previous occurrences of the failure signature. Agents record their diagnosis,
intervention, comparison and outcome in `recovery.json`. A diagnosed failure
without accepted progress triggers a handoff after `diagnosis_seconds`; this
does not interrupt a long native solve merely because it is slow. The agent
should check for `SESSION_ENDING` between bounded operations.

An unfinished agent response is resumed automatically. Three consecutive CLI
errors stop with a receipt. A reported subscription reset inside the original
window is waited out, with a distinct waiting status and an interruptible stop.
Disk, session and assignment limits remain in force. Authentication failures or
unknown resource limits are reported, not bypassed.

## Completion and research record

`result.json` uses `status: "achieved"` only when the requested trajectory exists.
It names `native_history`, `validation`, and at least two evidence files with
SHA-256 hashes. The completion contract in `control.json` specifies identity
fields and required true checks in the validation file. Optional temperature
criteria can allow a separately demonstrated residual-heating equilibrium.
Missing checks, a wrong mass or changed evidence prevent completion.

These receipts support review; they do not independently prove the physics.
Keep the commands, input/source identities, compact history, selected accepted
checkpoints and comparison results needed to reproduce the calculation. Large
regenerable diagnostics belong in the study workspace, not in the package.

Regression tests use real subprocesses with synthetic agent and stellar outputs:

```sh
python3 -B tests/test_agent_supervisor.py -v
```
