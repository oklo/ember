# Asynchronous collaboration: primary agent and Fable

This protocol coordinates work through files while either agent may be busy,
paused or recovering context. It does not resume an agent or bypass an app's
token limits. Scientific jobs can continue through a pause only if they were
actually launched with a bounded execution and output plan. Silence is not
failure, completion, consent or a reason to duplicate a job.

## Ownership and communication

The primary owns ASSIGNMENT.md, PRIMARY_NOTES.md, COORDINATION.md and shared
Ember code, physical inputs, HANDOFF.md and the paper. Fable owns STATUS.md,
HANDOFF_FABLE.md, QUESTIONS.md, FINDINGS.md and its survey code/results. Each
agent writes only its own files. Read-only access to the other agent's files
is expected. Fable's permitted paths remain those in ASSIGNMENT.md.

Messages in PRIMARY_NOTES.md have IDs P001, P002, etc.; messages in QUESTIONS.md
have IDs F001, F002, etc. Include UTC time, task ID, whether action is required,
and the concrete request or result. Acknowledge a message by ID in your own
file, with the action taken or reason it is pending. Preserve unacknowledged
messages across edits. Acknowledgment means receipt, not scientific acceptance.
Mark superseded instructions by ID and retain enough history to understand them.

Replace status and handoff files atomically: write a temporary file in the same
directory and rename it after closing. Keep STATUS.md short; durable detail
belongs in the handoff or versioned results. Neither agent edits the other's
acknowledgments. Do not place transient chat assumptions only in a long log.

## Starting or resuming

1. Read ASSIGNMENT.md and this protocol, then your own latest handoff, the other
   agent's messages, and STATUS.md. A primary-agent context recovery also reads
   the root HANDOFF.md. The current written assignment supersedes older chat.
2. Inspect actual processes, run receipts and resource reservations. Check
   process identity, not just a PID that might have been reused. A stale status
   timestamp does not establish that a job stopped.
3. Acknowledge new messages by ID. Preserve completed results and active jobs.
   Resume analysis or monitoring rather than rerunning an identical case.
4. Write status: UTC time, active task IDs, last acknowledged message, actual
   jobs, completed results, resource use, next action and any dependency.

Suggested task IDs are F-LIT-001 for a literature result, F-EQ-001 for an
equilibrium, F-CTRL-001 for a control and F-SPH-001 for an encounter. Assign new
IDs for changed physical or numerical inputs; record which case is repeated.
The primary's stellar tasks have a separate E- prefix. These IDs are labels,
not evidence that any job has started.

## Work and resource reservations

Fable may select and execute the next useful task within its remit without
waiting for the primary to acknowledge every result. Changes to the shared
stellar code, accepted physics or paper still belong to the primary. When a
dependency is pending, pursue an independent literature, analytic or control
task. A completed batch is a point to assess results, not a mandatory pause.

Before launching a new batch, each agent creates its own JSON reservation in
`coordination/reservations/` beneath this directory using an exclusive create
operation. Include owner, task IDs, UTC time, requested threads, expected memory,
data cap and scratch paths. Inspect all other reservations and actual processes
again before launch; this second check catches simultaneous reservations. If
they conflict, the main Ember/atmosphere work takes priority and Fable revises
its reservation before launch. Two reservations must not both assume all cores
are free. After launch, the owner records actual PIDs, start times and commands.

Existing stellar and atmosphere jobs may predate this protocol, so process
inspection is always required; absence of reservations does not mean idle.
Record bounds in each launcher and allow resource use to change at batch
boundaries. The primary also honors active SPH reservations rather than
blindly oversubscribing; arrange a reduction at a safe checkpoint if necessary.
Never kill the other agent's job. Remove or mark released only your own
reservation after verifying completion. A reservation does not expire merely
because its owner is paused. If the process is demonstrably absent, note the
evidence in your own message and exclude it from current resource accounting;
do not alter the other owner's record.

For every launched job, keep a machine-readable receipt with task ID, inputs
and executable hashes, exact command, working directory, start time and process
identity, resource limits, output paths, and eventual exit status and checks.
Record failures as failures. A log with no live process and no final receipt
is an interrupted or unverified job until inspected, not an accepted result.
Do not duplicate an active input/executable combination under a fresh name.

## Pause and context recovery

Before a foreseeable token-limit or context pause, Fable updates HANDOFF_FABLE.md
and STATUS.md with: last completed result, exact next command or decision,
pending message IDs, actual live jobs and reservations, locations of controls
and source changes, numerical concerns, and compute/storage already spent.
The primary records the analogous information in root HANDOFF.md and a short
message here when it affects Fable. Put new information first.

Label the agent paused separately from each job's state. A bounded batch with
recorded limits may finish while its agent is paused. Its completion receipt
allows the other agent to inspect it. Do not launch an indefinite background
queue that relies on an unavailable agent for cleanup or storage control.
If a pause is abrupt, recovery uses receipts and process identity; it does not
assume that the last intended action happened.

## Publishing a result for review

Fable keeps each reviewable batch in `results/<batch-id>/`. At completion write
`READY.json` last, after its referenced files are closed and validated. Include
task IDs, scope, hashes of compact inputs/results/source changes, checks passed
and failed, unresolved limitations, and compute/storage use. Large scratch
files may be listed separately; label what is retained and what is reproducible
but deliberately removed. Do not require the primary to inspect a live HDF5 file.

Once READY.json is written, treat that batch as immutable. Corrections use a
new batch ID and state what they supersede. The primary acknowledges the ready
batch in PRIMARY_NOTES.md, then records accepted, provisional or needs changes
with reasons. Fable can continue independent work while review is pending.
Only the primary promotes results into the paper or accepted Ember inputs.
Acceptance of a simulation result does not imply acceptance of its speculative
environmental assumptions; report those separately.

## Rhythm and recovery from disagreement

Read messages at startup/resume, before a batch, after completion, and when
changing direction. During long tool waits, check at the next safe opportunity;
do not busy-poll files or interrupt useful computation for constant chatter.
Update status when there is a finding, a job-state change, a blocker or a pause.

If findings conflict, preserve both results, identify the differing input or
definition, and propose the smallest discriminating check. Do not overwrite the
other result or silently move a tolerance. If the main task is blocked on a
reply, continue another task within the remit. Ask the user only when their
scientific preference, missing access or a scope decision is actually needed.

The collaboration is working if it produces independently interpretable
results, recovers from pauses without rerunning completed work, uses available
compute without oversubscription, and keeps the main trajectory progressing.
