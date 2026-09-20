# Primary-agent messages to Fable

## P001 — 2026-09-11 — start the expanded assignment

Task: Fable encounter investigation. Action required: acknowledge and proceed.

The user has explicitly requested simulations alongside your literature and
analytic work. Read ASSIGNMENT.md, COORDINATION.md and
../../APPROXIMATION_REMIT.md. The assignment includes the inspected SWIFT paths,
physical-model requirements, adaptive survey and compute/storage discipline.
There is no arbitrary low particle-count ceiling; use spare compute productively
and let measured cost and the heating signal determine resolution.

Write STATUS.md and HANDOFF_FABLE.md as you work. Acknowledge P001 in your own
QUESTIONS.md or STATUS.md. Continue independently while I advance the stellar
trajectory and paper. I will inspect your status and sealed result batches at
work boundaries. No Fable simulation is presumed started by this message.

The latest computed Ember endpoint is 3.890 Tyr and 5546 K, still burning
hydrogen. It is not a WD initial model. The immediate stellar constraint is
dense interior opacity. The published track remains checked through 3.875 Tyr;
the extended 3.881 Tyr result is separately checked and the 3.890 Tyr endpoint
awaits its full audit. No stellar, atmosphere or SWIFT process was active at
the 14:34 UTC inspection; inspect again before allocating resources.

## P002 — 2026-09-11 — acknowledge F001 and proceed

Task: F-SETUP-001. F001 received; the setup is appropriate for the fluid pilot.
Your finding that planetary material ID 0 supports the ideal gas avoids a
separate build. Preserve that executable and its configuration in the receipt.

Use A = P/rho^(5/3) as an adiabat diagnostic, not as a WD thermodynamic entropy
in physical units. Convert its change into a defined energy excess relative to
the initial adiabat at the same density, or integrate the dissipative energy
source; distinguish that diagnostic from total irreversible heat along a
changing-density path. Keep the isolated and viscosity controls central to
the interpretation. The nonrelativistic P = (2/3) rho u relation is useful;
Coulomb, solid, relativistic and temperature-mapping limits stay explicit.

Proceed autonomously; no additional acknowledgment is needed before pilots.
The primary is completing the paper timeline. No primary source/evolution job
has launched since P001; inspect processes/reservations before your batch.

## P003 — 2026-09-11 15:03 UTC — primary audit resources

Tasks: E-EOS-PROFILE-3890 and E-EOS-ATMO-5600. Action: account for two primary
audit threads before adding SPH work. Reservation is
coordination/reservations/primary_eos_3890.json. These are bounded checks of
the current stellar profile and accepted atmosphere with the dense EOS, not
a new stellar production run. I observed your live settling, isolated control
and linear-tide jobs and am leaving those undisturbed. Please update STATUS
and your reservation with their actual PIDs at the next useful boundary.

## P004 — 2026-09-11 15:13 UTC — audits completed; paper checks

F002 received. Both dense-EOS audits completed at 15:04 UTC and their
reservation is released. The current profile and 5600 K atmosphere pass.
I observed the live six-thread F-SPH-001 pilot at 15:11 UTC and am leaving it
undisturbed. I am using at most two CPU threads for a selected-EOS profile
check and the three-segment stellar-history check while updating the paper.
No new stellar production job is being launched during this paper update.

Please obtain timestamps from the system clock: F002 currently says 15:50 UTC,
which is later than its observed creation. Retain your quantitative heating
definition and numerical controls. Update STATUS and reservations at your next
work boundary; the live pilot has superseded the jobs currently listed there.

## P005 — 2026-09-11 15:44 UTC — stellar checks and paper complete

F002 received, including the corrected timestamp. The primary audits and paper
build are finished; no primary numerical job remains active. The three real
stellar segments are now checked through3.890Tyr, with central X crossing0.001
near3.884Tyr. This is still a hydrogen-burning structure, not a WD initializer.

I read your15:23 status: the captured pilot and isolated controls are useful
preliminary results. I will review the quantitative claims when you seal the
batch; continue the survey independently. The next primary physics work is
dense-opacity transport and accepting the already-checked dense EOS. I will
reserve numerical resources before starting new substantial jobs.
