"""Read native run evidence and give an agent a bounded scientific recovery task.

The controller does not choose new physics or weaken solver tolerances. It
preserves the failed state and asks the agent to implement and test a remedy.
"""
import hashlib
import json
import re
from pathlib import Path


def read_json(path, default=None):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return default


def last_record(path):
    try:
        with path.open("rb") as stream:
            stream.seek(0, 2)
            stream.seek(max(0, stream.tell() - 65536))
            for line in reversed(stream.read().splitlines()):
                try:
                    return json.loads(line)
                except ValueError:
                    continue
    except OSError:
        pass
    return {}


RECOVERY = {
    "table_boundary": (
        "Record the failing T, density, gravity and composition. Search the retained "
        "input catalog and nearby successful tracks first. Check mixture conventions, "
        "phase and overlap values/derivatives before selecting an existing extension. "
        "If needed, build the smallest local extension or bounded analytic continuation "
        "and test it in the overlap. Resume the same accepted checkpoint with explicit "
        "input-identity changes. Do not clip a query to a table edge or repeatedly shrink "
        "the timestep to approach an uncovered state."
    ),
    "convergence": (
        "Locate the first failing equation, zone and phase, using the saved residuals. "
        "Compare nearby successful cases and try a short controlled solver experiment "
        "(predictor, Jacobian, damping or timestep). Keep the final equations and "
        "conservation checks. Repeated rejection of the same state calls for a different "
        "diagnosis, not more identical restarts or a looser acceptance threshold."
    ),
    "conservation": (
        "Preserve the rejected trial and last accepted checkpoint. Reconstruct the "
        "energy/species ledger, including mixing, phase changes, envelope reservoirs "
        "and timestep integration. Test a smaller-step control and repair the specific "
        "accounting/coupling error. Never promote the failed trial as accepted."
    ),
    "segment_limit": (
        "This is a segment budget, not a scientific endpoint. Check the final paired "
        "checkpoint and continue from it. Review whether larger adaptive steps or "
        "existing acceleration controls are justified by the accepted errors."
    ),
    "process_failure": (
        "Read stderr, exit status and the exact launch receipt. Repair the invocation "
        "or resource problem and restart from the last accepted paired checkpoint. "
        "Do not regenerate the star from its initial state."
    ),
    "unclassified": (
        "Read the native report and failing state, classify the obstacle, and run a "
        "small discriminating test. Consult prior successful neighboring calculations. "
        "Implement and validate the next useful step rather than ending with a proposal."
    ),
}


def classify(reason):
    text = reason.lower()
    # A rejection wrapper often contains the more informative underlying cause.
    if any(s in text for s in ("outside source", "outside grid", "outside table",
                               "table boundary", "degeneracy range", "source mask",
                               "outside the", "outside assessed", "out of domain", "not covered")):
        return "table_boundary"
    if any(s in text for s in ("energy audit", "first law", "first-law", "inventory",
                               "conservation", "mass defect", "thermal-reservoir")):
        return "conservation"
    if any(s in text for s in ("newton", "convergence", "singular", "step rejection",
                               "iteration limit", "line search")):
        return "convergence"
    if any(s in text for s in ("cpu budget", "step budget", "wall limit", "requested age")):
        return "segment_limit"
    if any(s in text for s in ("signal", "killed", "exception", "not found", "exit code")):
        return "process_failure"
    return "unclassified"


def observe(work):
    """Read the declared primary case, or the latest launch when none is recorded.

    A report describes that case only. The agent must verify ancestry before
    promoting it to the primary lineage; control runs can have their own ages.
    """
    candidates = []
    for launch in (work / "cases").glob("*/launch.json"):
        case = launch.parent
        if (case / "evolution").is_dir():
            candidates.append((launch.stat().st_mtime, case))
    if not candidates:
        return {"state": "no_native_case", "fingerprint": "no_native_case"}
    _, latest = max(candidates)
    progress = read_json(work / "progress.json", {})
    primary = progress.get("current_case") if isinstance(progress, dict) else None
    # A declaration selects a known launch; it cannot introduce an arbitrary path.
    case = next((p for _, p in candidates if primary in (p.name, f"cases/{p.name}")), latest)
    evolution = case / "evolution"
    report = read_json(evolution / "report.json", {})
    completion = read_json(case / "completion.json", {})
    last = last_record(evolution / "history.jsonl")
    reason = report.get("stop_reason") or completion.get("stop_reason") or ""
    if not reason and completion.get("exit_code") not in (None, 0):
        reason = f"exit code {completion['exit_code']}"
    kind = classify(reason) if reason else None
    checkpoint = None
    for name in ("final.checkpoint", "latest.checkpoint"):
        p = evolution / name
        if p.is_file() and Path(str(p) + ".groups").is_file():
            checkpoint = str(p.relative_to(work))
            break
    normalized = re.sub(r"[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?", "#", reason)
    signature = hashlib.sha256((str(kind) + normalized).encode()).hexdigest()[:16]
    return {
        "case": str(case.relative_to(work)),
        "selection": "declared_primary" if primary in (case.name, f"cases/{case.name}") else "latest_launch",
        "latest_case": str(latest.relative_to(work)),
        "state": "needs_review" if reason else "native_case_open",
        "reason": reason,
        "kind": kind,
        "fingerprint": signature,
        "checkpoint": checkpoint,
        "accepted": {key: last[key] for key in
                     ("years", "Teff_K", "central_T_K", "central_X", "accepted") if key in last},
        "report": str((evolution / "report.json").relative_to(work)) if report else None,
    }


def recovery_prompt(observation, prior):
    kind = observation.get("kind")
    repeats = sum(entry.get("native", {}).get("fingerprint") == observation.get("fingerprint")
                  for entry in prior)
    text = "Latest native evidence (case membership is not proof of primary ancestry):\n"
    text += json.dumps(observation, indent=2) + "\n"
    text += RECOVERY.get(kind, RECOVERY["unclassified"])
    if repeats:
        text += (f"\nThis failure signature has appeared in {repeats} prior agent turns. "
                 "Read their recovery receipts. Identify what changed and why the next "
                 "attempt can succeed; do not repeat a failed intervention unchanged.")
    return text


def operator_stop(work):
    intent = read_json(work / ".supervision/operator.json", {})
    return intent.get("action") == "stop"


CONTROL_INSTRUCTIONS = """
CONTROL PROTOCOL: The operator stops the whole assignment with the supervisor's
`stop` command; .supervision/operator.json belongs to the operator. Never write it
or signal the supervisor yourself. PAUSE_NATIVE pauses stellar jobs while you
diagnose and repair them. A legacy STOP file also pauses native jobs, but does not
stop this agent supervisor. If you created such a native-pause file, remove it
when ready to resume. SESSION_ENDING is a supervisor handoff request: finish the
current bounded work, preserve paired checkpoints, and write NEXT.md. Reconcile
all registered cases before launching replacements. Never leave unregistered
long-lived processes. Register native jobs under cases/<name>/launch.json with
PID and exact argv, then write completion.json when they exit.

SCIENTIFIC RECOVERY: Table coverage and convergence failures are repair tasks,
not evidence of a physical endpoint. Use the retained research record and nearby
successful trajectories. Keep source/configuration hashes and preserve accepted
state and inventories. Numerical convergence, physical validity and input-domain
coverage are separate checks. A new approximation must be named, bounded and
tested; do not turn off a heating source or weaken a scientific audit to finish.
Update progress.json with the actual state, and recovery.json with diagnosis,
evidence, intervention, comparison, outcome and next action. If an agent turn
ends unfinished, the supervisor will resume it within the original limits.
"""
