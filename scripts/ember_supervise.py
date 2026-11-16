#!/usr/bin/env python3
"""Supervise an Ember research agent through native failures and bounded recovery.

Uses Claude Max CLI today. Native jobs retain their own scientific checks; this
process supplies durable intent, recovery context, resource limits and handoffs.
"""
from datetime import datetime, timezone
import argparse
from ember_recovery import observe, recovery_prompt, operator_stop, CONTROL_INSTRUCTIONS
from pathlib import Path
import fcntl
import hashlib
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import time


def now():
    return datetime.now(timezone.utc).isoformat()


def read_json(path, default=None):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return default


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def file_in_work(work, value):
    path = Path(value)
    path = (path if path.is_absolute() else work / path).resolve()
    if not path.is_relative_to(work) or not path.is_file():
        raise ValueError(f"Missing evidence or path outside study: {path}")
    return path


def terminal_event(path):
    result = None
    with path.open() as handle:
        for line in handle:
            try:
                event = json.loads(line)
                if event.get("type") == "result":
                    result = event
            except ValueError:
                pass
    return result


def subscription_auth(claude):
    """Require the user's Max login; never fall back to API billing."""
    reply = subprocess.run([claude, "auth", "status"], capture_output=True,
                           text=True, timeout=20)
    try:
        auth = json.loads(reply.stdout)
    except ValueError as error:
        raise RuntimeError("Claude authentication status unavailable; no API fallback") from error
    if (reply.returncode or not auth.get("loggedIn")
            or auth.get("authMethod") != "claude.ai"
            or auth.get("subscriptionType") != "max"):
        raise RuntimeError("Claude Max subscription login required; no API fallback")
    return {key: auth[key] for key in ("authMethod", "subscriptionType")}


def estimated_session_cost(work, session_id):
    """Claude reports a cumulative list-price estimate, not a Max bill."""
    estimates = [0.0]
    for path in (work / "sessions").glob("*/completion.json"):
        result = read_json(path, {}).get("claude_result") or {}
        if result.get("session_id") == session_id:
            value = result.get("total_cost_usd")
            if isinstance(value, (int, float)) and math.isfinite(value) and value >= 0:
                estimates.append(value)
    return max(estimates)


def provider_usage(path):
    """Read provider limits separately from its unused overage allowance."""
    latest = {}
    error = None
    with path.open() as handle:
        for line in handle:
            try:
                event = json.loads(line)
            except ValueError:
                continue
            if event.get("type") == "rate_limit_event":
                latest = event.get("rate_limit_info") or {}
                if latest.get("status") in ("allowed", "allowed_warning"):
                    error = None
            # Also handle a limit error when no structured usage event arrived.
            if (event.get("type") == "assistant" and event.get("error") == "rate_limit"
                    or event.get("type") == "result" and event.get("is_error")):
                detail = json.dumps(event.get("message", event.get("errors", event.get("result", ""))))
                if (event.get("error") == "rate_limit"
                        or re.search(r"rate[_ -]?limit|usage limit|hit your limit", detail, re.I)):
                    error = detail[:2000]
    return latest, error


def usage_stop_reason(usage, error):
    if usage.get("isUsingOverage"):
        return "subscription usage credits active; no further continuation"
    windows = usage.get("unifiedWindows") or {}
    if (usage.get("status") == "rejected"
            or any(window.get("status") == "rejected" for window in windows.values())
            or error):
        return "Claude subscription usage limit; see provider_usage for reset time"
    # overageStatus=rejected alone means paid overflow is disabled, not that
    # the included subscription allowance has been exhausted.
    return None


def completion_decision(work, config):
    result = read_json(work / "result.json", {})
    if result.get("status") != "achieved":
        return False, "scientific objective unfinished"
    try:
        contract = config["completion"]
        for key, expected in contract.get("expected", {}).items():
            actual = result.get(key)
            if actual != expected:
                raise ValueError("wrong result identity: " + key)
        history = file_in_work(work, result["native_history"])
        if history.stat().st_size == 0:
            raise ValueError("empty native history")
        validation = read_json(file_in_work(work, result["validation"]), {})
        required = contract["checks"]
        if not required:
            raise ValueError("completion requires scientific checks")
        for key in required:
            if validation.get(key) is not True:
                raise ValueError("missing assessment: " + key)
        if "temperature_max_K" in contract:
            teff = float(result["endpoint_Teff_K"])
            if not math.isfinite(teff) or teff <= 0:
                raise ValueError("invalid endpoint temperature")
            equilibrium = contract.get("allow_equilibrium", False) and validation.get("equilibrium_supported") is True
            if teff > contract["temperature_max_K"] and not equilibrium:
                raise ValueError("endpoint or supported residual equilibrium not reached")
        evidence = result.get("evidence", [])
        if len({e["path"] for e in evidence}) < 2:
            raise ValueError("need at least two evidence artifacts")
        for e in evidence:
            if hashlib.sha256(file_in_work(work, e["path"]).read_bytes()).hexdigest() != e["sha256"]:
                raise ValueError("evidence hash mismatch")
        return True, "scientific completion receipts verified; scientific review still required"
    except (OSError, ValueError, KeyError, TypeError) as error:
        return False, str(error)


def continuation_prompt(work, config):
    prior = []
    for path in sorted((work / "sessions").glob("*/recovery-context.json")):
        prior.append(read_json(path, {}))
    return (
        "Continue the authorized scientific assignment under CURRENT brief.txt. "
        "Read USER_UPDATE.md, NEXT.md, progress.json, recovery.json and the last "
        "session receipt. Reconcile existing jobs before restarting.\nObjective: "
        + config["objective"] + "\n" + CONTROL_INSTRUCTIONS + "\n"
        + recovery_prompt(observe(work), prior)
        + "\nCompletion assessment: " + completion_decision(work, config)[1]
    )


def provider_reset(usage, current_time):
    """Return the latest rejected-window reset; disabled paid overflow is not a limit."""
    resets = []
    windows = list((usage.get("unifiedWindows") or {}).values()) + [usage]
    for window in windows:
        if window.get("status") == "rejected":
            t = window.get("resetsAt")
            if isinstance(t, (int, float)) and math.isfinite(t) and t > current_time:
                resets.append(t)
    return max(resets) if resets else None


def stop_child(child):
    if child.poll() is None:
        try:
            os.killpg(child.pid, signal.SIGTERM)
            child.wait(timeout=20)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGKILL)
            child.wait()
        except ProcessLookupError:
            child.wait()


def stop_registered_cases(work):
    # Only stop a PID whose live command still matches a launch receipt here.
    for receipt in (work / "cases").glob("*/launch.json"):
        if (receipt.parent / "completion.json").exists():
            continue
        data = read_json(receipt, {})
        pid = data.get("pid")
        command = data.get("command", [])
        if not pid or not command:
            continue
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            continue
        live = subprocess.run(["/bin/ps", "-p", str(pid), "-o", "command="],
                              capture_output=True, text=True).stdout.strip()
        expected = command if isinstance(command, str) else " ".join(command)
        if live == expected:
            try:
                os.kill(pid, signal.SIGTERM)
                write_json(receipt.parent / "supervisor-stop.json",
                           {"utc": now(), "pid": pid, "reason": "session ended"})
            except ProcessLookupError:
                pass


def work_bytes(work):
    return sum(path.stat().st_size for path in work.rglob("*")
               if not path.is_symlink() and path.is_file())


def main(work, resume=False, wait_lock=False, config_name="control.json"):
    config_path = work / config_name
    config = read_json(config_path)
    lock = (work / "supervisor.lock").open("a")
    while True:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            break
        except BlockingIOError:
            if not wait_lock:
                raise RuntimeError("Another supervisor owns this study")
            if operator_stop(work):
                return
            time.sleep(min(5, config["poll_seconds"]))
    resume = resume or (work / "RESUME_SUPERVISOR.json").exists()
    if (work / "completion.json").exists() and not resume:
        raise RuntimeError("Finished assignment: prepare a new receipt before relaunching")
    previous = read_json(work / "supervisor.json", {}) if resume else {}
    elapsed = (datetime.now(timezone.utc) - datetime.fromisoformat(previous["started_utc"])).total_seconds() if previous else 0.0
    started = time.monotonic() - elapsed
    estimate = estimated_session_cost(work, config["session_id"])
    first_invocation = max((int(path.parent.name) for path in (work / "sessions").glob("*/launch.json")), default=0) + 1
    if resume:
        receipt = read_json(work / "RESUME_SUPERVISOR.json", {})
        receipt.update(resumed_utc=now(), elapsed_seconds=elapsed,
                       estimated_api_cost_usd=estimate, billing_mode="max_subscription")
        write_json(work / "supervisor-resume-receipt.json", receipt)
        (work / "RESUME_SUPERVISOR.json").unlink(missing_ok=True)
        if (work / "completion.json").exists():
            (work / "completion.json").replace(work / "supervisor-before-update.json")
    failures = 0
    stop = "session limit"
    observation = observe(work)
    observation_key = None
    observation_changed = time.monotonic()
    replan_requested = None
    child = None
    guard = None
    requested_stop = False
    state = {"started_utc": previous.get("started_utc", now()), "supervisor_pid": os.getpid(),
             "session_id": config["session_id"], "phase": "starting",
             "billing_mode": "max_subscription", "estimated_api_cost_usd": estimate,
             "cost_estimate_is_bill": False}

    def signal_stop(_signum, _frame):
        nonlocal requested_stop
        requested_stop = True

    signal.signal(signal.SIGTERM, signal_stop)
    signal.signal(signal.SIGINT, signal_stop)
    caffeinate = shutil.which("caffeinate")
    if config.get("caffeinate", True) and caffeinate:
        guard = subprocess.Popen([caffeinate, "-i", "-w", str(os.getpid())])
        state["guard_pid"] = guard.pid

    try:
        for number in range(first_invocation, config["max_sessions"] + 1):
            # User corrections are live; resource usage is carried across supervisor updates.
            config = read_json(config_path, config)
            if completion_decision(work, config)[0]:
                stop = "scientific completion receipts verified"
                break
            if requested_stop or operator_stop(work):
                stop = "operator stop" if operator_stop(work) else "supervisor termination requested"
                break
            if time.monotonic() - started >= config["wall_seconds"]:
                stop = "assignment wall limit"
                break
            if shutil.disk_usage(work).free < config["minimum_free_bytes"]:
                stop = "disk reserve"
                break
            if work_bytes(work) >= config["max_work_bytes"]:
                stop = "study disk allowance"
                break
            state["authentication"] = subscription_auth(config["claude"])
            directory = work / "sessions" / f"{number:03d}"
            directory.mkdir()
            (work / "SESSION_ENDING").unlink(missing_ok=True)
            observation = observe(work)
            prompt = (((work / "brief.txt").read_text() + "\n" + CONTROL_INSTRUCTIONS)
                      if number == 1 else continuation_prompt(work, config))
            write_json(directory / "recovery-context.json", {"utc": now(), "native": observation})
            replan_requested = None
            observation_key = None
            observation_changed = time.monotonic()
            (directory / "prompt.txt").write_text(prompt + "\n")
            command = [config["claude"], "--print", "--model", config["model"],
                       "--session-id" if number == 1 else "--resume", config["session_id"],
                       "--permission-mode", "auto", "--permission-prompts", "none",
                       "--restricted", "--tools", "Read,Glob,Grep,Edit,Write,Bash",
                       "--add-dir", config.get("reference_root", str(work)), "--strict-mcp-config",
                       "--mcp-config", '{"mcpServers":{}}', "--output-format", "stream-json",
                       "--verbose"]
            session_start = time.monotonic()
            interrupted = None
            disk_check = 0.0
            native_check = 0.0
            with (directory / "prompt.txt").open() as inp, (directory / "stream.jsonl").open("x") as out, (directory / "stderr.log").open("x") as err:
                child = subprocess.Popen(command, cwd=work, stdin=inp, stdout=out,
                                         stderr=err, start_new_session=True)
                write_json(directory / "launch.json", {"utc": now(), "pid": child.pid,
                           "supervisor_pid": os.getpid(), "command": command})
                while child.poll() is None:
                    config = read_json(config_path, config)
                    elapsed = time.monotonic() - started
                    session_elapsed = time.monotonic() - session_start
                    state.update(utc=now(), phase="working", invocation=number,
                                 child_pid=child.pid, wall_seconds=elapsed,
                                 estimated_api_cost_usd=estimate, stream=str(directory / "stream.jsonl"))
                    if time.monotonic() - native_check >= config.get("native_poll_seconds", 30):
                        native_check = time.monotonic()
                        observation = observe(work)
                        state["native"] = observation
                        key = (observation.get("case"), observation.get("fingerprint"),
                               observation.get("accepted", {}).get("years"))
                        if key != observation_key:
                            observation_key, observation_changed = key, time.monotonic()
                            replan_requested = None
                        stagnant = time.monotonic() - observation_changed
                        # Only interrupt repeated diagnosed failures. A long advancing
                        # native solve or a source build is not itself a failed task.
                        if (observation.get("kind") in ("table_boundary", "convergence", "conservation", "process_failure")
                                and stagnant > config.get("diagnosis_seconds", 1800)):
                            if replan_requested is None:
                                replan_requested = time.monotonic()
                                write_json(work / "REPLAN.json", {"utc": now(), "native": observation,
                                           "reason": "same native failure without accepted progress"})
                                (work / "SESSION_ENDING").touch()
                    write_json(work / "supervisor.json", state)
                    if requested_stop or operator_stop(work):
                        interrupted = "operator stop" if operator_stop(work) else "supervisor termination requested"
                    elif elapsed >= config["wall_seconds"]:
                        interrupted = "assignment wall limit"
                    elif session_elapsed >= config["session_wall_seconds"]:
                        interrupted = "session wall limit; automatic continuation"
                    elif replan_requested is not None and time.monotonic() - replan_requested > config.get("handoff_seconds", 60):
                        interrupted = "recovery handoff; automatic continuation"
                    if min(config["wall_seconds"] - elapsed,
                           config["session_wall_seconds"] - session_elapsed) < 300:
                        (work / "SESSION_ENDING").touch()
                    if time.monotonic() - disk_check >= 60:
                        disk_check = time.monotonic()
                        if shutil.disk_usage(work).free < config["minimum_free_bytes"]:
                            interrupted = "disk reserve"
                        elif work_bytes(work) >= config["max_work_bytes"]:
                            interrupted = "study disk allowance"
                    if interrupted:
                        stop_child(child)
                        break
                    time.sleep(config["poll_seconds"])
            stop_registered_cases(work)
            result = terminal_event(directory / "stream.jsonl")
            cost = result.get("total_cost_usd") if result else None
            if isinstance(cost, (int, float)) and math.isfinite(cost) and cost >= 0:
                estimate = max(estimate, cost)
            usage, usage_error = provider_usage(directory / "stream.jsonl")
            state.update(provider_usage=usage, provider_limit_error=usage_error)
            write_json(directory / "completion.json", {"utc": now(), "exit_code": child.returncode,
                       "wall_seconds": time.monotonic() - session_start, "stop_reason": interrupted or "turn ended",
                       "estimated_session_api_cost_usd": cost,
                       "cost_estimate_is_bill": False, "provider_usage": usage,
                       "provider_limit_error": usage_error, "claude_result": result})
            config = read_json(config_path, config)
            achieved, explanation = completion_decision(work, config)
            state["endpoint_receipt_check"] = explanation
            if achieved:
                stop = "scientific completion receipts verified"
                break
            if interrupted and "automatic continuation" not in interrupted:
                stop = interrupted
                break
            limit_reason = usage_stop_reason(usage, usage_error)
            if limit_reason:
                reset = provider_reset(usage, time.time())
                remaining = config["wall_seconds"] - (time.monotonic() - started)
                if reset is None or reset - time.time() >= remaining or usage.get("isUsingOverage"):
                    stop = limit_reason
                    break
                state.update(phase="waiting_for_subscription_reset", reset_utc=datetime.fromtimestamp(reset, timezone.utc).isoformat())
                while time.time() < reset + 2:
                    state["utc"] = now()
                    write_json(work / "supervisor.json", state)
                    if requested_stop or operator_stop(work) or time.monotonic() - started >= config["wall_seconds"]:
                        break
                    time.sleep(min(5, max(.1, reset + 2 - time.time())))
                failures = 0
                continue
            error = child.returncode != 0 or not result or result.get("is_error", False)
            failures = failures + 1 if error else 0
            if failures >= 3:
                stop = "three consecutive CLI failures; inspect session receipts"
                break
            state.update(phase="continuing", utc=now(), estimated_api_cost_usd=estimate)
            write_json(work / "supervisor.json", state)
            # Short, interruptible retry backoff; no daily blind authentication retries.
            for _ in range(6 if error else 1):
                if requested_stop or operator_stop(work):
                    break
                time.sleep(config.get("retry_seconds", 5))
    except Exception as error:
        stop = f"supervisor error: {type(error).__name__}: {error}"
        raise
    finally:
        if child is not None:
            stop_child(child)
        stop_registered_cases(work)
        if guard is not None:
            guard.terminate()
            guard.wait(timeout=10)
        state.update(utc=now(), phase="finished", stop_reason=stop,
                     wall_seconds=time.monotonic() - started, estimated_api_cost_usd=estimate)
        write_json(work / "supervisor.json", state)
        state["native"] = observe(work)
        write_json(work / "completion.json", state)
        lock.close()
        print(json.dumps(state), flush=True)


def cli():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("run", "status", "stop"))
    parser.add_argument("study", type=Path)
    parser.add_argument("--resume", action="store_true", help="resume within the original assignment window")
    parser.add_argument("--wait-for-handover", action="store_true", help="wait for the previous supervisor to release its lock")
    parser.add_argument("--config", default="control.json", help="study-relative control file; useful during a controller handover")
    parser.add_argument("--reason", default="operator request")
    args = parser.parse_args()
    work = args.study.resolve()
    if not work.is_dir():
        parser.error("study directory does not exist")
    if args.action == "stop":
        (work / ".supervision").mkdir(exist_ok=True)
        write_json(work / ".supervision/operator.json", {"action": "stop", "utc": now(), "reason": args.reason})
        return
    if args.action == "status":
        state = read_json(work / "supervisor.json", {})
        state["native"] = observe(work)
        stamp = state.get("utc")
        state["heartbeat_stale"] = bool(stamp and (datetime.now(timezone.utc) - datetime.fromisoformat(stamp)).total_seconds() > 60 and state.get("phase") != "finished")
        print(json.dumps(state, indent=2))
        return
    main(work, resume=args.resume, wait_lock=args.wait_for_handover, config_name=args.config)


if __name__ == "__main__":
    cli()
