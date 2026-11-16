"""Exercise supervisor recovery with actual subprocesses and synthetic stellar receipts."""
import hashlib
import importlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
supervisor = importlib.import_module("ember_supervise")
recovery = importlib.import_module("ember_recovery")

FAKE_AGENT = r'''
import hashlib,json,os,sys,time
from pathlib import Path
w=Path(__file__).resolve().parent
mode=(w/'mode').read_text()
if sys.argv[1:2]==['auth']:
 print(json.dumps(dict(loggedIn=True,authMethod='api' if mode=='api' else 'claude.ai',subscriptionType='max')))
 sys.exit(0)
prompt=sys.stdin.read()
n=len(list(w.glob('invoked-*')))+1
(w/f'invoked-{n}').write_text(prompt)
print(json.dumps(dict(type='system',subtype='init',model='test')),flush=True)
if mode=='pause':
 (w/'STOP').touch()
 (w/'PAUSE_NATIVE').touch()
 time.sleep(.25)
elif mode=='wait':
 while True:time.sleep(.05)
elif mode=='limit' and n==1:
 print(json.dumps(dict(type='rate_limit_event',rate_limit_info=dict(status='rejected',resetsAt=time.time()+.2))),flush=True)
 print(json.dumps(dict(type='result',is_error=True,result='usage limit')),flush=True)
 sys.exit(0)
elif mode=='repair' and n==1:
 p=w/'cases/a/evolution';p.mkdir(parents=True)
 (p.parent/'launch.json').write_text(json.dumps(dict(pid=os.getpid(),command=['fixture'])))
 (p/'report.json').write_text(json.dumps(dict(stop_reason='repeated step rejection: Teff=2599.95 outside source grid')))
 (p/'history.jsonl').write_text(json.dumps(dict(years=380691998,Teff_K=2600.2))+'\n')
 (p/'final.checkpoint').write_text('accepted physical state')
 (p/'final.checkpoint.groups').write_text('accepted composition')
 print(json.dumps(dict(type='result',is_error=False,result='Needs atmosphere repair')),flush=True)
 sys.exit(0)
(w/'history.jsonl').write_text('native history\n')
(w/'validation.json').write_text(json.dumps(dict(ancestry=True,energy=True)))
r=dict(status='achieved',mass_Msun=.08,native_history='history.jsonl',validation='validation.json',endpoint_Teff_K=4,
 evidence=[dict(path=f,sha256=hashlib.sha256((w/f).read_bytes()).hexdigest()) for f in ['history.jsonl','validation.json']])
(w/'result.json').write_text(json.dumps(r))
print(json.dumps(dict(type='result',is_error=False,result='Complete')),flush=True)
'''


class SupervisorTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="ember-supervisor-")
        self.work = Path(self.temp.name).resolve()
        (self.work / "sessions").mkdir()
        (self.work / "cases").mkdir()
        (self.work / "brief.txt").write_text("Evolve the star, preserving its accepted state.")
        self.agent = self.work / "agent"
        self.agent.write_text("#!" + sys.executable + "\n" + FAKE_AGENT)
        self.agent.chmod(0o755)
        self.config = dict(claude=str(self.agent), model="test", session_id="fixture",
                           objective="test recovery", wall_seconds=40, session_wall_seconds=20,
                           max_sessions=4, poll_seconds=.02, native_poll_seconds=.02,
                           minimum_free_bytes=0, max_work_bytes=10000000, caffeinate=False,
                           retry_seconds=.01, completion=dict(expected=dict(mass_Msun=.08),
                           checks=["ancestry", "energy"], temperature_max_K=5))
        self.process = None

    def tearDown(self):
        if self.process and self.process.poll() is None:
            self.process.terminate()
            self.process.wait(timeout=10)
        self.temp.cleanup()

    def start(self, mode, extra=()):
        (self.work / "mode").write_text(mode)
        (self.work / "control.json").write_text(json.dumps(self.config))
        self.process = subprocess.Popen([sys.executable, str(ROOT / "scripts/ember_supervise.py"),
                                         "run", str(self.work), *extra],
                                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    def finish(self):
        text, _ = self.process.communicate(timeout=20)
        self.assertEqual(self.process.returncode, 0, text)
        return json.loads((self.work / "completion.json").read_text())

    def test_unfinished_turn_resumes_with_failure_and_checkpoint(self):
        self.start("repair")
        result = self.finish()
        self.assertEqual(result["stop_reason"], "scientific completion receipts verified")
        second = (self.work / "invoked-2").read_text()
        self.assertIn("table_boundary", second)
        self.assertIn("final.checkpoint", second)
        self.assertIn("overlap", second)
        self.assertEqual((self.work / "cases/a/evolution/final.checkpoint").read_text(),
                         "accepted physical state")

    def test_agent_native_stop_does_not_kill_supervisor(self):
        self.start("pause")
        result = self.finish()
        self.assertEqual(result["stop_reason"], "scientific completion receipts verified")

    def test_operator_stop_is_authoritative(self):
        self.start("wait")
        deadline = time.monotonic() + 5
        while not (self.work / "invoked-1").exists():
            self.assertLess(time.monotonic(), deadline)
            time.sleep(.02)
        subprocess.run([sys.executable, str(ROOT / "scripts/ember_supervise.py"), "stop", str(self.work)], check=True)
        self.assertEqual(self.finish()["stop_reason"], "operator stop")

    def test_controller_handover_preserves_session_and_deadline(self):
        self.start("wait")
        deadline = time.monotonic() + 5
        while not (self.work / "invoked-1").exists():
            self.assertLess(time.monotonic(), deadline)
            time.sleep(.02)
        original = json.loads((self.work / "supervisor.json").read_text())["started_utc"]
        duplicate = subprocess.run([sys.executable, str(ROOT / "scripts/ember_supervise.py"),
                                    "run", str(self.work)], capture_output=True, text=True)
        self.assertNotEqual(duplicate.returncode, 0)
        self.assertIn("Another supervisor", duplicate.stderr)
        (self.work / "handover.json").write_text(json.dumps(self.config))
        replacement = subprocess.Popen([sys.executable, str(ROOT / "scripts/ember_supervise.py"),
                                        "run", str(self.work), "--resume", "--wait-for-handover",
                                        "--config", "handover.json"],
                                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        (self.work / "mode").write_text("pause")
        self.process.terminate(); self.process.communicate(timeout=10)
        self.process = replacement
        result = self.finish()
        self.assertEqual(result["stop_reason"], "scientific completion receipts verified")
        self.assertEqual(result["started_utc"], original)
        self.assertTrue((self.work / "invoked-2").exists())

    def test_subscription_reset_resumes_without_paid_overflow(self):
        self.start("limit")
        result = self.finish()
        self.assertEqual(result["stop_reason"], "scientific completion receipts verified")
        self.assertTrue((self.work / "invoked-2").exists())

    def test_api_auth_never_launches_agent(self):
        self.start("api")
        text, _ = self.process.communicate(timeout=10)
        self.assertNotEqual(self.process.returncode, 0)
        self.assertIn("no API fallback", text)
        self.assertFalse((self.work / "invoked-1").exists())

    def test_original_deadline_survives_restart(self):
        start = "2000-01-01T00:00:00+00:00"
        (self.work / "supervisor.json").write_text(json.dumps(dict(started_utc=start)))
        (self.work / "completion.json").write_text("{}")
        self.start("wait", ("--resume",))
        self.assertEqual(self.finish()["stop_reason"], "assignment wall limit")
        self.assertFalse((self.work / "invoked-1").exists())

    def test_completion_requires_intact_evidence(self):
        self.start("pause"); self.finish()
        self.assertTrue(supervisor.completion_decision(self.work, self.config)[0])
        (self.work / "history.jsonl").write_text("changed")
        self.assertFalse(supervisor.completion_decision(self.work, self.config)[0])

    def test_wrong_mass_cannot_complete(self):
        self.start("pause"); self.finish()
        p = self.work / "result.json"; result = json.loads(p.read_text())
        result["mass_Msun"] = .1; p.write_text(json.dumps(result))
        self.assertFalse(supervisor.completion_decision(self.work, self.config)[0])

    def test_new_case_supersedes_old_failure(self):
        for name, reason in [("old", "outside source grid"), ("new", None)]:
            p = self.work / "cases" / name / "evolution"; p.mkdir(parents=True)
            (p.parent / "launch.json").write_text("{}")
            if reason: (p / "report.json").write_text(json.dumps(dict(stop_reason=reason)))
            time.sleep(.01)
        self.assertEqual(recovery.observe(self.work)["state"], "native_case_open")
        self.assertEqual(recovery.observe(self.work)["case"], "cases/new")

    def test_table_error_inside_rejection_is_not_solver_failure(self):
        self.assertEqual(recovery.classify("repeated step rejection: electron density outside pair-table degeneracy range"), "table_boundary")
        self.assertEqual(recovery.classify("repeated step rejection: quantum ion EOS: outside assessed ionization/quantum range"), "table_boundary")

    def test_later_control_does_not_hide_primary_failure(self):
        for name, reason, age in [("main", "outside assessed ionization/quantum range", 106e9),
                                  ("control", "accepted-step budget", 792e6)]:
            p = self.work / "cases" / name / "evolution"
            p.mkdir(parents=True)
            (p.parent / "launch.json").write_text("{}")
            (p / "report.json").write_text(json.dumps(dict(stop_reason=reason)))
            (p / "history.jsonl").write_text(json.dumps(dict(years=age)) + "\n")
            time.sleep(.01)
        for name in ["main", "cases/main"]:
            (self.work / "progress.json").write_text(json.dumps(dict(current_case=name)))
            observed = recovery.observe(self.work)
            self.assertEqual(observed["case"], "cases/main")
            self.assertEqual(observed["kind"], "table_boundary")
            self.assertEqual(observed["accepted"]["years"], 106e9)
            self.assertEqual(observed["latest_case"], "cases/control")

    def test_unknown_primary_uses_known_launch(self):
        p = self.work / "cases" / "main" / "evolution"
        p.mkdir(parents=True)
        (p.parent / "launch.json").write_text("{}")
        for name in ["missing", "../../outside"]:
            (self.work / "progress.json").write_text(json.dumps(dict(current_case=name)))
            observed = recovery.observe(self.work)
            self.assertEqual(observed["case"], "cases/main")
            self.assertEqual(observed["selection"], "latest_launch")

    def test_disabled_overflow_is_not_exhausted_subscription(self):
        self.assertIsNone(supervisor.usage_stop_reason(dict(status="allowed", overageStatus="rejected"), None))
        self.assertIsNone(supervisor.provider_reset(dict(status="allowed", overageStatus="rejected"), time.time()))


if __name__ == "__main__":
    unittest.main()
