#!/usr/bin/env python3
"""Control-flow tests for solve_column.py. Acceptance uses import_nongrey_grid and
is not re-tested here; the fake source never produces an acceptable column."""
import hashlib, json, os, shutil, subprocess, sys, tempfile, unittest
from pathlib import Path
from unittest import mock
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
SOLVER = ROOT / 'scripts/solve_column.py'; FAKE = HERE / 'fixtures/fake_tlusty.py'
SHA = hashlib.sha256(FAKE.read_bytes()).hexdigest()
ND = 20

def column(tmp, newton, damped, sleep=0.3):
    c = Path(tmp) / 'col'; c.mkdir()
    (c / 'fort.5').write_text('4000. 6.0\nT F\n')
    (c / 'tas').write_text(f'IOPTAB=-1,IFRSET=5000\nND={ND},NITER=200,CHMAX=1e-08\nORELAX=0.3\n')
    (c / 'ember-masses.dat').write_text('1.0\n'); (c / 'fort.15').write_text("'opacity.bin' 1\n")
    (c / 'opacity.bin').write_text('x'); (c / 'data').mkdir()
    spec = dict(hydrogen=[0.98], helium3=[0.0], metals=[0.0]*5, teff_K=[4000.0], log_g=[6.0], tau=100,
                log_temperature=[3, 3.0, 4.0], log_density=[3, -10.0, -1.0])
    (c / 'specification.json').write_text(json.dumps(spec))
    (c / 'provenance.json').write_text(json.dumps({'executables': {'tlusty': SHA}}))
    (c / 'opacity.sha256').write_text(hashlib.sha256(b'x').hexdigest()+'\n')
    (c / 'fake.json').write_text(json.dumps(dict(newton=newton, damped=damped, nd=ND, sleep=sleep)))
    return c

def solve(c, *extra, exe=str(FAKE), cwd=None):
    return subprocess.run([sys.executable, str(SOLVER), str(c), '--executable', exe, '--executable-sha256', SHA, *extra],
                          capture_output=True, text=True, cwd=cwd)

class T(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(); self.addCleanup(shutil.rmtree, self.tmp)
        os.environ['EMBER_ROOT'] = str(ROOT)
    def test_refuses_existing_attempts(self):
        c = column(self.tmp, [1e-2], [1e-2]); (c / 'attempts').mkdir(); (c / 'attempts' / 'keep').write_text('x')
        r = solve(c); self.assertNotEqual(r.returncode, 0); self.assertIn('attempts/ already exists', r.stderr + r.stdout)
        self.assertTrue((c / 'attempts' / 'keep').exists())
    def test_prepared_source_mismatch_is_rejected(self):
        c = column(self.tmp, [1], [1])
        (c / 'opacity.bin').write_text('changed')
        r = solve(c)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn('opacity differs', r.stderr+r.stdout)
        (c / 'opacity.bin').write_text('x')
        (c / 'provenance.json').write_text(json.dumps({'executables': {'tlusty': 'wrong'}}))
        r = solve(c)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn('prepared source physics', r.stderr+r.stdout)
    def test_relative_executable_is_normalized(self):
        c = column(self.tmp, [1e-2, 1e-3], [1e-2])
        rel = os.path.relpath(FAKE, HERE)
        r = solve(c, '--max-phases', '1', exe=rel, cwd=HERE)
        res = json.loads((c / 'result.json').read_text())
        self.assertEqual(res['attempts'][0]['returncode'], 0, r.stderr)   # the fake ran from inside attempts/00_*
    def test_missing_metadata_fails_before_source_launch(self):
        c = column(self.tmp, [1], [1])
        scripts = Path(self.tmp) / 'snapshot' / 'scripts'; scripts.mkdir(parents=True)
        for name in ['solve_column.py', 'generate_nongrey_grid.py', 'import_nongrey_grid.py',
                     'prepare_nongrey_sources.py', 'nongrey_opacity.py']:
            shutil.copyfile(SOLVER.parent / name, scripts / name)
        env = {k:v for k,v in os.environ.items() if k != 'EMBER_ROOT'}
        result = subprocess.run([sys.executable, str(scripts / 'solve_column.py'), str(c),
                                 '--executable', str(FAKE), '--executable-sha256', SHA],
                                capture_output=True, text=True, env=env)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('synple-elements.json', result.stderr)
        self.assertFalse((c / 'attempts').exists())
    def test_invalid_budgets_rejected(self):
        for args in (['--total-cpu', '0'], ['--attempt-cpu', 'nan'], ['--max-phases', '0'], ['--chmax', '1e-3']):
            c = column(tempfile.mkdtemp(dir=self.tmp), [1], [1])
            self.assertNotEqual(solve(c, *args).returncode, 0, args)
    def test_too_small_budget_writes_unaccepted_result(self):
        c = column(self.tmp, [1e-2], [1e-2]); r = solve(c, '--total-cpu', '10')
        res = json.loads((c / 'result.json').read_text())
        self.assertFalse(res['accepted']); self.assertIn('budget', res['failure']); self.assertEqual(res['attempts'], [])
    def test_caller_log_is_not_part_of_source_checksums(self):
        c = column(self.tmp, [1e-2], [1e-2])
        with (c / 'controller.log').open('w') as log:
            subprocess.run([sys.executable, str(SOLVER), str(c), '--executable', str(FAKE),
                            '--executable-sha256', SHA, '--total-cpu', '10'], stdout=log, check=True)
            log.write('caller finished\n')
        hashes = json.loads((c / 'result.json').read_text())['retained_sha256']
        self.assertNotIn('controller.log', hashes)
        for name in ('fort.5', 'tas', 'specification.json', 'provenance.json', 'opacity.sha256'):
            self.assertEqual(hashes[name], hashlib.sha256((c / name).read_bytes()).hexdigest())
    def test_last_phase_is_damped(self):
        grow = [1e-2, 5e-3, 1e-1, 5.0, 50.0]     # undamped growth -> monitor stop
        for phases in (1, 2, 3, 4, 5):
            c = column(tempfile.mkdtemp(dir=self.tmp), grow, [1e-1, 5e-2], sleep=0.25)
            solve(c, '--max-phases', str(phases))
            kinds = [a['kind'] for a in json.loads((c / 'result.json').read_text())['attempts']]
            self.assertEqual(kinds[-1], 'damped', (phases, kinds))
            self.assertLessEqual(len(kinds), phases)
    def test_stop_and_reap_tolerates_exited_process(self):
        sys.path.insert(0, str(SOLVER.parent)); import solve_column as sc
        p = subprocess.Popen(['/usr/bin/true'], start_new_session=True)
        import time; time.sleep(0.2)
        with mock.patch('os.killpg', side_effect=ProcessLookupError):
            pid, status, usage = sc.stop_and_reap(p.pid)
        p.returncode = os.waitstatus_to_exitcode(status)
        self.assertEqual(pid, p.pid)

if __name__ == '__main__':
    unittest.main(verbosity=2)
