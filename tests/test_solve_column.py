#!/usr/bin/env python3
"""Control-flow tests for solve_column.py. Acceptance uses import_nongrey_grid and
is not re-tested here; the fake source never produces an acceptable column."""
import hashlib, json, math, os, shutil, subprocess, sys, tempfile, unittest
from pathlib import Path
from unittest import mock
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
SOLVER = ROOT / 'scripts/solve_column.py'; FAKE = HERE / 'fixtures/fake_tlusty.py'
SHA = hashlib.sha256(FAKE.read_bytes()).hexdigest()
ND = 20
PROV = {'executables': {'tlusty': SHA, 'synspec': 's'}, 'data_sha256': 'd', 'line_list_sha256': ['l'],
        'opacity_method': 'sampling'}

def structure(nd, teff):
    mass = [1e-6 * 10 ** (i / 4) for i in range(nd)]
    rows = [f'{teff * (1 + i / nd):.6e} {1e14:.6e} {1e-8 * (1 + i):.6e} {1e17:.6e}' for i in range(nd)]
    return f'{nd} -4\n' + '\n'.join(f'{m:.6e}' for m in mass) + '\n' + '\n'.join(rows) + '\n'

def donor(tmp, accepted=True, alpha=1.9, teff=3900.0, logg=6.0):
    d = Path(tmp) / 'donor'; d.mkdir()
    spec = dict(hydrogen=[0.98], helium3=[0.0], metals=[0.0]*5, teff_K=[teff], log_g=[logg], tau=100, alpha=alpha,
                wavelength_A=[900, 300000], microturbulence_km_s=1.0, line_threshold=1e-4)
    (d / 'specification.json').write_text(json.dumps(spec)); (d / 'provenance.json').write_text(json.dumps(PROV))
    (d / 'fort.7').write_text(structure(24, teff))
    sums = {n: hashlib.sha256((d / n).read_bytes()).hexdigest() for n in ('fort.7', 'specification.json')}
    (d / 'result.json').write_text(json.dumps(dict(accepted=accepted, retained_sha256=sums)))
    return d

def column(tmp, newton, damped, sleep=0.3):
    c = Path(tmp) / 'col'; c.mkdir()
    (c / 'fort.5').write_text('4000. 6.0\nT F\n')
    (c / 'tas').write_text(f'IOPTAB=-1,IFRSET=5000\nND={ND},NITER=200,CHMAX=1e-08\nORELAX=0.3\n')
    (c / 'ember-masses.dat').write_text('1.0\n'); (c / 'fort.15').write_text("'opacity.bin' 1\n")
    (c / 'opacity.bin').write_text('x'); (c / 'data').mkdir()
    spec = dict(hydrogen=[0.98], helium3=[0.0], metals=[0.0]*5, teff_K=[4000.0], log_g=[6.0], tau=100,
                log_temperature=[3, 3.0, 4.0], log_density=[3, -10.0, -1.0], alpha=1.9, wavelength_A=[900, 300000],
                microturbulence_km_s=1.0, line_threshold=1e-4)
    (c / 'specification.json').write_text(json.dumps(spec))
    (c / 'provenance.json').write_text(json.dumps(PROV))
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
        self.assertEqual(res['solver_sha256'], hashlib.sha256(SOLVER.read_bytes()).hexdigest())
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

    def test_cpu_limited_damping_reuses_guess_within_budget_and_still_checks_source(self):
        sys.path.insert(0, str(SOLVER.parent)); import solve_column as sc
        for budget, count in [(100, 4), (40, 2)]:
            c = column(tempfile.mkdtemp(dir=self.tmp), [1], [1])
            calls = []
            def fake_run(att, exe, cap, nd, rule):
                i = len(calls)
                if i == 2:
                    self.assertEqual((att/'fort.8').read_text(), structure(ND, 4001.))
                calls.append(cap)
                (att/'best.fort.7').write_text(structure(ND, 4000.+i))
                (att/'run.log').write_text('incomplete source output\n')
                (att/'fort.9').write_text('')
                return dict(returncode=[-sc.signal.SIGKILL, -sc.signal.SIGXCPU, 0, 0][i],
                            CPU_seconds=10 if i == 0 else 30, wall_seconds=0,
                            monitor_stop='undamped growth' if i == 0 else None,
                            iterations=5, history=[.15, .09, .006, .016, .006])
            args = [str(SOLVER), str(c), '--executable', str(FAKE),
                    '--executable-sha256', SHA, '--attempt-cpu', '60',
                    '--total-cpu', str(budget), '--max-phases', '4']
            with mock.patch.object(sys, 'argv', args), mock.patch.object(sc, 'run', fake_run), \
                    mock.patch.object(sc, 'source_inputs', wraps=sc.source_inputs) as check:
                sc.main()
            result = json.loads((c/'result.json').read_text())
            self.assertEqual(len(calls), count)
            self.assertFalse(result['accepted'])
            self.assertLessEqual(result['CPU_seconds'], budget)
            self.assertEqual(check.call_count, max(0, count-2))
            self.assertIn('next_guess', result['attempts'][1])

    def test_damped_recovery_rejects_stalls_crashes_and_incomplete_guesses(self):
        sys.path.insert(0, str(SOLVER.parent)); import solve_column as sc
        att = Path(self.tmp)
        saved = att/'best.fort.7'; saved.write_text(structure(ND, 4000.))
        base = dict(returncode=-sc.signal.SIGXCPU, history=[.15, .09, .006, .016, .006])
        self.assertTrue(sc.recover_damped_budget(base, att, ND))
        self.assertFalse(sc.recover_damped_budget(dict(base, returncode=-sc.signal.SIGSEGV), att, ND))
        self.assertFalse(sc.recover_damped_budget(dict(base, history=[.008]*6), att, ND))
        saved.write_text('20 -4\n1.0\n')
        self.assertFalse(sc.recover_damped_budget(base, att, ND))


    def test_donor_must_be_accepted(self):
        c = column(self.tmp, [1e-2], [1e-2]); d = donor(tempfile.mkdtemp(dir=self.tmp), accepted=False)
        r = solve(c, '--initial-from', str(d)); self.assertNotEqual(r.returncode, 0); self.assertIn('not accepted', r.stderr + r.stdout)
        self.assertFalse((c / 'attempts').exists())
    def test_donor_physics_must_match(self):
        c = column(self.tmp, [1e-2], [1e-2]); d = donor(tempfile.mkdtemp(dir=self.tmp), alpha=1.5)
        r = solve(c, '--initial-from', str(d)); self.assertNotEqual(r.returncode, 0); self.assertIn('physics differs', r.stderr + r.stdout)
    def test_changed_donor_output_rejected(self):
        c = column(self.tmp, [1e-2], [1e-2]); d = donor(tempfile.mkdtemp(dir=self.tmp))
        (d / 'fort.7').write_text(structure(24, 3950.0))
        r = solve(c, '--initial-from', str(d)); self.assertNotEqual(r.returncode, 0); self.assertIn('donor output changed', r.stderr + r.stdout)
    def test_donor_guess_is_scaled_and_recorded(self):
        sys.path.insert(0, str(ROOT / 'scripts'))
        from generate_nongrey_grid import continuation_structure, resample_initial_structure
        c = column(self.tmp, [1e-2, 1e-3], [1e-2]); d = donor(tempfile.mkdtemp(dir=self.tmp))
        solve(c, '--initial-from', str(d), '--max-phases', '1')
        expected = resample_initial_structure(continuation_structure((d / 'fort.7').read_text(), 3900.0, 6.0, 4000.0, 6.0), ND, allow_coarsen=True)
        attempt = next((c / 'attempts').iterdir())
        self.assertEqual((attempt / 'fort.8').read_text(), expected)
        res = json.loads((c / 'result.json').read_text())
        self.assertEqual((res['initial_from']['teff_K'], res['initial_from']['log_g']), (3900.0, 6.0))

    def test_upper_guess_keeps_mesh_and_requires_full_flux_failure_pattern(self):
        sys.path.insert(0, str(SOLVER.parent)); import solve_column as sc
        mass = [math.exp(i/5)*1e-5 for i in range(ND)]
        state = [[2700., 1e6, 1e-8, 1e16] for _ in mass]
        state[3][0] = 2900.
        def inputs(tau=1e-7, other_bad=False):
            rows = [[i+1, m, tau*(i+1), v[0], v[1], v[2], 3000., 0., 1.,
                     2. if i == 3 else 0., 3. if i == 3 else 1.]
                    for i, (m,v) in enumerate(zip(mass,state))]
            if other_bad: rows[10][10] = 1.1
            text = f'{ND} -4\n'+'\n'.join(map(str,mass))+'\n'+'\n'.join(' '.join(map(str,r))for r in state)+'\n'
            log = 'FINAL MODEL ATMOSPHERE\n'+'\n'.join(' '.join([str(int(r[0])),*map(str,r[1:])])for r in rows)
            return text, log
        text, log = inputs()
        guess, depths = sc.upper_layer_guess(text, log, ND)
        self.assertEqual(depths, [4])
        before, after = list(map(float,text.split()[2:])), list(map(float,guess.split()[2:]))
        self.assertEqual(before[:ND], after[:ND])
        self.assertAlmostEqual(after[ND+4*3], 2700., places=9)
        self.assertEqual(before[ND+4*4:], after[ND+4*4:])
        self.assertIsNone(sc.upper_layer_guess(*inputs(tau=.1), ND))
        self.assertIsNone(sc.upper_layer_guess(*inputs(other_bad=True), ND))
        self.assertIsNone(sc.upper_layer_guess(text, log.replace('2900.0', '2950.0'), ND))
        self.assertIsNone(sc.upper_layer_guess(text, log.rsplit('\n', 1)[0], ND))

if __name__ == '__main__':
    unittest.main(verbosity=2)
