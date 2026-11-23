import json
from pathlib import Path
import sys
import tempfile
import time
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from ember_campaign_recovery import save, tick

SOLVER = """
import json
from pathlib import Path
spec = json.loads(Path('spec.json').read_text())
guess = json.loads(Path(spec['resume_guess']).read_text())
assert spec['teff'] == 2850 and guess['temperature_K'][0] == 110
assert spec['resume_region_top'] == 3
Path('hhe_arm_a_t.json').write_text(json.dumps({'spec':spec,'quality':{'usable':True}}))
"""


class RecoveryTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.column = self.root / 'original'
        self.column.mkdir()
        (self.column / 'solve.py').write_text(SOLVER)
        save(self.column / 'spec.json', dict(anchor='a', tag='_t', levels=6, teff=2850))
        initial = dict(pressure_bar=list(range(1, 7)), temperature_K=[100]*6, rcb_guess=3)
        save(self.column / 'start_a_t.json', initial)
        self.guess = dict(initial, temperature_K=[110]*6, accepted=False, outer_calls=2)
        save(self.column / 'latest_guess_a_t.json', self.guess)
        self.job = dict(name='column', trigger='blocked', outputs=[], wall_seconds=5,
                        resume_column=dict(directory=str(self.column), max_restarts=1,
                                           cpu_seconds=10, python=sys.executable))
        save(self.root / 'recovery-plan.json', [self.job])
        d = self.root / 'recovery' / 'column'
        d.mkdir(parents=True)
        save(d / 'result.json', dict(status='failed', exit_code=-24))

    def status(self, slots=1, allow=True):
        return tick(self.root, {9: ('blocked', None)}, slots, allow)[1]['column']

    def test_budget_interrupt(self):
        self.assertEqual(self.status(allow=False), 'ready')
        self.assertEqual(self.status(slots=0), 'ready')
        self.assertEqual(self.status(), 'running')
        deadline = time.monotonic() + 5
        while self.status() == 'running' and time.monotonic() < deadline:
            time.sleep(.03)
        self.assertEqual(self.status(), 'passed')
        winner = self.root / 'recovery/column/selected-column.json'
        self.assertTrue(Path(json.loads(winner.read_text())['receipt']).is_file())

    def test_refuse_unchanged_or_invalid_guess(self):
        for temps in ([100]*6, [float('nan')]*6, [110]):
            save(self.column / 'latest_guess_a_t.json', dict(self.guess, temperature_K=temps))
            self.assertEqual(self.status(), 'failed')

    def test_scientific_failure(self):
        save(self.column / 'hhe_arm_a_t.json', {'quality': {'usable': False}})
        self.assertEqual(self.status(), 'failed')

    def test_retry_limit(self):
        child = self.root / 'recovery/column-resume-01'
        child.mkdir()
        save(child / 'input.json', dict(name=child.name, column_directory=str(self.column),
                                        seed=dict(profile_sha256='previous')))
        save(child / 'result.json', dict(status='failed', exit_code=-24))
        self.assertEqual(self.status(), 'failed')


if __name__ == '__main__':
    unittest.main()
