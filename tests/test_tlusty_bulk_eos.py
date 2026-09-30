"""Exercise the Fortran bulk-EOS reader with both table formats."""
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


@unittest.skipUnless(shutil.which('gfortran'), 'gfortran is required')
class NativeBulkEos(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.work.cleanup)
        cls.directory = Path(cls.work.name)
        driver = cls.directory/'probe.f90'
        driver.write_text('''program probe
use ember_atmosphere_bulk
implicit none
real(dp)::t,p,r,s,e
read(*,*)t,p
call bulk_value(t,p,r,s,e)
write(*,'(3es26.17)')r,s,e
end program
''')
        source = Path(__file__).resolve().parents[1]/'scripts/tlusty_bulk_eos.f90'
        cls.exe = cls.directory/'probe'
        subprocess.run(['gfortran', '-O2', '-fcheck=all', str(source), str(driver),
                        '-o', str(cls.exe)], cwd=cls.directory, check=True,
                       capture_output=True, text=True)

    def run_query(self, version, lt, lp, mask='0 1', pressures=(1, 2)):
        # Exact log-linear fields: density=P/T, entropy=T, energy=T**2.
        lines = [f'EMBER_ATMOSPHERE_BULK_EOS {version}',
                 f'3 {len(pressures)} .98', '2 3 4',
                 ' '.join(map(str, pressures))]
        for p in pressures:
            for t in [2, 3, 4]:
                lines.append(f'{p-t} -1 1 0 {t} 1 0 0 {2*t} 2 0 0')
        if version == 2:
            lines += ['cells', mask]
        table = self.directory/'bulk.dat'
        table.write_text('\n'.join(lines)+'\n')
        return subprocess.run([str(self.exe)], input=f'{math.exp(lt)} {math.exp(lp)}\n',
                              env=dict(os.environ, EMBER_ATMOSPHERE_BULK_EOS=str(table)),
                              capture_output=True, text=True)

    def test_warm_pressure_domain_and_interpolation_are_preserved(self):
        for version in [1, 2]:
            result = self.run_query(version, 3.5, 1.8)
            self.assertEqual(result.returncode, 0, result.stderr)
            got = list(map(float, result.stdout.splitlines()[-1].split()))
            expected = [math.exp(1.8-3.5), math.exp(3.5), math.exp(7)]
            for a, b in zip(got, expected):
                self.assertAlmostEqual(a/b, 1, places=13)

    def test_unsupported_cell_is_rejected(self):
        result = self.run_query(2, 2.5, 1.5)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('masked cell', result.stdout)
        self.assertEqual(self.run_query(1, 2.5, 1.5).returncode, 0)

    def test_shared_corner_can_use_either_supported_neighbor(self):
        for mask in ['0 1 0 0', '0 0 1 0', '0 0 0 1']:
            result = self.run_query(2, 3, 2, mask=mask, pressures=(1, 2, 3))
            self.assertEqual(result.returncode, 0, result.stderr)
            got = list(map(float, result.stdout.splitlines()[-1].split()))
            for a, b in zip(got, [math.exp(-1), math.exp(3), math.exp(6)]):
                self.assertAlmostEqual(a/b, 1, places=13)

    def test_invalid_mask_is_rejected(self):
        result = self.run_query(2, 3.5, 1.5, mask='0 2')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('invalid cell mask', result.stderr)

    def test_outside_axes_is_rejected(self):
        result = self.run_query(2, 1.5, 1.5)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('outside table', result.stdout)


if __name__ == '__main__':
    unittest.main()
