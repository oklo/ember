"""Independent limits and input checks for the scattering-inclusive mean."""
import sys
import unittest
import tempfile
from pathlib import Path

import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from gas_rosseland import isotherm, population_export, rosseland, write_opacity_table


class GasRosselandTests(unittest.TestCase):
    def test_grey_and_harmonic_limits(self):
        nu = np.geomspace(1e10, 1e17, 20000)
        np.testing.assert_allclose(rosseland(nu, 3000, np.full((len(nu), 2), [2., 7.])), [2., 7.])
        opacity = np.where(nu < 1e14, .01, 100.)
        self.assertLess(rosseland(nu, 3000, opacity), 1.)
        self.assertGreater(rosseland(nu, 3000, opacity), .01)

    def test_particle_count_is_not_atomic_hydrogen(self):
        nu = np.geomspace(1e12, 1e16, 1000)
        table = dict(shape=(len(nu),2,1), frequency=nu[::-1],
                     log_temperature=[np.log(2000.)], log_density=np.log([.1,1.]),
                     log_electron_density=np.log([1.,1.]), log_opacity=np.log(np.full(len(nu)*2,.1)))
        p = np.array([[2000,.1,1,0,1e20,1e18,1e23], [2000,1,1,0,1e21,1e19,1e24]])
        a, _ = isotherm(table,p)
        p[:,6] *= 100
        b, _ = isotherm(table,p)
        np.testing.assert_array_equal(a,b)
        self.assertTrue(np.all(a > .1))
        p[0,1] *= 2
        with self.assertRaises(ValueError): isotherm(table,p)

    def test_invalid_spectrum(self):
        with self.assertRaises(ValueError): rosseland([2.,1.],3000,[1.,1.])
        with self.assertRaises(ValueError): rosseland([1.,2.],3000,[1.,0.])

    def test_density_support_is_explicit(self):
        rows = {(1., t, r): 10**r for t in [100.,200.,300.,400.]
                for r in range(4 if t < 300 else 5)}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'table.dat'
            write_opacity_table(rows, path)
            lines = path.read_text().splitlines()
            self.assertEqual(lines[:2], ['EMBER_OPACITY_TABLE 2',
                                        '1 4 5 HHe_absorption_electron_Rayleigh'])
            self.assertEqual([int(line.split()[0]) for line in lines[5:]], [4,4,5,5])
            del rows[1., 300., 2]
            with self.assertRaises(ValueError): write_opacity_table(rows, path)
            del rows[1., 100., 0]
            with self.assertRaises(ValueError): write_opacity_table(rows, path)

    def test_rectangular_format_is_retained(self):
        rows = {(1., t, r): 10**r for t in [100.,200.,300.,400.] for r in range(4)}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'table.dat'
            write_opacity_table(rows, path)
            lines = path.read_text().splitlines()
            self.assertEqual(lines[0], '1 4 4 HHe_absorption_electron_Rayleigh')
            self.assertEqual(lines[4:], ['0 1 2 3']*4)

    def test_export_is_diagnostic_and_idempotent(self):
        source = '      subroutine ougrid(abso)\n      if (nfreq.le.3) return \n      end\n'
        changed = population_export(source)
        self.assertEqual(population_export(changed), changed)
        self.assertTrue(changed.endswith('      end\n'))
        self.assertEqual(changed.count('write(184,*)'),1)


if __name__ == '__main__': unittest.main()
