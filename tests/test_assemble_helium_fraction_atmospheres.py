"""Composition-domain and missing-coverage checks for the source assembler."""
import itertools
import math
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from assemble_helium_fraction_atmospheres import physical_helium3, render_table


class FractionAssemblyTest(unittest.TestCase):
    def setUp(self):
        self.metals = [1e-16, 0, 0, 0, 0]
        self.plan = dict(hydrogen=[.87, .9955], helium3_fraction=[0, .95],
                         teff_K=[4500., 4600., 4700.], log_g=[5.55, 5.8],
                         tau=100., metal_tolerance=1e-15, source='analytic test only')
        self.states = {key: dict(T=key[2]*1.2, Pgas=10**key[3]*100)
                       for key in itertools.product(*(self.plan[k] for k in
                           ['hydrogen', 'helium3_fraction', 'teff_K', 'log_g']))}

    def test_enclosing_compositions_remain_physical(self):
        # A rectangular absolute-He3 grid with these extremes would be unphysical.
        for h, f in itertools.product(self.plan['hydrogen'], self.plan['helium3_fraction']):
            x3 = physical_helium3(h, f, self.metals)
            x4 = 1-h-sum(self.metals)-x3
            self.assertGreaterEqual(x4, 0)
            self.assertAlmostEqual(x3/(x3+x4), f)
        self.assertGreater(.9955+physical_helium3(.87,.95,self.metals), 1)

    def test_missing_corner_removes_only_incident_cells(self):
        del self.states[(.87, 0, 4700., 5.55)]
        text, cells = render_table(self.plan, self.states, self.metals, 'test')
        self.assertEqual(len(cells), 1)
        self.assertEqual(cells[0][2], [4500., 4600.])
        rows = text.split('data\n')[1].splitlines()
        self.assertEqual(len(rows), 24)
        self.assertEqual(rows.count('0'), 1)
        self.assertEqual(rows[4], '0')
        # Independent expected ordering: H, fraction, Teff, gravity, with gravity fastest.
        t, pg = map(float, rows[3].split()[1:])
        self.assertAlmostEqual(10**t, 4600*1.2)
        self.assertAlmostEqual(math.log10(10**pg/100), 5.8)

    def test_zero_weight_missing_nodes_cannot_supply_a_cell(self):
        for key in list(self.states):
            if key[1] == .95:
                del self.states[key]
        with self.assertRaisesRegex(ValueError, 'no complete interpolation cell'):
            render_table(self.plan, self.states, self.metals, 'test')

    def test_invalid_values_and_coordinates_rejected(self):
        for h, f in [(1., 0), (.9, -1e-4), (.9, 1.001)]:
            with self.assertRaises(ValueError):
                physical_helium3(h, f, self.metals)
        wrong = dict(self.states)
        wrong[(.87, .1235, 4500., 5.55)] = dict(T=5000, Pgas=1e7)
        with self.assertRaisesRegex(ValueError, 'unexpected source coordinate'):
            render_table(self.plan, wrong, self.metals, 'test')
        for value in [-1, float('nan'), float('inf')]:
            bad = dict(self.states)
            bad[next(iter(bad))] = dict(T=value, Pgas=1e7)
            with self.assertRaisesRegex(ValueError, 'invalid physical matching state'):
                render_table(self.plan, bad, self.metals, 'test')


if __name__ == '__main__':
    unittest.main()
