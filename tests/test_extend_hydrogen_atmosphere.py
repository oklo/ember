import itertools
import math
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from extend_hydrogen_atmosphere import append_hydrogen_rows


class HydrogenExtension(unittest.TestCase):
    def setUp(self):
        self.axes = [[.45, .7], [0., .12], [math.log10(3500), math.log10(3700)], [4.9, 5.15]]
        self.lines = ['hydrogen 2 0.45000000000000001 0.69999999999999996', 'data']
        self.lines += ['1 3.7000000000000002 6.7999999999999998', *['0']*15]
        self.state = {'T': 5200., 'Pgas': 6.4e6}

    def test_old_rows_and_new_missing_corners_retained(self):
        key = (.85, 0., 3500., 4.9)
        lines, axes, added = append_hydrogen_rows(self.lines, self.axes, [(key, self.state)])
        self.assertEqual(lines[2:18], self.lines[2:])
        self.assertEqual(axes[0], [.45, .7, .85])
        self.assertEqual(lines[0].split()[2:4], self.lines[0].split()[2:])
        self.assertEqual(len(lines)-2, math.prod(map(len, axes)))
        rows = dict(zip(itertools.product(*axes), lines[2:], strict=True))
        self.assertEqual(sum(v.startswith('1 ') for v in rows.values()), 2)
        self.assertEqual(rows[(.85, .12, math.log10(3700), 5.15)], '0')
        self.assertEqual(added, {(.85, 0., math.log10(3500), 4.9)})

    def test_existing_or_interior_hydrogen_rejected(self):
        for x in [.45, .6, .7]:
            with self.assertRaisesRegex(ValueError, 'only higher'):
                append_hydrogen_rows(self.lines, self.axes, [((x, 0., 3500., 4.9), self.state)])

    def test_other_axes_cannot_change(self):
        with self.assertRaisesRegex(ValueError, 'other source axes'):
            append_hydrogen_rows(self.lines, self.axes, [((.85, 0., 3600., 4.9), self.state)])

    def test_duplicate_sources_rejected(self):
        entry = ((.85, 0., 3500., 4.9), self.state)
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            append_hydrogen_rows(self.lines, self.axes, [entry, entry])


if __name__ == '__main__':
    unittest.main()
