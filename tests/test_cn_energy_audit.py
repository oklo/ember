"""The independent fuel audit must remain useful when nuclear burning vanishes."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from evolve_cn_transport import nuclear_mass_energy_error


class NuclearEnergyAudit(unittest.TestCase):
    def test_burning_budget_is_retained(self):
        nuclear, photon = 1e31, 9e30
        for fractional_defect in [-3e-6, -1e-6, 0, 1e-6, 3e-6]:
            release = nuclear * (1 + fractional_defect)
            old = release / nuclear - 1
            new = nuclear_mass_energy_error(release, nuclear, photon)
            self.assertAlmostEqual(new, old, delta=2e-16)
            self.assertEqual(abs(old) <= 2e-6, abs(new) <= 2e-6)

    def test_zero_burning_and_finite_cooling(self):
        self.assertEqual(nuclear_mass_energy_error(0, 0, 1e23), 0)

    def test_cooling_still_detects_an_energy_defect(self):
        photon = 1e23
        for nuclear in [0, 1e-200, 1e12]:
            for sign in [-1, 1]:
                small = nuclear_mass_energy_error(nuclear + sign * 1e-6 * photon,
                                                 nuclear, photon)
                large = nuclear_mass_energy_error(nuclear + sign * 3e-6 * photon,
                                                 nuclear, photon)
                self.assertLess(abs(small), 2e-6)
                self.assertGreater(abs(large), 2e-6)

    def test_power_units_and_cold_luminosity_do_not_set_an_artificial_floor(self):
        for power in [1e-200, 1, 1e200]:
            self.assertAlmostEqual(nuclear_mass_energy_error(3e-6 * power, 0, power),
                                   3e-6, delta=1e-20)

    def test_invalid_and_unscaled_powers_reject(self):
        for values in [(0, 0, 0), (0, -1, 1), (0, 1, -1),
                       (float('nan'), 1, 1), (0, float('inf'), 1),
                       (0, 0, float('nan'))]:
            with self.assertRaises(ValueError):
                nuclear_mass_energy_error(*values)


if __name__ == '__main__':
    unittest.main()
