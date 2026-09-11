"""Analytic controls for the independent spectral-averaging diagnostic.

Run with a Python environment containing NumPy and SciPy.
"""
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_tops_spectral_means import Integrals


class SpectralMeanTests(unittest.TestCase):
    def test_grey_absorption_and_scattering_at_any_cutoff(self):
        u = np.geomspace(1e-7, 800., 10001)
        data = np.column_stack([u, np.full_like(u, 7.), np.full_like(u, 5.), np.full_like(u, 2.)])
        integral = Integrals(data, 1., 16)
        for cutoff in [.1, 1.3, 40.]:
            value = integral.above(cutoff)
            self.assertAlmostEqual(value[0]/value[1], 7., places=12)
            self.assertAlmostEqual(value[3]/value[2], 5., places=12)
            full = (4*np.pi**4/15)/value[1]
            self.assertAlmostEqual(full*value[0]/(4*np.pi**4/15), 7., places=12)
            self.assertGreater((4*np.pi**4/15)/integral.refractive(cutoff)[1], full)

    def test_power_law_against_planck_integrals(self):
        u = np.geomspace(1e-7, 800., 10001)
        data = np.column_stack([u, 3*u**2, 2*u**2, u**2])
        value = Integrals(data, 1., 16).tail[0]
        rosseland = (4*np.pi**4/15)/value[1]
        planck = value[3]/(np.pi**4/15)
        self.assertLess(abs(rosseland/(3*.8*np.pi**2)-1), 4e-8)
        self.assertLess(abs(planck/(2*40*np.pi**2/21)-1), 1e-12)


if __name__ == '__main__':
    unittest.main()
