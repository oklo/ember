"""Independent normalization checks for highly degenerate source integrals."""
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.special import expit

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from electron_electron_energy_modes_v2 import PolynomialBasis, _three_panel_rule
from electron_electron_collision import energy_basis


class PairSourceTests(unittest.TestCase):
    def test_basis_at_extended_degeneracy(self):
        # Integrate independently with panelled Legendre rules, rather than
        # the Jacobi rule used to construct the orthogonal polynomials.
        for eta in [64.,96.,128.]:
            basis=PolynomialBasis(eta,9)
            x,w=_three_panel_rule(768,[0.,eta-12.,eta+12.,eta+80.])
            q=basis.values(x)
            q[1]*=np.sqrt(basis.base.variance)
            weight=w*x**1.5*expit(eta-x)*expit(x-eta)/basis.base.normalization
            self.assertLess(np.max(abs((q*weight)@q.T-np.eye(10))),2e-8)

    def test_finer_normalization_is_needed_at_high_degeneracy(self):
        with self.assertRaisesRegex(ValueError,'orthogonality'):
            PolynomialBasis(128.,9,points=384)
        self.assertLess(PolynomialBasis(128.,9).gram_error,2e-9)
        self.assertEqual(PolynomialBasis(20.,9).points,384)

    def test_source_domain_and_composite_quadrature(self):
        for eta in [float('nan'),-31.,129.]:
            with self.assertRaises(ValueError):energy_basis(eta)
        x,w=_three_panel_rule(40,[-1.,-.02,.02,1.])
        self.assertTrue(np.all(w>0))
        for power in range(6):
            exact=0. if power%2 else 2/(power+1)
            self.assertAlmostEqual(float(w@x**power),exact,places=13)


if __name__=='__main__':unittest.main()
