import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts/solar_grid'))
from cooling_interval import compare


class CoolingIntervalTest(unittest.TestCase):
    def test_clock_uses_common_temperature_and_elapsed_time(self):
        # The restored photosphere can be a few microkelvin warmer than the forecast anchor.
        native = [dict(Teff_K=t, years=age, central_T_K=1e6, radius_Rsun=.1, nuclear_fraction=.002)
                  for t,age in [(900.00001,1e10), (890,1.1e10), (880,1.2e10)]]
        reduced = [dict(Teff_K=t, age_yr=age, Tc_K=1e6, R_ph_cm=6.957e9, Lnuc_over_L=.002)
                   for t,age in [(900,3e10), (890,3.13e10), (880,3.26e10)]]
        result = compare(native,reduced,900,890)
        self.assertAlmostEqual(result['errors']['elapsed_relative'],.30,places=4)
        self.assertEqual(result['errors']['core_relative'],0)

    def test_insufficient_resolution_cannot_qualify_wide_spacing(self):
        native = [dict(Teff_K=t, years=y, central_T_K=1, radius_Rsun=1, nuclear_fraction=0)
                  for t,y in [(850,1e10),(800,2e10)]]
        reduced = [dict(Teff_K=r['Teff_K'],age_yr=r['years'],Tc_K=1,R_ph_cm=6.957e10,Lnuc_over_L=0)
                   for r in native]
        result = compare(native,reduced,850,825)
        self.assertEqual(result['errors']['elapsed_relative'],0)
        self.assertIsNone(result['elapsed_error_bracket'])
