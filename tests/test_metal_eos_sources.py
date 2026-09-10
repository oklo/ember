#!/usr/bin/env python3
"""Absolute element numbers must remain valid when hydrogen is exhausted."""
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from metal_eos_composition import mixture
from plan_metal_eos_refinement import observables
from import_metal_eos import fraction_label


class ZeroHydrogen(unittest.TestCase):
    def test_refinement_response_with_gas_and_radiation(self):
        # Equal gas/radiation pressures, monatomic gas. This exercises the
        # compressibility factor in grad_ad, absent from an ideal-gas-only test.
        p, rho, temperature = 2.e12, .1, 2.e6
        scale = p/(rho*temperature)
        values = observables([p, 1., 6.75*scale, .5*p, 2.5*p], temperature, rho)
        self.assertAlmostEqual(values[3]/scale, 19.25)
        self.assertAlmostEqual(values[4], 5/19.25)

    def test_subpermille_source_names_are_unique(self):
        fractions = [0, .001, .002, .002125, .00325, .004375, .0055, .00775]
        labels = [fraction_label(v, 1000) for v in fractions]
        self.assertEqual(len(set(labels)), len(fractions))
        self.assertEqual(labels[:3], ['000', '001', '002'])
        for x, label in zip(fractions, labels):
            self.assertEqual(float(label.replace('p', '.'))/1000, x)

    def test_true_zero_hydrogen_and_isotope_number_counts(self):
        reference=mixture(.7,0)
        for he3 in [0,.004,.12]:
            state=mixture(0,he3)
            numbers=[v*state['source_mass_scale'] for v in state['eps']]
            self.assertEqual(numbers[0],0.)
            self.assertAlmostEqual(numbers[1],he3/3+(.98-he3)/4,places=15)
            for i in range(2,len(numbers)):
                self.assertAlmostEqual(numbers[i],reference['eps'][i]*reference['source_mass_scale'],places=17)
            self.assertEqual(state['potassium_number_per_baryon_mass'],reference['potassium_number_per_baryon_mass'])

    def test_continuity_as_hydrogen_goes_to_zero(self):
        zero=mixture(0,.004);small=mixture(1e-10,.004)
        for a,b in zip(zero['eps'],small['eps']):self.assertAlmostEqual(a,b,delta=2e-10)
        self.assertAlmostEqual(zero['source_mass_scale'],small['source_mass_scale'],delta=2e-10)

    def test_invalid_compositions_rejected(self):
        for x,y in [(float('nan'),0),(0,float('inf')),(-1,0),(0,-1),(0,1)]:
            with self.assertRaises(ValueError):mixture(x,y)


if __name__=='__main__':unittest.main()
