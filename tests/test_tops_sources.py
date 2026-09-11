#!/usr/bin/env python3
"""Regression checks for refinement archives and selected input provenance."""
import json
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from fetch_tops_composition import fraction_label, cached_result
from import_tops_composition import validate_mixture, read as read_tops
from run_evolution_snapshot import input_data
from audit_tops_heldout import active_top_weight, KEV_TO_K


class SourceTests(unittest.TestCase):
    def test_zero_hydrogen_response_omits_h_but_retains_every_metal(self):
        # Normalized composition from the actual TOPS zero-H response. This
        # excerpt tests composition identity, not missing opacity cell data.
        root=Path(__file__).resolve().parents[1]
        text=(root/'tests/fixtures/tops_zero_hydrogen_composition.txt').read_text()
        request=json.loads((root/'data/opacity/sources/tops_gs98_x070_z020.request.json').read_text())
        tokens=request['mixture'].split(' he ',1)[1].split()
        metals={tokens[i+1].capitalize():float(tokens[i]) for i in range(0,len(tokens),2)}
        expected={'X':0.,'Z':.02,'metals':metals}
        elements=validate_mixture(text,expected)
        self.assertNotIn('H',elements)
        self.assertEqual(set(elements),set(metals)|{'He'})
        for x in [1e-8,.025]:
            with self.assertRaisesRegex(ValueError,'wrong source mixture'):
                validate_mixture(text,{**expected,'X':x})
        with self.assertRaisesRegex(ValueError,'wrong source mixture'):
            validate_mixture(text.replace('9.8000E-01','9.7000E-01'),expected)
        with self.assertRaisesRegex(ValueError,'wrong source metal abundance'):
            validate_mixture(text.replace('3.4367E-03','3.0000E-03'),expected)
        missing='\n'.join(row for row in text.splitlines() if ' Fe ' not in row).replace('materials =  20','materials =  19')
        with self.assertRaisesRegex(ValueError,'wrong source mixture'):
            validate_mixture(missing,expected)
        with self.assertRaisesRegex(ValueError,'wrong source mixture'):
            validate_mixture(text.replace('materials =  20','materials =  21'),expected)

    def test_density_dimensions_require_an_explicit_request(self):
        root=Path(__file__).resolve().parents[1]
        text=(root/'tests/fixtures/tops_zero_hydrogen_composition.txt').read_text()
        request=json.loads((root/'data/opacity/sources/tops_gs98_x070_z020.request.json').read_text())
        tokens=request['mixture'].split(' he ',1)[1].split()
        expected={'X':0.,'Z':.02,'metals':{
            tokens[i+1].capitalize():float(tokens[i]) for i in range(0,len(tokens),2)}}
        import re
        changed=re.sub(r'(Number of rho\s*=\s*)71',r'\g<1>11',text)
        self.assertNotEqual(changed,text)
        with self.assertRaisesRegex(ValueError,'unexpected dimensions'):
            validate_mixture(changed,expected)
        self.assertEqual(validate_mixture(changed,expected,dimensions=(50,11)),
                         validate_mixture(text,expected))
        with self.assertRaisesRegex(ValueError,'unexpected dimensions'):
            validate_mixture(changed,expected,dimensions=(50,12))
        with self.assertRaisesRegex(ValueError,'wrong source metal abundance'):
            validate_mixture(changed.replace('3.4367E-03','3.0000E-03'),expected,
                             dimensions=(50,11))

    def test_supported_isotherm_without_a_warning_block(self):
        import re
        root=Path(__file__).resolve().parents[1]
        text=(root/'tests/fixtures/tops_zero_hydrogen_composition.txt').read_text()
        text=re.sub(r'(Number of T\s*=\s*)50',r'\g<1>1',text)
        text=re.sub(r'(Number of rho\s*=\s*)71',r'\g<1>2',text)
        # Synthetic values isolate the response format from opacity physics.
        text+='Density     Ross opa    Planck opa  No. Free    Av Sq Free  T=  .225\n1.E5 2.E1 3.E1 1. 1.\n2.E5 3.E1 4.E1 1. 1.\n'
        request=json.loads((root/'data/opacity/sources/tops_gs98_x070_z020.request.json').read_text())
        tokens=request['mixture'].split(' he ',1)[1].split()
        expected={'X':0.,'Z':.02,'metals':{
            tokens[i+1].capitalize():float(tokens[i]) for i in range(0,len(tokens),2)}}
        with tempfile.TemporaryDirectory() as work:
            path=Path(work)/'source.txt';path.write_text(text)
            expected['sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
            tt,rr,cells,excluded=read_tops(path,expected,dimensions=(1,2))
            self.assertEqual((tt,rr),([.225],[1e5,2e5]))
            self.assertEqual(cells[.225,2e5],30.)
            self.assertFalse(excluded)
            path.write_text(text.replace('Normalized composition','Unknown density substitution warning\nNormalized composition'))
            expected['sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
            with self.assertRaisesRegex(ValueError,'unrecognized source density warning'):
                read_tops(path,expected,dimensions=(1,2))

    def test_active_domain_includes_cool_tops_and_respects_both_rectangles(self):
        rectangles={'low':(1e-10,251.19),'high':(1e-10,1e4)}
        self.assertIsNone(active_top_weight(.002,1,rectangles,[-8,6]))
        self.assertEqual(active_top_weight(.003,1,rectangles,[-8,6]),1)
        self.assertIsNone(active_top_weight(.025,1000,rectangles,[-8,6]))
        self.assertEqual(active_top_weight(.05,1000,rectangles,[-8,6]),1)
        self.assertAlmostEqual(active_top_weight(10**4.45/KEV_TO_K,.01,rectangles,[-8,6]),.5,places=13)
        self.assertIsNone(active_top_weight(10**4.45/KEV_TO_K,100,rectangles,[-8,6]))

    def test_refinement_names_do_not_alias(self):
        samples=[.1+i*.000125 for i in range(401)]
        # Decimal inputs represent the exact compositions in the request.
        labels=[fraction_label(f'{x:.6f}',100) for x in samples]
        self.assertEqual(len(labels),len(set(labels)))
        self.assertEqual(fraction_label(.1,100),'010')
        self.assertEqual(fraction_label(.125,100),'012p5')
        self.assertEqual(fraction_label(.1125,100),'011p25')
        for x in [float('nan'),float('inf'),-.01,1e-11]:
            with self.assertRaises(ValueError):fraction_label(x,100)

    def test_partial_archive_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);source=d/'sample.txt';source.write_bytes(b'original source response')
            with self.assertRaisesRegex(ValueError,'incomplete existing'):
                cached_result(d,'sample',.125,.02,{})
            self.assertEqual(source.read_bytes(),b'original source response')

    def test_snapshot_follows_external_source_references(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);selected=root/'selected';selected.mkdir()
            plane=root/'external plane.dat';plane.write_text('fixture')
            for name in ['aesopus21_gs98_mixture.dat','tops_gs98_mixture_low.dat','tops_gs98_mixture_high.dat']:
                (selected/name).write_text('EMBER_OPACITY_MIXTURE 1 2 logRho fixture\n'
                                           '0.01 "../external plane.dat"\n0.02 "../external plane.dat"\n')
            found=input_data(['--opacity-directory',str(selected)])
            self.assertIn(plane.resolve(),found)
            self.assertEqual(sum(p==plane.resolve() for p in found),1)


if __name__=='__main__':unittest.main()
