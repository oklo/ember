#!/usr/bin/env python3
"""Physical acceptance checks using the independently archived cold model.

Run with the offline source-analysis Python environment (numpy required).
These alter physical labels and inputs, never the archived original files.
"""
import gzip
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from validate_condensate_model import validate
from audit_condensates import closure
from compare_condensate_control import gas_state
from audit_condensate_material import thermal_sensitivity
from archive_condensate_model import archive_model


class CondensateAcceptance(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.archive=ROOT/'data/atmosphere/sources/condensation_atmosphere_x700_2600_g515'
        cls.provenance=json.loads(gzip.decompress((cls.archive/'provenance.json.gz').read_bytes()))
        cls.spec=cls.provenance['specification'];cls.source=cls.provenance['prepared']

    def test_original_archive_replays_all_physical_checks(self):
        result=validate(self.archive,self.spec,self.source,[.7,0,2600,5.15])
        self.assertLess(result['diagnostics']['flux_error'],.002)
        self.assertTrue(result['chemistry']['grain_enthalpy_supported'])

    def test_rejects_different_family_physics(self):
        for key,value in [('metals',[.004,.0001,.001,.0093,.0056]),('opacity_frequencies',15000),('alpha',2.)]:
            with self.subTest(key=key),self.assertRaisesRegex(ValueError,'physical inputs'):
                validate(self.archive,{**self.spec,key:value},self.source,[.7,0,2600,5.15])

    def test_rejects_different_composition_or_opacity(self):
        with self.assertRaisesRegex(ValueError,'coordinate'):
            validate(self.archive,self.spec,self.source,[.45,0,2600,5.15])
        with self.assertRaisesRegex(ValueError,'opacity plane'):
            validate(self.archive,self.spec,self.source,[.7,0,2600,5.15],'0'*64)

    def test_rejects_gas_control_relabelled_as_depleted(self):
        with tempfile.TemporaryDirectory() as temp:
            target=Path(temp)/'model'
            shutil.copytree(self.archive,target,ignore=shutil.ignore_patterns('gas-control'))
            file=target/'ember-condensates.cfg.gz'
            content=gzip.decompress(file.read_bytes()).replace(b'equilibrium\n',b'gas\n')
            file.write_bytes(gzip.compress(content,mtime=0))
            with self.assertRaisesRegex(ValueError,'input fingerprint'):
                validate(target,self.spec,self.source,[.7,0,2600,5.15])

    def test_chemistry_must_belong_to_actual_profile(self):
        chemistry=json.loads(gzip.decompress((self.archive/'chemistry-audit-source.json.gz').read_bytes()))
        log=gzip.decompress((self.archive/'run.log.gz').read_bytes()).decode()
        rows=[]
        for line in log.rsplit('FINAL MODEL ATMOSPHERE',1)[1].splitlines():
            tokens=line.replace('D','E').split()
            if len(tokens)==11 and tokens[0].isdigit():rows.append(list(map(float,tokens)))
        rows.reverse();rows[0][3]*=1.01
        with self.assertRaisesRegex(ValueError,'does not match atmosphere'):
            closure(chemistry,rows,{}, {})

    def test_rejects_nonfinite_or_nonpositive_nuclei_density(self):
        chemistry=json.loads(gzip.decompress((self.archive/'chemistry-audit-source.json.gz').read_bytes()))
        # Reconstruct the source T/P sequence; other profile fields do not
        # enter this guard. Invalid conserved density must fail before any
        # relative-error comparisons can hide a NaN.
        rows=[[0,0,0,r['T_K'],0,0,r['P_bar']*1e6] for r in chemistry['rows']]
        for density in [float('nan'),float('inf'),0.,-1.]:
            chemistry['rows'][0]['total_element_density']=density
            with self.subTest(density=density),self.assertRaisesRegex(ValueError,'invalid conserved nuclei'):
                closure(chemistry,rows,{}, {})

    def test_matched_control_replays_from_archived_originals(self):
        archive=ROOT/'data/atmosphere/sources/condensation_atmosphere_x700_2800_g515'
        result=validate(archive)
        control=gas_state(archive/'gas-control',result)
        self.assertLess(control['flux_error'],.002)
        self.assertGreater(result['diagnostics']['Pgas'],control['Pgas'])

    def test_control_comparison_rejects_different_opacity_physics(self):
        archive=ROOT/'data/atmosphere/sources/condensation_atmosphere_x700_2800_g515'
        result=validate(archive)
        original=json.loads(gzip.decompress((archive/'gas-control/opacity-control.json.gz').read_bytes()))
        for key,value in [('mode','equilibrium'),('XH',.45)]:
            opacity=json.loads(json.dumps(original));opacity['provenance'][key]=value
            with self.subTest(key=key),self.assertRaisesRegex(ValueError,'opacity source physics'):
                gas_state(archive/'gas-control',result,opacity)
        opacity=json.loads(json.dumps(original))
        opacity['provenance']['specification']['opacity_frequencies']=15000
        with self.assertRaisesRegex(ValueError,'opacity source physics'):
            gas_state(archive/'gas-control',result,opacity)

    def test_diagnostic_grain_transport_does_not_accept_a_production_model(self):
        archive=ROOT/'data/atmosphere/sources/condensation_material_refinement/refined-2800-diagnostic'
        with self.assertRaisesRegex(ValueError,'grain enthalpy'):
            validate(archive)
        diagnostic=validate(archive,diagnostic_grain_transport=True)
        self.assertFalse(diagnostic['chemistry']['grain_enthalpy_supported'])
        self.assertGreater(diagnostic['chemistry']['maximum_convective_flux_fraction_in_condensing_layers'],1e-8)
        with tempfile.TemporaryDirectory() as temp,self.assertRaisesRegex(ValueError,'grain enthalpy'):
            archive_model(archive,Path(temp)/'production')

    def test_thermal_controls_replay_original_source_outputs(self):
        archive=ROOT/'data/atmosphere/sources/condensation_material_refinement/refined-2800-diagnostic'
        reference=validate(archive,diagnostic_grain_transport=True)
        result=thermal_sensitivity(ROOT/'docs/results/condensation_warm_heat_capacity_controls.json',reference,archive)
        self.assertEqual(len(result['controls']),2)
        self.assertLess(abs(result['controls'][1]['relative_boundary_difference']['Pgas']),1e-5)
        self.assertFalse(reference['chemistry']['grain_enthalpy_supported'])

    def test_thermal_diagnostic_requires_matching_reference_and_both_controls(self):
        archive=ROOT/'data/atmosphere/sources/condensation_material_refinement/refined-2800-diagnostic'
        reference=validate(archive,diagnostic_grain_transport=True)
        original=json.loads((ROOT/'docs/results/condensation_warm_heat_capacity_controls.json').read_text())
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'control.json'
            altered={**original,'reference_source_sha256':'0'*64};path.write_text(json.dumps(altered))
            with self.assertRaisesRegex(ValueError,'different reference'):
                thermal_sensitivity(path,reference,archive)
            altered={**original,'controls':original['controls'][:1]};path.write_text(json.dumps(altered))
            with self.assertRaisesRegex(ValueError,'both specified'):
                thermal_sensitivity(path,reference,archive)

    def test_hot_layers_have_an_explicit_independently_checked_join(self):
        archive=ROOT/'data/atmosphere/sources/condensation_validation/vaporized_x300_3200_g490'
        result=validate(archive);join=result['chemistry']['vaporized_continuation']
        self.assertGreater(join['actual_temperature_range_K'][0],6000)
        self.assertEqual(join['independently_queried_temperature_K'],6000)
        self.assertEqual(join['maximum_join_condensate_particles_per_nucleus'],0)
        self.assertTrue(result['chemistry']['grain_enthalpy_supported'])
        chemistry=json.loads(gzip.decompress((archive/'chemistry-audit-source.json.gz').read_bytes()))
        self.assertTrue(all(r['T_K']<=6000 for r in chemistry['rows']))
        self.assertEqual(len(chemistry['rows'])+join['layers'],300)

    def test_hot_continuation_cannot_omit_its_independent_join(self):
        archive=ROOT/'data/atmosphere/sources/condensation_validation/vaporized_x300_3200_g490'
        with tempfile.TemporaryDirectory() as temp:
            target=Path(temp)/'model';shutil.copytree(archive,target)
            (target/'chemistry-hot-join-source.json.gz').unlink()
            with self.assertRaises(FileNotFoundError):validate(target)


if __name__=='__main__':unittest.main()
