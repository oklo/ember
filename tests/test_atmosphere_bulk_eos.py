"""Prevent different atmosphere EOS treatments from sharing a source identity."""
import copy
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from atmosphere_bulk_eos import bulk_identity, validate_bulk_eos
from assemble_nongrey_grid import physical_identity


class BulkEosIdentity(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        self.table = self.directory/'bulk.dat'
        self.table.write_text('EMBER_ATMOSPHERE_BULK_EOS 1\n2 2 .98\n7 8\n12 13\n'
                              + ('0 1 1 0 '*3+'\n')*4)
        self.spec = dict(hydrogen=[.98], helium3=[0.], metals=[0.]*5, alpha=1.9,
                         tau=100, wavelength_A=[900, 300000], microturbulence_km_s=1.,
                         line_threshold=.0001)
        self.prepared = dict(executables=dict(tlusty='t', synspec='s'), data_sha256='d',
                             line_list_sha256=['l'], opacity_method='sampling')
        self.bulk = dict(format=1, family=dict(model='test bulk EOS',
                        source_sha256={'H': 'a'*64}, recipe_sha256={'interpolation': 'b'*64}),
                        table='bulk.dat', table_sha256=hashlib.sha256(self.table.read_bytes()).hexdigest())
        self.prepared['bulk_eos'] = self.bulk
        self.log = ' EMBER NONIDEAL BULK EOS: /original/location/bulk.dat 0.98\n'

    def test_family_identity_excludes_composition_specific_filename(self):
        changed = copy.deepcopy(self.prepared)
        changed['bulk_eos']['table'] = 'other.dat'
        changed['bulk_eos']['table_sha256'] = 'c'*64
        self.assertEqual(physical_identity(self.spec, self.prepared), physical_identity(self.spec, changed))
        changed['bulk_eos']['family']['source_sha256']['H'] = 'd'*64
        self.assertNotEqual(physical_identity(self.spec, self.prepared), physical_identity(self.spec, changed))

    def test_ordinary_identity_is_preserved_and_distinct(self):
        ordinary = {k: v for k, v in self.prepared.items() if k != 'bulk_eos'}
        identity = physical_identity(self.spec, ordinary)
        self.assertNotIn('bulk_eos', identity)
        self.assertNotEqual(identity, physical_identity(self.spec, self.prepared))
        self.assertIsNone(validate_bulk_eos(self.directory, self.spec, ordinary, 'ordinary atmosphere'))
        with self.assertRaisesRegex(ValueError, 'unrecorded'):
            validate_bulk_eos(self.directory, self.spec, ordinary, self.log)
        ordinary['bulk_eos_experiment'] = {'enabled': True}
        with self.assertRaisesRegex(ValueError, 'reviewed source identity'):
            bulk_identity(ordinary)

    def test_loaded_table_and_composition_are_verified(self):
        result = validate_bulk_eos(self.directory, self.spec, self.prepared, self.log)
        self.assertEqual(result['hydrogen'], .98)
        for spec, log in [(dict(self.spec, hydrogen=[.99]), self.log),
                          (dict(self.spec, helium3=[.001]), self.log),
                          (dict(self.spec, metals=[.01]*5), self.log),
                          (self.spec, ''), (self.spec, self.log.replace('0.98', '0.99'))]:
            with self.assertRaises(ValueError):
                validate_bulk_eos(self.directory, spec, self.prepared, log)
        self.table.write_text(self.table.read_text().replace('12 13', '12 14'))
        with self.assertRaisesRegex(ValueError, 'checksum'):
            validate_bulk_eos(self.directory, self.spec, self.prepared, self.log)

    def test_masked_bulk_table_format(self):
        text = self.table.read_text().replace('BULK_EOS 1', 'BULK_EOS 2')
        self.table.write_text(text+'cells\n1\n')
        self.bulk['table_sha256'] = hashlib.sha256(self.table.read_bytes()).hexdigest()
        self.assertEqual(validate_bulk_eos(self.directory, self.spec, self.prepared,
                                          self.log)['hydrogen'], .98)
        for suffix in ['', 'wrong\n1\n', 'cells\n2\n', 'cells\n1 0\n']:
            self.table.write_text(text+suffix)
            self.bulk['table_sha256'] = hashlib.sha256(self.table.read_bytes()).hexdigest()
            with self.assertRaisesRegex(ValueError, 'cell mask'):
                validate_bulk_eos(self.directory, self.spec, self.prepared, self.log)

    def test_bad_table_is_rejected_even_with_a_matching_checksum(self):
        for text in ['bad header\n', 'EMBER_ATMOSPHERE_BULK_EOS 1\n2 2 .98\n7 8\n12 13\n',
                     self.table.read_text().replace('7 8', '8 7'),
                     self.table.read_text().replace('0 1 1 0', 'nan 1 1 0')]:
            self.table.write_text(text)
            self.bulk['table_sha256'] = hashlib.sha256(self.table.read_bytes()).hexdigest()
            with self.assertRaises(ValueError):
                validate_bulk_eos(self.directory, self.spec, self.prepared, self.log)


if __name__ == '__main__':
    unittest.main()
