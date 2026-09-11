"""Check source-row retention and refusal of incompatible density additions."""
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from metal_eos_composition import mixture


class DensityMergeTests(unittest.TestCase):
    def sources(self, root, numerical=True):
        options = [3, 223, -2] if numerical else [3, 1, -2]
        ts = [5.+.0125*i for i in range(5)]
        q = [[2.45+.0125*i for i in range(5)], [2.5+.0125*i for i in range(5)]]
        # Use an exact common endpoint; all spacing is checked by the importer.
        q[0][-1] = q[1][0]
        roots = [root/'previous', root/'addition']
        for directory, axis in zip(roots, q):
            directory.mkdir()
            spec = {'hydrogen': [.1], 'helium3': [0.], 'logT': ts, 'logQ': axis,
                    'probe_sha256': 'same-probe', 'source_archive_sha256': 'same-archive'}
            if numerical:
                spec.update(source_options=options, source_radiation_included=False)
            (directory/'specification.json').write_text(json.dumps(spec))
            plane = directory/'plane-000'; plane.mkdir()
            rows = [[0.]+[t+value+i for i in range(21)] for t in ts for value in axis]
            raw = {**mixture(.1, 0.), 'version': 'FreeEOS 3.0.0', 'options': options,
                   'logT': ts, 'logQ': axis, 'probe_sha256': 'same-probe',
                   'source_archive_sha256': 'same-archive', 'data': rows}
            (plane/'source.json.gz').write_bytes(gzip.compress(json.dumps(raw).encode()))
        return roots

    def run_merge(self, roots, output):
        return subprocess.run([sys.executable, '-B', str(ROOT/'scripts/merge_metal_eos_density_sources.py'),
                               *map(str, roots), str(output)], capture_output=True, text=True)

    def test_retains_old_rows_and_appends_each_isotherm(self):
        for numerical in [False, True]:
            with self.subTest(numerical=numerical), tempfile.TemporaryDirectory() as work:
                root = Path(work); roots = self.sources(root, numerical)
                original = [json.loads(gzip.decompress((p/'plane-000/source.json.gz').read_bytes())) for p in roots]
                result = self.run_merge(roots, root/'merged')
                self.assertEqual(result.returncode, 0, result.stderr)
                merged = json.loads(gzip.decompress((root/'merged/plane-000/source.json.gz').read_bytes()))
                self.assertEqual(merged['options'], original[0]['options'])
                self.assertEqual(len(merged['data']), 45)
                for i in range(5):
                    self.assertEqual(merged['data'][9*i:9*i+5], original[0]['data'][5*i:5*i+5])
                    self.assertEqual(merged['data'][9*i+5:9*i+9], original[1]['data'][5*i+1:5*i+5])

    def test_changed_shared_heat_capacity_is_rejected(self):
        with tempfile.TemporaryDirectory() as work:
            root = Path(work); roots = self.sources(root)
            path = roots[1]/'plane-000/source.json.gz'
            raw = json.loads(gzip.decompress(path.read_bytes()))
            raw['data'][10][13] *= 1.001
            path.write_bytes(gzip.compress(json.dumps(raw).encode()))
            result = self.run_merge(roots, root/'merged')
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('density overlap differs', result.stderr)
            self.assertFalse((root/'merged/plane-000/source.json.gz').exists())

    def test_changed_electron_treatment_is_rejected(self):
        with tempfile.TemporaryDirectory() as work:
            root = Path(work); roots = self.sources(root)
            path = roots[1]/'specification.json'
            spec = json.loads(path.read_text()); spec['source_options'] = [3, 1, -2]
            path.write_text(json.dumps(spec))
            result = self.run_merge(roots, root/'merged')
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('source physics or radiation treatment differs', result.stderr)
            self.assertFalse((root/'merged').exists())


if __name__ == '__main__':
    unittest.main()
