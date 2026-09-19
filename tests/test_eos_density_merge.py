"""Check source-row retention and refusal of incompatible density additions."""
import gzip
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from metal_eos_composition import mixture
from eos_source_coverage import absent_source_rows, inconsistent_source_rows


class DensityMergeTests(unittest.TestCase):
    def sources(self, root, numerical=True, hot_subset=False):
        options = [3, 223, -2] if numerical else [3, 1, -2]
        ts = [5.+.0125*i for i in range(9 if hot_subset else 5)]
        q = [[2.45+.0125*i for i in range(5)], [2.5+.0125*i for i in range(5)]]
        # Use an exact common endpoint; all spacing is checked by the importer.
        q[0][-1] = q[1][0]
        roots = [root/'previous', root/'addition']
        for directory, axis in zip(roots, q):
            directory.mkdir()
            thermal = ts[-5:] if hot_subset and directory == roots[1] else ts
            spec = {'hydrogen': [.1], 'helium3': [0.], 'logT': thermal, 'logQ': axis,
                    'probe_sha256': 'same-probe', 'source_archive_sha256': 'same-archive'}
            if numerical:
                spec.update(source_options=options, source_radiation_included=False)
            (directory/'specification.json').write_text(json.dumps(spec))
            plane = directory/'plane-000'; plane.mkdir()
            rows = [[0.]+[t+value+i for i in range(21)] for t in thermal for value in axis]
            if directory == roots[1]:
                for row in rows:
                    row[1] += 7  # Iteration count may differ; physics must not.
            raw = {**mixture(.1, 0.), 'version': 'FreeEOS 3.0.0', 'options': options,
                   'logT': thermal, 'logQ': axis, 'probe_sha256': 'same-probe',
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

    def test_hot_addition_retains_precision_provenance(self):
        with tempfile.TemporaryDirectory() as work:
            root = Path(work); roots = self.sources(root, hot_subset=True)
            spec_path = roots[1]/'specification.json'
            spec = json.loads(spec_path.read_text())
            precision = {'probe_sha256': 'tighter-probe', 'fallback_relative_integral_error': 1e-11}
            spec['precision_fallback'] = precision
            spec_path.write_text(json.dumps(spec))
            path = roots[1]/'plane-000/source.json.gz'
            source = json.loads(gzip.decompress(path.read_bytes()))
            source['precision_fallback'] = precision
            source['precision_fallback_isotherms'] = [dict(
                temperature_index=1, logT=spec['logT'][1],
                actual_probe_input_sha256='source-request', nominal_failure_sha256='source-failure')]
            path.write_bytes(gzip.compress(json.dumps(source).encode()))
            result = self.run_merge(roots, root/'merged')
            self.assertEqual(result.returncode, 0, result.stderr)
            merged = json.loads(gzip.decompress((root/'merged/plane-000/source.json.gz').read_bytes()))
            merged_spec = json.loads((root/'merged/specification.json').read_text())
            self.assertEqual(merged['precision_fallback'], precision)
            self.assertEqual(merged_spec['precision_fallback'], precision)
            override, = merged['precision_fallback_isotherms']
            self.assertEqual(override['temperature_index'], 5)
            self.assertEqual(override['source_temperature_index'], 1)
            self.assertEqual(override['actual_probe_input_sha256'], 'source-request')
            self.assertEqual(override['retained_logQ_range'], [spec['logQ'][1], spec['logQ'][-1]])
            source['precision_fallback']['probe_sha256'] = 'unrecorded-probe'
            path.write_bytes(gzip.compress(json.dumps(source).encode()))
            result = self.run_merge(roots, root/'rejected')
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('precision differs', result.stderr)

    def test_density_and_temperature_offsets_preserve_exclusions(self):
        with tempfile.TemporaryDirectory() as work:
            root = Path(work); roots = self.sources(root, hot_subset=True)
            for directory, index in zip(roots, [12, 17]):
                path = directory/'specification.json'
                spec = json.loads(path.read_text()); spec['source_consistency_exclusion_limit'] = 1e-7
                path.write_text(json.dumps(spec))
                path = directory/'plane-000/source.json.gz'
                source = json.loads(gzip.decompress(path.read_bytes()))
                row = source['data'][index]
                rho, temperature, pressure = row[2:5]; chit, er, et, sr, st = row[8:13]
                defect = max(abs(rho*er/pressure+chit-1), abs(temperature*st/et-1),
                             abs(rho*temperature*sr/pressure+chit))
                source['source_consistency_exclusions'] = [dict(
                    index=index, reason='first_law_defect', criterion=1e-7,
                    maximum_absolute_defect=defect)]
                path.write_bytes(gzip.compress(json.dumps(source).encode()))
            result = self.run_merge(roots, root/'merged')
            self.assertEqual(result.returncode, 0, result.stderr)
            merged = json.loads(gzip.decompress((root/'merged/plane-000/source.json.gz').read_bytes()))
            self.assertEqual(inconsistent_source_rows(merged), {20, 69})
            spec = json.loads((root/'merged/specification.json').read_text())
            self.assertEqual(spec['source_consistency_exclusion_limit'], 1e-7)
            # An exclusion may not hide a coordinate whose actual response differs.
            merged['data'][20][11] *= 2
            with self.assertRaisesRegex(ValueError, 'does not match'):
                inconsistent_source_rows(merged)

    def test_hot_addition_has_declared_absence_and_retains_source(self):
        with tempfile.TemporaryDirectory() as work:
            root = Path(work); roots = self.sources(root, hot_subset=True)
            original = json.loads(gzip.decompress((roots[0]/'plane-000/source.json.gz').read_bytes()))
            result = self.run_merge(roots, root/'merged')
            self.assertEqual(result.returncode, 0, result.stderr)
            raw = json.loads(gzip.decompress((root/'merged/plane-000/source.json.gz').read_bytes()))
            self.assertEqual(sum(absent_source_rows(raw)), 16)
            for i in range(9):
                self.assertEqual(raw['data'][9*i:9*i+5], original['data'][5*i:5*i+5])
            for index, value in [(0, None), (5, [0.]*22)]:
                changed = json.loads(json.dumps(raw)); changed['data'][index] = value
                with self.assertRaisesRegex(ValueError, 'declared absent coverage'):
                    absent_source_rows(changed)
            del raw['source_coverage']
            with self.assertRaisesRegex(ValueError, 'without a coverage declaration'):
                absent_source_rows(raw)

    def test_absent_source_masks_full_potential_stencil(self):
        # Analytic ideal gas is a thermodynamically exact importer fixture.
        # The test checks that absent hot-extension cells cannot become usable
        # through neighboring derivative stencils, and preserves old jets.
        with tempfile.TemporaryDirectory() as work:
            root = Path(work)
            ts = [6.+.0125*i for i in range(11)]
            qs = [2.+.0125*i for i in range(11)]
            raw = {**mixture(.1, 0.), 'version': 'FreeEOS 3.0.0', 'options': [3, 223, -2],
                   'logT': ts, 'logQ': qs, 'data': []}
            for t in ts:
                for q in qs:
                    temperature = 10**t; rho = 10**(q+1.5*(t-6))
                    gas = 1.e8; cv = 1.5*gas; e = cv*temperature
                    entropy = cv*math.log(temperature)-gas*math.log(rho)
                    raw['data'].append([0., 1., rho, temperature, gas*rho*temperature,
                                       e, entropy, 1., 1., 0., e, -gas, cv, 2.5*gas,
                                       .4, 1., 1., 1., 0., 0., 0., 0.])
            source = root/'source.gz'; output = root/'potential.dat'
            def run_import():
                source.write_bytes(gzip.compress(json.dumps(raw).encode()))
                result = subprocess.run([sys.executable, '-B', str(ROOT/'scripts/import_freeeos_potential.py'),
                                         str(source), str(output)], capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                return output.read_text().split('data\n')[1].splitlines()
            complete = run_import()
            # A returned but inconsistent state is distinct from an absent
            # source. Its actual values and source flag must remain intact.
            k = 5*11+5
            raw['data'][k][11] *= 1.001
            r = raw['data'][k]
            defect = abs(r[2]*r[3]*r[11]/r[4]+r[8])
            raw['source_consistency_exclusions'] = [dict(
                index=k, reason='first_law_defect', criterion=1e-7,
                maximum_absolute_defect=defect)]
            self.assertEqual(inconsistent_source_rows(raw), {k})
            excluded = run_import()
            for i in range(2, 9):
                for j in range(2, 9):
                    index = (i-2)*7+j-2
                    valid = not ((i == 5 and abs(j-5) <= 2) or (j == 5 and abs(i-5) <= 2))
                    self.assertEqual(excluded[index].startswith('1 '), valid)
                    if valid:
                        self.assertEqual(excluded[index], complete[index])
            self.assertEqual(raw['data'][k][0], 0.)
            raw['source_consistency_exclusions'][0]['index'] = k+1
            with self.assertRaisesRegex(ValueError, 'does not match its defect'):
                inconsistent_source_rows(raw)
            del raw['source_consistency_exclusions']
            # Restore the analytic value exactly before testing absent cells.
            raw['data'][k][11] = -1.e8
            raw['source_coverage'] = dict(kind='hot_density_extension', original_logQ_max=qs[6],
                                          minimum_added_logT=ts[4])
            for i, t in enumerate(ts):
                for j, q in enumerate(qs):
                    if t < ts[4] and q > qs[6]:
                        raw['data'][i*11+j] = None
            masked = run_import()
            for i in range(2, 9):
                for j in range(2, 9):
                    index = (i-2)*7+j-2
                    expected = all(not (ts[i] < ts[4] and qs[j+d] > qs[6])
                                   and not (ts[i+d] < ts[4] and qs[j] > qs[6]) for d in range(-2, 3))
                    self.assertEqual(masked[index].startswith('1 '), expected)
                    if expected:
                        self.assertEqual(masked[index], complete[index])


if __name__ == '__main__':
    unittest.main()
