"""Verify an explicit, completed capacity replay before combining gas sources.

Actual executable identities remain in each source record. A replay establishes
capacity-only source changes, not grid convergence or general accuracy.
"""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from prepare_nongrey_sources import digest


class CapacityEquivalence:
    def __init__(self, reports=()):
        self.paths = [str(Path(p).resolve(strict=True)) for p in reports]
        self.inputs = {}
        self.variants = {}
        for name in self.paths:
            path = Path(name)
            report = json.loads(path.read_text())
            if report.get('passed') is not True or report.get('parameter') != 'MDEPTH':
                raise ValueError('a passing depth-capacity replay is required')
            inputs = report['input_sha256']
            prepared_paths = [Path(p) for p in inputs if Path(p).name == 'prepared.json']
            model_paths = [Path(p).parent for p in inputs if Path(p).name == 'provenance.json']
            if len(prepared_paths) != 1 or len(model_paths) != 1:
                raise ValueError('capacity replay must identify one preparation and reference model')
            for p, checksum in inputs.items():
                if digest(p) != checksum:
                    raise ValueError('capacity replay dependency changed: ' + p)
            # Reuse the full source/output check. This mode never launches TLUSTY
            # or rewrites an input; its sole output is a temporary check report.
            with tempfile.TemporaryDirectory() as temporary:
                checked = Path(temporary)/'rechecked.json'
                result = subprocess.run([
                    sys.executable, '-B', str(Path(__file__).with_name('check_atmosphere_capacity.py')),
                    str(prepared_paths[0]), str(model_paths[0]), report['work'],
                    str(checked), '--completed-replay'], capture_output=True, text=True)
                if result.returncode:
                    raise ValueError('completed capacity replay failed: ' + result.stderr)
                if json.loads(checked.read_text()) != report:
                    raise ValueError('capacity replay no longer reproduces its recorded checks')
            prepared = json.loads(prepared_paths[0].read_text())
            variant = prepared['array_capacity']
            base = variant['canonical_prepared']
            # Preparation may change only the executable path/hash and capacities.
            allowed = {'tlusty', 'executables', 'array_capacity'}
            if ({k: v for k, v in prepared.items() if k not in allowed}
                    != {k: v for k, v in base.items() if k not in allowed}
                    or set(prepared['executables']) != set(base['executables'])
                    or any(prepared['executables'][k] != v
                           for k, v in base['executables'].items() if k != 'tlusty')):
                raise ValueError('preparation changed beyond TLUSTY capacity')
            sha = prepared['executables']['tlusty']
            if (sha in self.variants or sha != report['new_executable_sha256']
                    or base['executables']['tlusty'] != report['original_executable_sha256']):
                raise ValueError('ambiguous capacity executable identity')
            self.variants[sha] = {'prepared': prepared, 'report': report, 'path': name}
            self.inputs.update(inputs)
            self.inputs[name] = digest(path)
            # Pin the actual Fortran files as well as the executable and replay.
            for source in (base, prepared):
                for p in Path(source['tlusty']).parent.iterdir():
                    if p.suffix.lower() in ('.f', '.for'):
                        self.inputs[str(p.resolve())] = digest(p)
        if self.paths:
            self.inputs[str(Path(__file__).resolve())] = digest(__file__)

    def canonical_preparation(self, spec, prepared):
        """Return the verified original preparation; reject an untested setup."""
        entry = self.variants.get(prepared['executables']['tlusty'])
        if entry is None or prepared != entry['prepared']:
            raise ValueError('source executable has no matching capacity replay')
        report = entry['report']
        if (spec['atmosphere_frequencies'] != report['frequencies']
                or not 2 <= spec['depths'] <= report['new_capacity']):
            raise ValueError('source sampling is outside the checked capacity setup')
        return prepared['array_capacity']['canonical_prepared']
