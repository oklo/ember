#!/usr/bin/env python3
"""Assemble validated gas atmospheres with explicit holes in source coverage.

Retains every original source recipe. Version-2 runtime tables reject any
interpolation or derivative stencil containing a missing source state.
This does not accept unfinished models, interpolate missing nodes, or
replace an independent source-grid refinement study.
"""
import argparse
from functools import lru_cache
import gzip
import itertools
import json
import math
from pathlib import Path
import tempfile

from generate_nongrey_grid import composition, temperatures, sequence, input_fingerprint
from import_nongrey_grid import GAS_CALCULATION, gas_model_state, import_grid
from nongrey_opacity import validate_table
from prepare_nongrey_sources import digest

KEYS = ['XH', 'X3', 'teff_K', 'log_g']
PHYSICS = ['metals', 'alpha', 'tau', 'wavelength_A', 'microturbulence_km_s', 'line_threshold']


@lru_cache(maxsize=64)
def _validate_shared_table(path, checksum, abundance, temperature, density):
    # Retain only successful validation, not the large arrays. The caller
    # checks the current bytes on every use, including a cache hit.
    validate_table(path, abundance, temperature, density)
    if digest(path) != checksum:
        raise ValueError('opacity changed during source validation')


def validate_shared_table(path, checksum, abundance, temperature, density):
    path = Path(path).resolve(strict=True)
    if digest(path) != checksum:
        raise ValueError('source opacity checksum mismatch')
    _validate_shared_table(str(path), checksum, tuple(abundance),
                           tuple(temperature), tuple(density))


def physical_identity(spec, prepared):
    if any(k in prepared for k in ['depletion', 'initialization_only']):
        raise ValueError('canonical gas source required')
    return {'settings': {k: spec[k] for k in PHYSICS},
            'executables': {k: prepared['executables'][k] for k in ['tlusty', 'synspec']},
            'data_sha256': prepared['data_sha256'],
            'line_list_sha256': prepared['line_list_sha256'],
            'opacity_method': prepared['opacity_method']}


def complete_cells(axes, states):
    result = []
    for base in itertools.product(*(range(len(a)-1) for a in axes)):
        choices = [a[i:i+2] for a, i in zip(axes, base)]
        if all(corner in states for corner in itertools.product(*choices)):
            result.append(choices)
    return result


def load_continuation(work, identity=None):
    """Revalidate archived canonical outputs without rerunning or rewriting them."""
    work = Path(work)
    inputs = {}
    provenance = work/'provenance.json'
    validation = work/'final/validated.json'
    receipt = work/'final/completed.json'
    record = json.loads(validation.read_text())
    source = json.loads(provenance.read_text())
    spec, prepared = source['specification'], source['prepared']
    actual_identity = physical_identity(spec, prepared)
    if identity is not None and actual_identity != identity:
        raise ValueError('source physics differs; do not silently mix atmosphere prescriptions')
    key = tuple(record[k] for k in KEYS)
    if any(record[k] != source[k] for k in KEYS):
        raise ValueError('continuation label mismatch')
    if record.get('continuation_provenance') != source:
        raise ValueError('final validation belongs to a different continuation')
    saved = json.loads(receipt.read_text())
    final = work/'final'
    if saved['input_sha256'] != input_fingerprint(prepared['executables']['tlusty'], final, archived=True):
        raise ValueError('canonical final source input fingerprint changed')
    if not {'fort.7', 'fort.9', 'run.log'}.issubset(saved['outputs']):
        raise ValueError('incomplete final source receipt')
    import hashlib
    for name, expected in saved['outputs'].items():
        path = final/name
        data = path.read_bytes() if path.exists() else gzip.decompress(path.with_name(name+'.gz').read_bytes())
        if hashlib.sha256(data).hexdigest() != expected:
            raise ValueError('final source output changed')
    opacity = final/'opacity.bin'
    if (final/'opacity.sha256').read_text().strip() != source['opacity_sha256']:
        raise ValueError('final atmosphere opacity identity changed')
    abundance, _ = composition(*key[:2], spec['metals'])
    validate_shared_table(opacity, source['opacity_sha256'], abundance,
                          temperatures(spec), sequence(spec['log_density']))
    own_spec = {**spec, 'opacity': {'temperature_K': temperatures(spec),
                                  'density_g_cm3': sequence(spec['log_density'])}}
    state = gas_model_state(work, own_spec, record)
    for path in [provenance, validation, receipt, opacity.resolve()]:
        inputs[str(path.resolve())] = digest(path)
    for kind in ['log', 'convergence', 'atmosphere_input', 'element_masses', 'parameters', 'initial_structure']:
        if kind in record:
            inputs[str((work/record[kind]).resolve())] = record[kind+'_sha256']
    for name, expected in record["initialization"]["files_sha256"].items():
        path = work/name
        if digest(path) != expected:
            raise ValueError("archived initializer changed")
        inputs[str(path.resolve())] = expected
    return key, state, spec, record, inputs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('base_manifest', type=Path)
    p.add_argument('output', type=Path)
    group=p.add_mutually_exclusive_group(required=True)
    group.add_argument('--continuations', type=Path, nargs='+',
                   help='completed run_nongrey_continuation work directories')
    group.add_argument('--continuation-plan', type=Path,
                       help='JSON with an explicit continuations list; no source selection or automatic adoption')
    a = p.parse_args()
    if a.continuation_plan:
        plan=json.loads(a.continuation_plan.read_text())
        if not isinstance(plan.get('continuations'),list) or not plan['continuations']:
            raise ValueError('assembly plan requires an explicit nonempty continuation list')
        a.continuations=[Path(v) for v in plan['continuations']]
    report_path = a.output.with_suffix('.manifest.json')
    if a.output.exists() or report_path.exists():
        raise FileExistsError('use a new candidate table and manifest')
    base = json.loads(a.base_manifest.read_text())
    if base.get('calculation') != GAS_CALCULATION:
        raise ValueError('base must be a canonical gas family')
    identity = physical_identity(base, base['provenance'])
    inputs = {str(a.base_manifest.resolve()): digest(a.base_manifest)}
    if a.continuation_plan:inputs[str(a.continuation_plan.resolve())]=digest(a.continuation_plan)
    with tempfile.TemporaryDirectory() as temporary:
        states = import_grid(a.base_manifest, Path(temporary)/'base.dat')
    base_keys = set(states)
    for model in base["models"]:
        for kind in ["log", "convergence", "atmosphere_input", "element_masses", "parameters", "initial_structure"]:
            if kind in model:
                path = a.base_manifest.parent/model[kind]
                inputs[str(path.resolve())] = model[kind+"_sha256"]
    accepted = []
    for work in a.continuations:
        key, state, spec, record, dependencies = load_continuation(work, identity)
        if key in states:
            raise ValueError("duplicate source coordinate")
        states[key] = state
        inputs.update(dependencies)
        validation = work/"final/validated.json"
        accepted.append({'coordinates': key, 'work': str(work.resolve()),
                         'state': states[key], 'source_specification': spec,
                         'validation_sha256': digest(validation)})
    axes = [sorted({k[i] for k in states}) for i in range(4)]
    if any(len(axis) < 2 for axis in axes) or axes[0][-1]+axes[1][-1]+sum(base['metals']) > 1+1e-12:
        raise ValueError('invalid output axes')
    cells = complete_cells(axes, states)
    if not cells:
        raise ValueError('no complete interpolation cells')
    lines = ['EMBER_COMPOSITION_ATMOSPHERE 2',
             'source '+json.dumps(base['source']+'; validated extension with explicit missing nodes'),
             'approximation '+json.dumps(base['approximation']), 'basis baryon_mass',
             f'tau {base["tau"]:.17g}', 'metals '+' '.join(f'{v:.17g}' for v in base['metals'])]
    for name, axis in zip(['hydrogen', 'helium3', 'log_teff', 'log_g'], axes):
        values = [math.log10(v) for v in axis] if name == 'log_teff' else axis
        lines.append(f'{name} {len(axis)} '+' '.join(f'{v:.17g}' for v in values))
    lines.append('data')
    for key in itertools.product(*axes):
        if key not in states:
            lines.append('0')
        else:
            r = states[key]
            lines.append(f'1 {math.log10(r["T"]):.17g} {math.log10(r["Pgas"]):.17g}')
    if any(digest(path) != expected for path, expected in inputs.items()):
        raise ValueError('source inputs changed during assembly')
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with a.output.open('x') as f:
        f.write('\n'.join(lines)+'\n')
    report = {'scope': __doc__, 'base_manifest': str(a.base_manifest.resolve()),
              'base_states': len(base_keys), 'accepted_states': len(states),
              'missing_states': math.prod(map(len, axes))-len(states),
              'axis_order': KEYS, 'axes': axes, 'complete_cells': cells,
              'physics': identity, 'extension': accepted, 'input_files_sha256': inputs,
              'table_sha256': digest(a.output), 'script_sha256': digest(__file__)}
    with report_path.open('x') as f:
        f.write(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: report[k] for k in ['base_states', 'accepted_states', 'missing_states', 'table_sha256']}))


if __name__ == '__main__':
    main()
