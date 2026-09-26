#!/usr/bin/env python3
"""Assemble independently checked H/helium atmosphere columns on fraction axes.

Physical source recipes retain XHe3. Only the runtime coordinate is
f3 = XHe3/(XHe3+XHe4). Missing columns remain explicit holes. An output table
requires at least one complete interpolation cell; this is not an independent
validation of interpolation accuracy or permission to select a stellar boundary.
"""
import argparse
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path

from assemble_nongrey_grid import complete_cells, physical_identity, validate_shared_table
from generate_nongrey_grid import composition, input_fingerprint, sequence, temperatures
from import_nongrey_grid import gas_model_state, source_diagnostics_match
from prepare_nongrey_sources import digest

NUMERICS = ['depths', 'atmosphere_frequencies', 'opacity_frequencies',
            'temperature_step_limit', 'newton_relaxation', 'convection_derivative_step',
            'tau_top', 'tau_bottom', 'log_temperature', 'log_density']


def axis(values, name, lower=None, upper=None):
    if (not isinstance(values, list) or len(values) < 2
            or any(not math.isfinite(v) for v in values)
            or any(a >= b for a, b in zip(values, values[1:]))
            or (lower is not None and values[0] < lower)
            or (upper is not None and values[-1] > upper)):
        raise ValueError(f'invalid {name} axis')
    return values


def physical_helium3(hydrogen, fraction, metals):
    helium = 1 - hydrogen - sum(metals)
    if not (0 <= hydrogen < 1 and 0 <= fraction <= 1 and helium > 0):
        raise ValueError('unphysical hydrogen/helium-fraction coordinates')
    return helium*fraction


def load_plane(work, hydrogen, fraction):
    """Reparse completed generator columns with exact original inputs/receipts."""
    work = Path(work).resolve(strict=True)
    spec_path, prepared_path = work/'specification.json', work/'provenance.json'
    spec, prepared = (json.loads(p.read_text()) for p in [spec_path, prepared_path])
    metals = spec['metals']
    if len(metals) != 5 or any(not math.isfinite(v) or v < 0 for v in metals):
        raise ValueError('invalid source metal vector')
    x3 = physical_helium3(hydrogen, fraction, metals)
    if spec['hydrogen'] != [hydrogen] or spec['helium3'] != [x3]:
        raise ValueError('plane recipe differs from physical composition; f3 is not XHe3')
    opacity = work/'plane-000/opacity/fort.63'
    opacity_hash = digest(opacity)
    abundance, _ = composition(hydrogen, x3, metals)
    validate_shared_table(opacity, opacity_hash, abundance,
                          temperatures(spec), sequence(spec['log_density']))
    source_spec = dict(spec, opacity={'temperature_K': temperatures(spec),
                                    'density_g_cm3': sequence(spec['log_density'])})
    inputs = {str(p): digest(p) for p in [spec_path, prepared_path, opacity]}
    identity = dict(physical_identity(spec, prepared),
                    numerical_settings={k: spec.get(k) for k in NUMERICS})
    states, records = {}, []
    for validation in sorted(work.glob('plane-000/model-*/validated.json')):
        directory = validation.parent
        record = json.loads(validation.read_text())
        if record['XH'] != hydrogen or record['X3'] != x3:
            raise ValueError('validated column has a different physical composition')
        teff, logg = record['teff_K'], record['log_g']
        if teff not in spec['teff_K'] or logg not in spec['log_g']:
            raise ValueError('column is not on its physical source axes')
        receipt_path = directory/'completed.json'
        receipt = json.loads(receipt_path.read_text())
        if receipt['input_sha256'] != input_fingerprint(
                prepared['executables']['tlusty'], directory, archived=True):
            raise ValueError('completed source input fingerprint mismatch')
        if not {'fort.7', 'fort.9', 'run.log'} <= receipt['outputs'].keys():
            raise ValueError('incomplete source output receipt')
        for name, expected in receipt['outputs'].items():
            path = directory/name
            if path.parent != directory:
                raise ValueError('unexpected source receipt path')
            if path.exists():
                contents = path.read_bytes()
            else:
                path = path.with_name(name+'.gz')
                contents = gzip.decompress(path.read_bytes())
            if hashlib.sha256(contents).hexdigest() != expected:
                raise ValueError('source output receipt checksum mismatch')
            inputs[str(path)] = digest(path)
        marker = directory/'opacity.sha256'
        if marker.read_text().strip() != opacity_hash:
            raise ValueError('column used a different opacity table')
        for kind in ['log', 'convergence', 'atmosphere_input', 'element_masses',
                     'parameters', 'initial_structure']:
            if kind in record:
                path = (work/record[kind]).resolve(strict=True)
                if not path.is_relative_to(work):
                    raise ValueError('source record leaves its archive')
                inputs[str(path)] = record[kind+'_sha256']
        state = gas_model_state(work, source_spec, record)
        if not source_diagnostics_match(record['diagnostics'], state):
            raise ValueError('recorded source diagnostics changed')
        key = (hydrogen, fraction, teff, logg)
        if key in states:
            raise ValueError('duplicate source column')
        states[key] = state
        records.append(dict(coordinates=key, physical_XHe3=x3, validation=str(validation),
                            state=state, source_record=record))
        for p in [validation, receipt_path, marker]:
            inputs[str(p)] = digest(p)
    return states, records, identity, inputs, spec


def render_table(plan, states, metals, approximation):
    axes = [axis(plan['hydrogen'], 'hydrogen', 0, 1),
            axis(plan['helium3_fraction'], 'helium fraction', 0, 1),
            axis(plan['teff_K'], 'effective temperature', 1),
            axis(plan['log_g'], 'log gravity')]
    physical_helium3(axes[0][-1], axes[1][-1], metals)
    tolerance = plan['metal_tolerance']
    if not math.isfinite(tolerance) or not 0 <= tolerance <= 1e-12:
        raise ValueError('explicit metal tolerance must be within [0,1e-12]')
    if not plan['source'] or not approximation or not math.isfinite(plan['tau']) or plan['tau'] <= 0:
        raise ValueError('source description, approximations and depth required')
    expected = set(itertools.product(*axes))
    if not states.keys() <= expected:
        raise ValueError('unexpected source coordinate')
    cells = complete_cells(axes, states)
    if not cells:
        raise ValueError('no complete interpolation cell; preserve columns and await missing work')
    lines = ['EMBER_COMPOSITION_ATMOSPHERE 3', 'source '+json.dumps(plan['source']),
             'approximation '+json.dumps(approximation), 'basis baryon_mass',
             f'tau {plan["tau"]:.17g}',
             'metals '+' '.join(f'{v:.17g}' for v in metals),
             f'metal_tolerance {tolerance:.17g}']
    for name, values in zip(['hydrogen', 'helium3_fraction', 'log_teff', 'log_g'], axes):
        a = [math.log10(v) for v in values] if name == 'log_teff' else values
        lines.append(f'{name} {len(a)} '+' '.join(f'{v:.17g}' for v in a))
    lines.append('data')
    for key in itertools.product(*axes):
        if key not in states:
            lines.append('0')
        else:
            t, pg = states[key]['T'], states[key]['Pgas']
            if not all(math.isfinite(v) and v > 0 for v in [t, pg]):
                raise ValueError('invalid physical matching state')
            lines.append(f'1 {math.log10(t):.17g} {math.log10(pg):.17g}')
    return '\n'.join(lines)+'\n', cells


def assemble(plan_path, output):
    plan_path, output = Path(plan_path).resolve(strict=True), Path(output)
    manifest = output.with_suffix('.manifest.json')
    if output.exists() or manifest.exists():
        raise FileExistsError('use a new candidate table and manifest')
    plan = json.loads(plan_path.read_text())
    if plan.get('format') != 1:
        raise ValueError('unrecognized assembly plan')
    states, records, inputs, plane_keys = {}, [], {str(plan_path): digest(plan_path)}, set()
    identity = approximation = metals = None
    for plane in plan['planes']:
        key = plane['hydrogen'], plane['helium3_fraction']
        if key in plane_keys:
            raise ValueError('duplicate composition plane')
        plane_keys.add(key)
        own, columns, physics, pins, spec = load_plane(plane['work'], *key)
        if identity is None:
            identity, approximation, metals = physics, spec['approximation'], spec['metals']
        if physics != identity or spec['approximation'] != approximation or spec['tau'] != plan['tau']:
            raise ValueError('inconsistent atmosphere physics or numerical resolution')
        if set(states) & own.keys():
            raise ValueError('duplicate physical source')
        states.update(own); records.extend(columns); inputs.update(pins)
    if plane_keys != set(itertools.product(plan['hydrogen'], plan['helium3_fraction'])):
        raise ValueError('assembly plan must declare every composition plane')
    text, cells = render_table(plan, states, metals, approximation)
    if any(digest(p) != h for p, h in inputs.items()):
        raise ValueError('source changed during assembly')
    report = dict(production_selected=False, interpolation_accuracy_validated=False,
                  axis_order=['XH', 'He3/(He3+He4)', 'Teff_K', 'log_g'],
                  plan=plan, physical_identity=identity, columns=records,
                  complete_cells=cells, accepted_states=len(states),
                  missing_states=math.prod(len(plan[k]) for k in
                      ['hydrogen', 'helium3_fraction', 'teff_K', 'log_g'])-len(states),
                  input_files_sha256=inputs, script_sha256=digest(__file__))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x') as f:
        f.write(text)
    report['table_sha256'] = digest(output)
    with manifest.open('x') as f:
        json.dump(report, f, indent=2, allow_nan=False); f.write('\n')
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('plan', type=Path)
    p.add_argument('output', type=Path)
    args = p.parse_args()
    report = assemble(args.plan, args.output)
    print(json.dumps({k: report[k] for k in ['accepted_states', 'missing_states', 'table_sha256']}))
