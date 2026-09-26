#!/usr/bin/env python3
"""Check the two C/N atmosphere follow-ups from raw source outputs and material tables."""
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import generate_nongrey_grid as gen
from import_nongrey_grid import read_text, source_inputs, source_state, source_diagnostics_match
from nongrey_opacity import read_table


def main():
    output = Path('docs/results/fable_cn_atmosphere_review_v2.json')
    if output.exists(): raise ValueError('preserve prior review')
    sealed = Path('docs/research/fable/results/fable-cn-coupled-v2')
    ready = json.loads((sealed/'READY.json').read_text())
    hashes = {}
    def pin(p):
        p = Path(p); h = hashlib.sha256()
        with p.open('rb') as f:
            for b in iter(lambda: f.read(1024*1024), b''): h.update(b)
        hashes[str(p)] = h.hexdigest(); return h.hexdigest()
    for file, expected in ready['files_sha256'].items(): assert pin(sealed/file) == expected
    original = gen.composition
    cases = [(.25, 3200, 'f025-x020-y000-t3200-g515',
              '/private/tmp/ember-nongrey-warm-extension-shallow-v1/x020-y000-t3200-g515'),
             (.5, 5400, 'f050-x020-y000-t5400-g515',
              '/tmp/ember-nongrey-t5200-t6000-v1/x200-y000-t5400-g515')]
    results = []
    for fraction, temperature, label, baseline in cases:
        root = Path('/tmp/ember-fable-cn-coupled-v1')
        table_path = root/(label+'_plane')/'fort.63'; table = read_table(table_path)
        pin(table_path)
        states = {}
        for name, directory, f in [('reference', Path(baseline), 0.), ('converted', root/label, fraction)]:
            receipt = directory/'final/validated.json'; pin(receipt)
            validated = json.loads(receipt.read_text())
            spec = validated['continuation_provenance']['specification']
            inputs = {}
            for k in ['log', 'convergence', 'atmosphere_input', 'element_masses', 'parameters', 'initial_structure']:
                if k not in validated: continue
                p = directory/validated[k]
                assert pin(p) == validated[k+'_sha256'], (label, k)
                inputs[k] = read_text(p)
            def converted(x, y, metals):
                abundance, masses = original(x, y, metals)
                moved = f*abundance[5]
                abundance[5] -= moved; abundance[6] += moved; abundance[1] -= moved/2
                return abundance, masses
            try:
                gen.composition = converted
                source_inputs(inputs, spec, .2, 0., temperature, 5.15, inputs['log'])
            finally: gen.composition = original
            if name == 'converted':
                expected, _ = converted(.2, 0., spec['metals'])
                for (_, eos, absorption), abundance in zip(table['abundance'], expected):
                    assert math.isclose(eos, abundance, rel_tol=1e-12)
                    assert math.isclose(absorption, abundance, rel_tol=1e-12)
                assert hashes[str(table_path)] == validated['continuation_provenance']['opacity_sha256']
                temperatures = [math.exp(v) for v in table['log_temperature']]
                densities = [math.exp(v) for v in table['log_density']]
            else:
                # The baseline material source was already accepted; recheck its
                # actual atmosphere against the recorded source rectangle.
                temperatures = spec.get('temperature_K')
                if temperatures is None:
                    n, lo, hi = spec['log_temperature']
                    temperatures = [10**(lo+i*(hi-lo)/(n-1)) for i in range(n)]
                n, lo, hi = spec['log_density']
                densities = [10**(lo+i*(hi-lo)/(n-1)) for i in range(n)]
            state = source_state(inputs['log'], inputs['convergence'], temperature, 5.15,
                                 {'temperature_K': temperatures, 'density_g_cm3': densities})
            assert source_diagnostics_match(validated['diagnostics'], state)
            states[name] = state
        relative = {k: states['converted'][k]/states['reference'][k]-1
                    for k in ['T', 'Pgas', 'source_density', 'column_mass']}
        results.append(dict(fraction=fraction, Teff=temperature, states=states, relative_changes=relative,
                            table_shape=table['shape'], decks_and_material_abundances_checked=True))
        print(label, relative, flush=True)
    pin(__file__)
    output.write_text(json.dumps(dict(created_utc=datetime.now(timezone.utc).isoformat(),
        outcome='passed_two_source_atmosphere_comparisons', cases=results, sha256=hashes,
        limitations=['Two composition/temperature points; not a full C/N atmosphere grid or a linearity bound.',
                     'Gas atmospheres with Rayleigh scattering; no condensates.',
                     'Reference material-source provenance reused; new material tables read and checked directly.']), indent=2)+'\n')


if __name__ == '__main__': main()
