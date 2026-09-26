#!/usr/bin/env python3
"""Recompute the warm X=.25 held-out interpolation errors from raw atmospheres."""
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from assemble_nongrey_grid import load_continuation
from atmosphere_capacity_identity import CapacityEquivalence


def main():
    output = Path('docs/results/fable_x025_atmosphere_review_v1.json')
    if output.exists(): raise ValueError('preserve previous review')
    manifest_path = Path('docs/research/fable/results/fable-atmosphere-x030-assembly-v1/ember-fable-nongrey-x030-candidate-v1.manifest.json')
    manifest = json.loads(manifest_path.read_text())
    capacity = CapacityEquivalence(manifest['capacity_replays']); inputs = dict(capacity.inputs)
    claimed_path = Path('docs/research/fable/transport/atmosphere/x025_heldout_interpolation_v1.json')
    claimed = json.loads(claimed_path.read_text()); rows = []
    for claim in claimed['rows']:
        t, g = claim['teff_K'], claim['log_g']; tag = f't{t}-g{round(100*g)}'
        paths = [f'/tmp/ember-nongrey-t5200-t6000-v1/x200-y000-{tag}',
                 f'/tmp/ember-fable-nongrey-x030-warm-v2/x250-y000-{tag}',
                 f'/tmp/ember-fable-nongrey-x030-warm-v1/x300-y000-{tag}']
        states = []
        for x, path in zip((.2, .25, .3), paths):
            key, state, spec, record, identity = load_continuation(path, manifest['physics'], capacity)
            assert key == (x, 0., t, g)
            states.append(state); inputs.update(identity)
        errors = {k: math.sqrt(states[0][k]*states[2][k])/states[1][k]-1
                  for k in ('T', 'Pgas', 'source_density', 'column_mass')}
        for k, error in errors.items(): assert abs(error-claim[k]['loglinear_error']) < 3e-14
        rows.append(dict(Teff=t, logg=g, sources=paths, states=states, interpolation_errors=errors))
        print(t, g, errors, flush=True)
    for p in [Path(__file__), manifest_path, claimed_path]: inputs[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
    output.write_text(json.dumps(dict(created_utc=datetime.now(timezone.utc).isoformat(),
        outcome='passed_raw_heldouts', cases=rows, sha256=inputs,
        scope='Two warm composition held-outs; not a bound over an entire atmosphere family.'), indent=2)+'\n')


if __name__ == '__main__': main()
