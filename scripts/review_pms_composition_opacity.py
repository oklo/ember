"""Check complete PMS composition opacity planes against every source isotherm."""
import datetime
import json
from pathlib import Path
import sys

import numpy as np

from generate_nongrey_grid import composition, input_fingerprint, temperatures, sequence
from nongrey_opacity import read_table, validate_table
from prepare_nongrey_sources import digest


def review(work):
    work = Path(work)
    plan = json.loads((work / 'plan.json').read_text())
    final = json.loads((work / 'final.json').read_text())
    assert final['all_isotherms_pass'] and final['completed'] == final['planned_isotherms'] == 51
    assert digest(work / 'run.py') == plan['script_sha256']
    executable = plan['prepared']['executables']['synspec']
    assert digest(plan['prepared']['synspec']) == executable
    spec = plan['specification']
    reference = Path(plan['reference_opacity'])
    assert digest(reference) == plan['reference_opacity_sha256']
    abundance, _ = composition(.7, 0., spec['metals'])
    base = validate_table(reference, abundance, temperatures(spec), sequence(spec['log_density']))
    base_values = np.asarray(base['log_opacity']).reshape(base['shape'])
    hashes = {str(p): digest(p) for p in [work / 'plan.json', work / 'run.py', work / 'final.json', reference]}
    records = []
    for c in plan['compositions']:
        plane = work / c['name']
        saved = json.loads((plane / 'validated.json').read_text())
        assert saved['plan_sha256'] == digest(work / 'plan.json')
        assert digest(plane / 'opacity.bin') == saved['opacity_sha256']
        abundance, _ = composition(c['XH'], c['X3'], spec['metals'])
        merged = validate_table(plane / 'opacity.bin', abundance, temperatures(spec), sequence(spec['log_density']))
        for key in ['shape', 'log_temperature', 'log_density', 'frequency', 'flags', 'ifmol', 'tmolim']:
            assert merged[key] == base[key], key
        values = np.asarray(merged['log_opacity']).reshape(merged['shape'])
        source_cpu = current_cpu = 0.
        reused = 0
        for i in range(len(temperatures(spec))):
            source = plane / f'temperature-{i:03d}'
            receipt = json.loads((source / 'completed.json').read_text())
            assert receipt['input_sha256'] == input_fingerprint(executable, source)
            for name, expected in receipt['outputs'].items():
                assert digest(source / name) == expected
            result = json.loads((source / 'validated.json').read_text())
            assert result['plan_sha256'] == digest(work / 'plan.json')
            original = read_table(source / 'fort.63')
            assert np.array_equal(values[:, :, i], np.asarray(original['log_opacity']).reshape(original['shape'])[:, :, 0])
            assert merged['log_electron_density'][i*len(merged['log_density']):(i+1)*len(merged['log_density'])] == original['log_electron_density']
            source_cpu += result['source_cpu_seconds']
            current_cpu += result['current_call_cpu_seconds']
            reused += int(result['reused_verified_isotherm'])
            for p in [source / 'completed.json', source / 'validated.json']:
                hashes[str(p)] = digest(p)
        assert abs(source_cpu - saved['source_cpu_seconds']) < 1e-8
        # A diagnostic of material sensitivity, not an atmosphere error bound.
        difference = values.astype(np.float64) - base_values
        largest = np.unravel_index(np.argmax(np.abs(difference)), difference.shape)
        records.append(dict(XH=c['XH'], XHe3=c['X3'], plane=str(plane),
                            source_cpu_seconds=source_cpu, current_call_cpu_seconds=current_cpu,
                            reused_isotherms=reused, shape=merged['shape'],
                            max_abs_log_opacity_change=float(np.max(np.abs(difference))),
                            max_change_indices=[int(v) for v in largest],
                            p99_abs_log_opacity_change=float(np.quantile(np.abs(difference), .99))))
        for p in [plane / 'opacity.bin', plane / 'validated.json']:
            hashes[str(p)] = digest(p)
    return dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                all_checks_pass=True, records=records, input_sha256=hashes,
                reviewer_sha256=digest(__file__),
                limitation='Gas opacity composition coverage only. Atmosphere columns, condensation feedback and stellar boundary interpolation require their own checks.')


if __name__ == '__main__':
    result = review(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(all_checks_pass=result['all_checks_pass'], records=result['records']), indent=2))
