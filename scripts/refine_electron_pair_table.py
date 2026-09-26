#!/usr/bin/env python3
"""Add pair-grid nodes without repeating completed sources; retain new controls."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path

from build_electron_pair_table import calculate, digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('specification', type=Path)
    args = parser.parse_args()
    spec = json.loads(args.specification.read_text())
    work = Path(spec['work'])
    target = Path(spec['report'])
    assert work.is_dir() and not target.exists()
    for path, sha in spec['input_sha256'].items():
        assert digest(path) == sha, path
    old = json.loads(Path(spec['previous_report']).read_text())
    retained = {}
    for entry in old['completed_sources']:
        assert digest(entry['path']) == entry['sha256']
        source = json.loads(Path(entry['path']).read_text())
        retained[source['eta'], source['b_thermal']] = entry
    source_dir = work/'sources'
    source_dir.mkdir()
    complete, controls, tasks, failed = [], [], [], []
    for i, eta in enumerate(spec['eta_grid']):
        for j, b in enumerate(spec['b_thermal_grid']):
            key = f'grid-e{i:02d}-b{j:02d}'
            if (eta, b) in retained:
                complete.append({**retained[eta, b], 'key': key, 'reused': True})
            else:
                tasks.append((key, dict(eta=eta, b_thermal=b, degree=9,
                    orders=(80, 40, 40, 48, 48), tail=50., normalization_tail=80.),
                    str(source_dir/(key+'.json'))))
    heat = json.loads(Path(spec['stellar_heat_report']).read_text())
    control_metadata = {}
    for zone in spec['control_zones']:
        record = heat['records'][zone]
        assert record['zone'] == zone
        for j, case in enumerate(record['cases']):
            eta, b = record['eta_nonrelativistic'], case['b_thermal']
            assert (eta, b) not in retained
            assert not (eta in spec['eta_grid'] and b in spec['b_thermal_grid'])
            key = f'control-zone-{zone:03d}-screen-{j}'
            control_metadata[key] = dict(zone=zone, screening=case['screening'])
            tasks.append((key, dict(eta=eta, b_thermal=b, degree=9,
                orders=(80, 40, 40, 48, 48), tail=50., normalization_tail=80.),
                str(source_dir/(key+'.json'))))
    # Controls enter neither grid construction nor interpolation selection.
    with ProcessPoolExecutor(max_workers=spec['workers']) as pool:
        futures = {pool.submit(calculate, t): t for t in tasks}
        for future in as_completed(futures):
            key = futures[future][0]
            try:
                entry = future.result()
                if key in control_metadata:
                    controls.append({**entry, **control_metadata[key]})
                else:
                    complete.append({**entry, 'reused': False})
                print(json.dumps(dict(completed_new=sum(not s['reused'] for s in complete)+len(controls),
                    total_new=len(tasks), **entry)), flush=True)
            except Exception as error:
                failed.append(dict(key=key, error=repr(error)))
                print(json.dumps(failed[-1]), flush=True)
            progress = dict(completed_sources=complete, independent_controls=controls,
                            failed_sources=failed, total_new=len(tasks))
            temp = work/'progress.json.tmp'
            temp.write_text(json.dumps(progress, indent=2)+'\n')
            temp.replace(work/'progress.json')
    result = dict(outcome='completed_unvalidated_table' if not failed else 'incomplete',
        accepted_for_stellar_evolution=False, eta_grid=spec['eta_grid'],
        b_thermal_grid=spec['b_thermal_grid'], mode_count=10,
        completed_sources=sorted(complete, key=lambda s: s['key']),
        independent_controls=sorted(controls, key=lambda s: s['key']), failed_sources=failed,
        reused_source_count=sum(s['reused'] for s in complete),
        new_source_cpu_seconds=sum(s['cpu_seconds'] for s in complete if not s['reused'])+
                               sum(s['cpu_seconds'] for s in controls),
        specification=str(args.specification), specification_sha256=digest(args.specification),
        input_sha256=spec['input_sha256'],
        limitations=['Prescribed nonrelativistic static Born/Pauli collision model.',
                     'Grid interpolation is not accepted by source generation.',
                     'Independent control matrices must not enter the interpolator.',
                     'Saved hot-star coverage does not establish full future-track coverage.'])
    with target.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(outcome=result['outcome'], grid_sources=len(complete),
        reused=result['reused_source_count'], independent_controls=len(controls),
        new_source_cpu_seconds=result['new_source_cpu_seconds'], report_sha256=digest(target))), flush=True)
    if failed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
