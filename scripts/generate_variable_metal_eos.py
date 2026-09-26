#!/usr/bin/env python3
"""Build physical FreeEOS planes in metal fraction and normalized H/He ratios.

Every source plane uses nonnegative mass fractions. Identical pure-hydrogen
endpoints share source work. Source generation is restartable; completed
family receipts are immutable. This does not select the family for evolution.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import time
from write_scientific_result import write_result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('probe', type=Path)
    p.add_argument('work', type=Path)
    p.add_argument('report', type=Path)
    p.add_argument('--grid-from', type=Path, required=True)
    p.add_argument('--metals', type=float, nargs='+', default=[0, .02, .05])
    p.add_argument('--hydrogen-share', type=float, nargs='+', default=[0, .01, .1, .3, .6, .85, .95, 1])
    p.add_argument('--helium3-share', type=float, nargs='+', default=[0, .1, .2])
    p.add_argument('--jobs', type=int, default=6)
    a = p.parse_args()
    if a.report.exists(): raise FileExistsError(a.report)
    for axis, low, high in [(a.metals, 3, 4), (a.hydrogen_share, 4, 100), (a.helium3_share, 3, 4)]:
        if (not low <= len(axis) <= high or sorted(set(axis)) != axis
                or not all(math.isfinite(v) and 0 <= v <= 1 for v in axis)):
            raise ValueError('invalid composition axis')
    if a.metals[-1] >= 1 or not 1 <= a.jobs <= 8: raise ValueError('invalid family or worker count')
    root = Path(__file__).resolve().parents[1]
    sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
    sources = [Path(__file__), root/'scripts/generate_metal_eos.py',
               root/'scripts/metal_eos_composition.py', root/'scripts/import_freeeos_potential.py',
               root/'scripts/eos_source_coverage.py', root/'scripts/stellar_composition.py',
               root/'scripts/generate_nongrey_grid.py', a.probe, a.grid_from,
               root/'data/atmosphere/sources/synple-elements.json',
               root/'data/opacity/sources/tops_gs98_x070_z020.request.json']
    hashes = {str(path.resolve()): sha(path) for path in sources}
    grid = json.loads(a.grid_from.read_text())
    step = grid['logT'][1]-grid['logT'][0]
    if min(len(grid[k]) for k in ['logT', 'logQ']) < 9:
        raise ValueError('need at least nine source coordinates to retain five potential nodes')
    a.work.mkdir(parents=True, exist_ok=True)
    plan = dict(metals=a.metals, hydrogen_share=a.hydrogen_share, helium3_share=a.helium3_share,
                input_sha256=hashes, grid=grid, step=step)
    plan_path = a.work/'plan.json'
    if plan_path.exists() and json.loads(plan_path.read_text()) != plan:
        raise ValueError('changed source settings require a new work directory')
    plan_path.write_text(json.dumps(plan, indent=2)+'\n')
    archive = a.work/'scripts'; archive.mkdir(exist_ok=True)
    for path in sources:
        if path.suffix == '.py': shutil.copy2(path, archive/path.name)
    jobs = []
    for iz, z in enumerate(a.metals):
        for iu, u in enumerate(a.hydrogen_share):
            x = (1-z)*u
            ys = [(1-z)*(1-u)*v for v in a.helium3_share] if u < 1 else [0.]
            jobs.append((iz, iu, x, ys, z, a.work/f'z{iz:02d}-u{iu:03d}'))
    started = time.monotonic()

    def run(job):
        iz, iu, x, ys, z, work = job
        work.mkdir(exist_ok=True)
        command = [sys.executable, str(root/'scripts/generate_metal_eos.py'),
                   str(a.probe.resolve()), str(work), '--hydrogen', repr(x),
                   '--helium3', *map(repr, ys), '--metallicity', repr(z),
                   '--step', repr(step), '--jobs', '1', '--electron-integrals', 'numerical',
                   '--grid-from', str(a.grid_from.resolve())]
        with (work/'generate.log').open('w') as out:
            subprocess.run(command, stdout=out, stderr=subprocess.STDOUT, check=True)
        planes = []
        for iv in range(len(a.helium3_share)):
            source_index = iv if len(ys) > 1 else 0
            raw = work/f'plane-{source_index:03d}/source.json.gz'
            potential = work/f'potential-{source_index:03d}.dat'
            if not potential.exists():
                with (work/f'import-{source_index:03d}.log').open('w') as out:
                    subprocess.run([sys.executable, str(root/'scripts/import_freeeos_potential.py'),
                                    str(raw), str(potential)], stdout=out, stderr=subprocess.STDOUT, check=True)
            data = json.loads(gzip.decompress(raw.read_bytes()))
            lines = potential.read_text().splitlines(); start = lines.index('data')+1
            planes.append(dict(index=[iz, iu, iv], composition=data['composition'],
                               raw=str(raw.resolve()), raw_sha256=sha(raw),
                               potential=str(potential.resolve()), potential_sha256=sha(potential),
                               raw_states=sum(r is not None for r in data['data']),
                               absent_source_states=sum(r is None for r in data['data']),
                               failed_source_states=sum(r is not None and r[0] != 0 for r in data['data']),
                               potential_nodes=len(lines)-start,
                               masked_nodes=sum(line.split()[0] == '0' for line in lines[start:])))
        print(f'Completed Z={z:.4g}, H share={a.hydrogen_share[iu]:.4g}', flush=True)
        return planes

    with ThreadPoolExecutor(a.jobs) as pool:
        planes = list(itertools.chain.from_iterable(pool.map(run, jobs)))
    for path, expected in hashes.items():
        if sha(path) != expected: raise ValueError(f'source changed during generation: {path}')
    manifest = a.work/'variable_metal.dat'
    lines = ['EMBER_VARIABLE_METAL_HELMHOLTZ 1']
    for key, axis in [('metals', a.metals), ('hydrogen_share', a.hydrogen_share), ('helium3_share', a.helium3_share)]:
        lines.append(f'{key} {len(axis)} '+' '.join(format(x, '.17g') for x in axis))
    lines.extend(json.dumps(str(Path(plane['potential']).relative_to(a.work.resolve()))) for plane in planes)
    manifest.write_text('\n'.join(lines)+'\n')
    unique = {plane['raw']: plane for plane in planes}
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(), outcome='generated_unselected_candidate',
                  selected_for_evolution=False, elapsed_seconds=time.monotonic()-started,
                  manifest=str(manifest.resolve()), manifest_sha256=sha(manifest), input_sha256=hashes,
                  planes=planes, unique_source_planes=len(unique), requested_planes=len(planes),
                  source_states=sum(plane['raw_states'] for plane in unique.values()),
                  absent_source_states=sum(plane['absent_source_states'] for plane in unique.values()),
                  failed_source_states=sum(plane['failed_source_states'] for plane in unique.values()),
                  masked_nodes=sum(plane['masked_nodes'] for plane in unique.values()),
                  limitations=['Source family only: interpolation and composition derivatives require separate controls.',
                               'Fixed relative GS98 metal abundances; no independent element pattern response.',
                               'No stellar evolution has used this family.'])
    write_result(a.report, report)
    print(json.dumps({k: v for k, v in report.items() if k not in ['planes', 'input_sha256']}, indent=2))


if __name__ == '__main__': main()
