#!/usr/bin/env python3
"""Assemble a conditional warm/hot connection from retained uncut TOPS means.

No source requests are made. Keep the full native density support allowed by
the retained ratio tables, including regions that fail independent comparisons.
Assembly does not select a stellar opacity or establish physical accuracy.
"""
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path

from assemble_tops_refractive_family_v2 import ratio_corners, ratio_limit
from audit_tops_electron_dispersion import KEV, KB
from audit_tops_group_factor_tables import query
from fetch_tops_composition import digest
from import_tops_composition import read
from reduce_tops_group_factors import add_inputs, verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    if args.work.exists() or args.report.exists():
        raise FileExistsError('preserve previous outputs')
    args.work.mkdir()
    inputs = {}

    def identify(path):
        path = Path(path).resolve()
        add_inputs(inputs, {str(path): digest(path)})
        return path

    for name in ('assemble_tops_bridge.py', 'assemble_tops_refractive_family_v2.py',
                 'audit_tops_factor_hydrogen.py', 'audit_tops_group_factor_tables.py',
                 'audit_tops_electron_dispersion.py', 'fetch_tops_composition.py',
                 'import_tops_composition.py', 'reduce_tops_group_factors.py'):
        identify(Path(__file__).with_name(name))
    plan_path = identify('data/opacity/sources/tops_native_uncut_family_batched_form_v4_specification.json')
    plan = json.loads(plan_path.read_text())
    baseline = identify(plan['baseline_manifest'])
    if digest(baseline) != plan['baseline_manifest_sha256']:
        raise ValueError('baseline changed')
    compositions = {(p['X'], p['Z']): p for p in json.loads(baseline.read_text())['planes']}
    xs = sorted({x for x, z in compositions})
    zs = sorted({z for x, z in compositions})
    if set(compositions) != {(x, z) for x in xs for z in zs}:
        raise ValueError('incomplete native composition grid')
    ratios = {}
    for x in ('000', '020', '050', '075'):
        for z in ('010', '030'):
            path = identify(f'docs/results/tops_cool_ratio_x{x}-z{z}_reduction_v1.json')
            r = json.loads(path.read_text())
            add_inputs(inputs, r['input_sha256'])
            add_inputs(inputs, r['output_sha256'])
            ratios[r['X'], r['Z']] = r
    rxs = sorted({x for x, z in ratios})
    rzs = sorted({z for x, z in ratios})
    ratio_tt = next(iter(ratios.values()))['temperatures_keV']
    if any(r['temperatures_keV'] != ratio_tt for r in ratios.values()):
        raise ValueError('inconsistent ratio temperatures')
    tt = [t for t in ratio_tt if .008 <= t <= .06]
    axes = []
    for plane in compositions.values():
        native_t, native_r, _, _ = read(identify(baseline.parent/plane['file']), plane)
        if any(t not in native_t for t in tt):
            raise ValueError('native temperature absent')
        axes.append(native_r)
    if any(axis != axes[0] for axis in axes):
        raise ValueError('inconsistent native density axes')
    rr = axes[0]
    corners = {key: ratio_corners(*key, rxs, rzs) for key in compositions}
    limits = {(x, z, t): min(ratio_limit(ratios[key], t) for key, w in corners[x, z])
              for x, z in compositions for t in tt}
    jobs = list(plan['requests'])
    for retained in plan['retained_requests']:
        path = identify(retained['plan'])
        if digest(path) != retained['plan_sha256']:
            raise ValueError('retained plan changed')
        jobs.append(json.loads(path.read_text())['requests'][retained['request_index']])
    if len({j['work'] for j in jobs}) != len(jobs):
        raise ValueError('duplicate source jobs')
    selected = [j for j in jobs if any(t in tt and min(j['densities_g_cm3']) <=
                limits[j['X'], j['Z'], t] for t in j['temperatures_keV'])]
    cells, masks, sources = {}, {}, []
    overlap = 0
    for job in selected:
        work = Path(job['work'])
        receipt = json.loads(identify(work/'receipt.json').read_text())
        request = json.loads(identify(work/'request.json').read_text())
        recipe = json.loads(identify(work/'recipe.json').read_text())
        x, z = job['X'], job['Z']
        if (recipe['job'] != job or (receipt['X'], receipt['Z']) != (x, z)
                or digest(work/'request.json') != receipt['request_sha256']
                or request['datype'] != 'gray' or request['plasnu'] != 'off'
                or request['lib'] != 'new'):
            raise ValueError('source request or composition mismatch')
        native_t, native_r, part, excluded = read(identify(work/'source.txt'), receipt,
            dimensions=(len(job['temperatures_keV']), len(job['densities_g_cm3'])))
        if native_t != job['temperatures_keV'] or any(abs(r/q-1) > 5e-5
                for r, q in zip(native_r, job['densities_g_cm3'], strict=True)):
            raise ValueError('source coordinates differ from request')
        used = 0
        for (t, rho), value in part.items():
            if t not in tt or rho not in rr or rho > limits[x, z, t]:
                continue
            key = x, z, t, rho
            flag = (t, rho) in excluded
            if key in cells:
                if masks[key] != flag or (not flag and cells[key] != value):
                    raise ValueError('source overlaps disagree')
                overlap += 1
            else:
                cells[key], masks[key] = value, flag
            used += 1
        add_inputs(inputs, recipe.get('input_sha256', {}))
        sources.append(dict(work=str(work),returned_states=len(part),
                            source_exclusions=len(excluded),used_states=used))
    requests = defaultdict(set)
    prefixes = {}
    omitted = []
    for x, z in sorted(compositions):
        for t in tt:
            prefix = []
            ended = False
            for rho in rr:
                key = x, z, t, rho
                planned = rho <= limits[x, z, t]
                if planned and key not in cells:
                    raise ValueError(f'missing planned source state {key}')
                excluded = planned and masks[key]
                if not planned or excluded:
                    ended = True
                    if excluded:
                        omitted.append(list(key))
                else:
                    if ended:
                        raise ValueError('valid native densities do not form a prefix')
                    prefix.append(rho)
                    for corner, weight in corners[x, z]:
                        requests[corner].add((t, rho))
            if len(prefix) < 4:
                raise ValueError('insufficient density support')
            prefixes[x, z, t] = prefix
    probe = Path('/tmp/ember-tops-spectral-runtime-probe-v1')
    reference = json.loads(identify('docs/results/tops_contraction_spectral_runtime_adaptive_v1.json').read_text())
    add_inputs(inputs, reference['input_sha256'])
    if inputs.get(str(probe.resolve())) != digest(probe):
        raise ValueError('unverified ratio runtime')
    # The retained binary was built before this shared header changed. Verify
    # its original bytes from a preserved source tree, never against a new hash.
    original_header = str(Path('include/ember/opacity_blend.hpp').resolve())
    retained_header = Path('/tmp/ember-screening-reuse-v2/source/include/ember/opacity_blend.hpp')
    relocations = []
    if original_header in inputs and digest(Path(original_header)) != inputs[original_header]:
        expected = inputs[original_header]
        if digest(retained_header) != expected:
            raise ValueError('original runtime header bytes unavailable')
        relocations.append(dict(original=original_header,retained=str(retained_header),sha256=expected))
        del inputs[original_header]
        add_inputs(inputs,{str(retained_header):expected})
    verify(inputs)
    factor = {}
    for key, points in sorted(requests.items()):
        points = sorted(points)
        table = next(Path(p) for p in ratios[key]['output_sha256'] if Path(p).name == 'factor.dat')
        replies = query(probe, table, points, key[0], args.work, f'ratio-x{key[0]}-z{key[1]}')
        if any(row[0] < 1-1e-12 for row in replies):
            raise ValueError('ratio below unity')
        factor[key] = {point: row[0] for point, row in zip(points, replies, strict=True)}
    output = args.work/'tables'
    output.mkdir()
    manifest = [f'EMBER_OPACITY_MIXTURE 1 {len(zs)} logRho ATOMIC_native_uncut_times_refractive_ratio']
    planes = []
    states = 0
    with (args.work/'assembled_states.jsonl').open('x') as stream:
        for z in zs:
            lines = [f'EMBER_OPACITY_TABLE 2 {len(xs)} {len(tt)} {len(rr)} Conditional TOPS bridge; native uncut mean times refractive ratio',
                     ' '.join(format(math.log10(r), '.17g') for r in rr),
                     ' '.join(format(math.log10(t*KEV/KB), '.17g') for t in tt)]
            for x in xs:
                lines.append(f'{x:.17g} {z:.17g}')
                counts = []
                for t in tt:
                    values = []
                    for rho in prefixes[x, z, t]:
                        logf = math.fsum(w*math.log(factor[c][t, rho]) for c,w in corners[x,z])
                        grey = cells[x,z,t,rho]
                        kappa = grey*math.exp(logf)
                        if not math.isfinite(kappa) or kappa <= 0:
                            raise ValueError('invalid assembled opacity')
                        values.append(math.log10(kappa))
                        stream.write(json.dumps(dict(X=x,Z=z,temperature_keV=t,density_atomic_g_cm3=rho,
                            uncut_rosseland=grey,refractive_ratio=math.exp(logf),kappa=kappa))+'\n')
                        states += 1
                    counts.append(len(values))
                    lines.append(str(len(values))+' '+' '.join(format(v,'.17g') for v in values))
                planes.append(dict(X=x,Z=z,density_prefix_sizes=counts))
            name = f'tops_bridge_z{round(z*1000):03d}.dat'
            (output/name).write_text('\n'.join(lines)+'\n')
            manifest.append(f'{z:.17g} "{name}"')
    (output/'tops_bridge.dat').write_text('\n'.join(manifest)+'\n')
    verify(inputs)
    result = dict(scope=__doc__,outcome='assembled_candidate',accepted_for_stellar_opacity=False,
        X=xs,Z=zs,temperatures_keV=tt,densities_atomic_g_cm3=rr,states=states,planes=planes,
        grey_source_jobs_used=len(selected),new_source_requests=0,exact_overlapping_states=overlap,
        unique_ratio_queries=sum(map(len,requests.values())),sources=sources,source_exclusions=omitted,
        ratio_source_recovery=[dict(X=x,Z=z,passed=r['uncut_recovery_check_passed'],
            maximum_relative_error=r['maximum_uncut_recovery_relative_error']) for (x,z),r in ratios.items()],
        limitations=['Independent grey/ratio interpolation and source joins require separate comparisons.',
            'Full supported densities are retained, including domains with known ratio-composition errors.',
            'Normalization to native means is not an independent validation of spectral transport.',
            'Absorption/scattering separation is not supplied; no atmosphere table is replaced.'],
        source_relocations=relocations,input_sha256=inputs,
        output_sha256={str(p.resolve()):digest(p) for p in args.work.rglob('*') if p.is_file()})
    args.report.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('outcome','states','grey_source_jobs_used','new_source_requests',
                                          'exact_overlapping_states','unique_ratio_queries')}),flush=True)


if __name__ == '__main__':
    main()
