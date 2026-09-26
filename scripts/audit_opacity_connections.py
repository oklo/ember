#!/usr/bin/env python3
"""Compare retained opacity sources across the current envelope, without selecting a blend."""
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess

from audit_mesa_oplib_native import read_native, evaluate
from compare_mesa_oplib_means_v2 import stencil

ROOT = Path(__file__).resolve().parents[1]
WORK = Path('/tmp/ember-opacity-connections-v1')
OUTPUT = ROOT/'docs/results/opacity_connections_v1.json'
PROBE = Path('/tmp/ember-opacity-family-probe-build-v1/opacity_family_probe')
SOURCE = Path('/tmp/ember-refractive-hot-family-refined-v1')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    if OUTPUT.exists() or WORK.exists():
        raise FileExistsError('preserve completed or partial connection work')
    WORK.mkdir()
    paths = [Path(__file__), PROBE, ROOT/'scripts/audit_mesa_oplib_native.py',
             ROOT/'scripts/compare_mesa_oplib_means_v2.py',
             ROOT/'docs/results/op_refined_compositions_v3.json',
             ROOT/'docs/results/op_warm_runtime_v1.json',
             ROOT/'docs/results/op_warm_envelope_v1.json',
             ROOT/'docs/reports/2026-09-11/evolution_latest_profile.csv']
    inputs = {str(p):digest(p) for p in paths}
    latest = json.loads(paths[4].read_text())
    original = json.loads(paths[5].read_text())
    envelope = json.loads(paths[6].read_text())
    with paths[7].open() as stream:
        profile = [{k:float(v) for k,v in row.items()} for row in csv.DictReader(stream)]
    selected = [(i,p) for i,p in enumerate(profile) if 1e4 <= p['temperature_K'] <= 1e6]
    queries = [(1,p['X'],p['Y3'],1-p['X']-p['Y3']-p['Y4'],p['temperature_K'],p['density_g_cm3'])
               for i,p in selected]
    manifests = {
        'op_ordinary':Path('/tmp/ember-op-warm-family-v4/density0025/op_gs98_ordinary.dat'),
        'op_absorption':Path('/tmp/ember-op-warm-family-v4/density0025/op_gs98_absorption.dat'),
        'molecular':SOURCE/'aesopus21_gs98_mixture.dat',
        'tops_low':SOURCE/'tops_gs98_mixture_low.dat',
        'tops_refractive_hot':SOURCE/'tops_gs98_mixture_high.dat'}
    retained = {name:{} for name in manifests}
    for row in latest['profile_samples']:
        for name,key in [('op_ordinary','runtime'),('op_absorption','absorption')]:
            retained[name][tuple(row[key]['query'])] = row[key]
    for row in original['profile_samples']:
        for name in ['molecular','tops_low']:
            retained[name][tuple(row[name]['query'])] = row[name]
    replies = {}; counts = {}
    for name,manifest in manifests.items():
        inputs[str(manifest)] = digest(manifest)
        for line in manifest.read_text().splitlines()[1:]:
            p = manifest.parent/line.split('"')[1];inputs[str(p)] = digest(p)
        pending = [q for q in queries if q not in retained[name]]
        text = ''.join(' '.join(str(v) for v in q)+'\n' for q in pending)
        (WORK/(name+'_queries.txt')).write_text(text)
        run = subprocess.run([str(PROBE),str(manifest)],input=text,text=True,capture_output=True,timeout=30)
        (WORK/(name+'_stdout.jsonl')).write_text(run.stdout)
        (WORK/(name+'_stderr.txt')).write_text(run.stderr)
        run.check_returncode()
        rows = [json.loads(line) for line in run.stdout.splitlines()]
        assert len(rows)==len(pending)
        for row,q in zip(rows,pending):
            assert tuple(row['query'])==q
            retained[name][q] = row
        replies[name] = [retained[name][q] for q in queries]
        counts[name] = dict(reused=len(queries)-len(pending),new=len(pending))
    tables = {}
    for entry in envelope['mean_extraction']:
        table = read_native(entry); tables[table['X'],table['Z']] = table
        inputs[entry['path']] = entry['sha256']
    xs = sorted({k[0] for k in tables});zs = sorted({k[1] for k in tables})
    assert len(tables)==16 and len(xs)==len(zs)==4
    c = envelope['composition'];scale = c['atomic_mass_scale']
    results = []
    for index,((zone,p),q) in enumerate(zip(selected,queries)):
        assert abs(p['X']-c['X'])<1e-13 and abs(p['Y3']-c['Y3'])<1e-13
        assert abs(p['Y4']-c['Y4'])<1e-13
        row = dict(zone=zone,query=q,sources={name:values[index] for name,values in replies.items()})
        logt = math.log10(p['temperature_K'])
        logr = math.log10(scale*p['density_g_cm3'])-3*logt+18
        means = []
        for composition_order,thermal_order in [(4,4),(2,4),(4,2)]:
            try:
                ii,wx = stencil(xs,c['atomic_X'],composition_order)
                jj,wz = stencil(zs,c['atomic_Z'],composition_order)
                logk = math.fsum(a*b*evaluate(tables[xs[i],zs[j]],logt,logr,thermal_order)
                                 for i,a in zip(ii,wx) for j,b in zip(jj,wz))
                means.append(dict(composition_order=composition_order,thermal_order=thermal_order,
                                  kappa_baryonic=scale*10**logk))
            except ValueError as error:
                means.append(dict(composition_order=composition_order,thermal_order=thermal_order,
                                  unavailable_reason=str(error)))
        row['mesa_original_means'] = means
        if 'kappa_baryonic' in means[0]:
            reference = means[0]['kappa_baryonic']
            row['relative_to_mesa'] = {name:s['kappa']/reference-1 for name,s in row['sources'].items() if s['covered']}
        results.append(row)
    def coverage(name):
        rows = [r for r in results if r['sources'][name]['covered']]
        return dict(zones=[r['zone'] for r in rows],minimum_T=min((r['query'][4] for r in rows),default=None),
                    maximum_T=max((r['query'][4] for r in rows),default=None))
    gaps = [r for r in results if not any(r['sources'][name]['covered'] for name in
                                        ['molecular','op_ordinary','tops_refractive_hot'])]
    bridge = dict(zones=[r['zone'] for r in gaps],
                  mesa_original_covered=[r['zone'] for r in gaps if 'kappa_baryonic' in r['mesa_original_means'][0]],
                  minimum_T=min((r['query'][4] for r in gaps),default=None),
                  maximum_T=max((r['query'][4] for r in gaps),default=None))
    for p,h in inputs.items():
        assert digest(p)==h,p
    report = dict(outcome='completed_source_connection_comparison',accepted_for_stellar_opacity=False,
                  input_sha256=inputs,reused_runtime_queries=counts,rows=results,
                  source_coverage={name:coverage(name) for name in manifests},
                  gap_without_old_tops_low=bridge,
                  limitations=['Coverage is evaluated only on the saved envelope composition and structure.',
                    'MESA provides ordinary mean opacity; its plasma convention and spectral correction remain unresolved.',
                    'OP, MESA, AESOPUS and TOPS have different atomic/mixture assumptions; differences are not interpolation errors.',
                    'The corrected hot TOPS family is compared directly, without its old temperature blend.',
                    'No new spectrum, atmosphere, source blend, stellar step or selected table is produced.'])
    OUTPUT.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(report=str(OUTPUT),sha256=digest(OUTPUT),layers=len(results),bridge=bridge)))


if __name__ == '__main__':
    main()
