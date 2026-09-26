#!/usr/bin/env python3
"""Bounded integration controls using a candidate whole-star opacity connection.

The saved trajectory is an immutable input. These are independent steps from
that profile, with diffusion omitted and the cooler opacity still unaccepted.
"""
import argparse, csv, hashlib, json, shlex, subprocess, time
from datetime import datetime, timezone
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def table_inputs(path):
    path = Path(path).resolve()
    result = {path}
    with path.open() as f:
        header = f.readline()
        if header.startswith(('EMBER_METAL_HELMHOLTZ ', 'EMBER_OPACITY_MIXTURE ')):
            for line in f:
                for token in shlex.split(line):
                    if token.endswith('.dat'):
                        result.update(table_inputs(path.parent / token))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('probe', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--scratch', type=Path, required=True)
    p.add_argument('--case', action='append', required=True, help='bridge factor:Julian years:abundance tolerance')
    a = p.parse_args()
    assert not a.output.exists() and not a.scratch.exists()
    a.scratch.mkdir()
    root = Path(__file__).resolve().parents[1]
    profile = root / 'docs/reports/2026-09-11/evolution_latest_profile.csv'
    history = root / 'docs/reports/2026-09-11/evolution_latest.csv'
    eos = Path('/tmp/ember-smooth-eos-refined-family-v1/freeeos300_gs98_z020.dat')
    olddir = Path('/tmp/ember-refractive-hot-family-refined-v1')
    tables = [eos, olddir / 'aesopus21_gs98_mixture.dat',
              Path('/tmp/ember-op-warm-family-v4/density0025/op_gs98_ordinary.dat'),
              Path('/tmp/ember-tops-bridge-v2/tables/tops_bridge.dat'),
              olddir / 'tops_gs98_mixture_high.dat',
              root / 'data/conduction/condtab21wd_metals.dat',
              root / 'data/atmosphere/nongrey_gs98_z020_exhaustion_t6400_g600_v1.dat']
    identities = {str(path): digest(path) for path in set().union(*(table_inputs(t) for t in tables))}
    other = [Path(__file__), root / 'scripts/full_star_connection_probe.cpp', a.probe,
             Path('/tmp/ember-physical-core-build-v5/src/libember.a'), profile, history]
    other += list((root / 'src').glob('*.cpp')) + list((root / 'src').glob('*.hpp'))
    other += list((root / 'include/ember').glob('*.hpp'))
    identities.update({str(path.resolve()): digest(path) for path in other})
    assert digest(profile) == '28ab0b571dbc495ff2f48a6c6e7cd79f09e5711f3b725a2f478d4b0390d7ef2f'
    assert digest(history) == 'fcf2cd4251ea3f9e839e7f0c84e3ffd5a412d48373fa140175307c9f541bfd6e'
    assert digest(eos) == '7d6bfa0b26264834307dad651852bc5ff693f91f3861b4934e582dee047fd7f5'
    assert digest(tables[-1]) == '88ee0d061b07c41e51ae8cabd73253ab5e99174c8a71024407598a24138f0345'
    with profile.open() as f:
        rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]
    assert len(rows) == 512
    fields = ['mass_g', 'radius_cm', 'density_g_cm3', 'temperature_K', 'luminosity_erg_s', 'X', 'Y3']
    model = ' '.join(format(row[k], '.17g') for row in rows for k in fields)
    cases = []
    for case in a.case:
        factor, years, tolerance = map(float, case.split(':'))
        cases.append(dict(factor=factor, years=years, seconds=years * 365.25 * 86400, abundance_tolerance=tolerance))
    query = a.scratch / 'queries.txt'
    query.write_text(''.join(f"{c['factor']:.17g} {c['seconds']:.17g} {c['abundance_tolerance']:.17g} 512 {model}\n" for c in cases))
    stdout, stderr = a.scratch / 'responses.jsonl', a.scratch / 'stderr.txt'
    command = [str(a.probe)] + list(map(str, tables))
    input_manifest = a.scratch / 'inputs.json'
    input_manifest.write_text(json.dumps(dict(command=command, input_sha256=identities), indent=2) + '\n')
    begin = time.monotonic()
    code, timeout = None, False
    with query.open() as inp, stdout.open('x') as out, stderr.open('x') as err:
        try:
            code = subprocess.run(command, stdin=inp, stdout=out, stderr=err, timeout=600).returncode
        except subprocess.TimeoutExpired:
            timeout = True
    elapsed = time.monotonic() - begin
    replies = [json.loads(line) for line in stdout.read_text().splitlines()]
    failures = []
    if timeout or code != 0:
        failures.append(dict(process_exit=code, timed_out=timeout))
    preflights = [r for r in replies if r.get('kind') == 'preflight']
    steps = [r for r in replies if r.get('kind') == 'step']
    if len(preflights) != len(cases) or len(steps) != len(cases):
        failures.append('Missing complete preflight or step records; inspect retained replies.')
    failures.extend(r for r in replies if 'error' in r or (r.get('kind') == 'preflight' and r['errors']))
    for r in steps:
        matching = [c for c in cases if (r['factor'],r['dt_seconds'],r['abundance_tolerance']) ==
                    (c['factor'],c['seconds'],c['abundance_tolerance'])]
        assert len(matching) == 1
        c = matching[0]
        if not r['converged']:
            failures.append(dict(case=c, message=r['message']))
            continue
        assert r['factor'] == c['factor'] and r['dt_seconds'] == c['seconds']
        if (not r['input_preserved'] or r['elapsed_age_seconds'] != c['seconds']
            or abs(r['luminosity_balance']) > 2e-8 or abs(r['nuclear_mass_balance']) > 2e-6
            or r['abundance_residual'] > 1e-12):
            failures.append(dict(case=c, message='Preservation, age or conservation check failed.'))
    assert all(digest(name) == value for name, value in identities.items())
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  outcome='passed_conditional_whole_star_controls' if not failures else 'failed_conditional_whole_star_controls',
                  accepted_for_stellar_evolution=False, selected_star_unchanged=True,
                  microscopic_diffusion_enabled=False, points=512, cases=cases,
                  elapsed_seconds=elapsed, failures=failures, preflights=preflights,
                  steps=[{k:v for k,v in r.items() if k != 'model'} for r in steps],
                  input_manifest=dict(path=str(input_manifest), sha256=digest(input_manifest)),
                  artifacts_sha256={str(p):digest(p) for p in a.scratch.iterdir()},
                  limitations=['Cooler bridge density sampling still fails the existing 0.5% acceptance target.',
                               'Independent implicit steps from the saved star, not an accepted trajectory or a timestep accuracy test.',
                               'Smooth native EOS and corrected opacity change physics relative to the saved state.',
                               'No microscopic diffusion; original conduction mixture, Ledoux mixing and plasma-neutrino prescriptions.',
                               'A coherent bridge-opacity multiplier is a sensitivity experiment, not a physical error bound.'])
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k:report[k] for k in ['outcome','elapsed_seconds','failures','preflights','steps']}, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
