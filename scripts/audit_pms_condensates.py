"""Condensation audit of an accepted PMS gas source, with a stated hot join.

This evaluates chemistry on a fixed atmosphere. It does not claim cloud-opacity
feedback, and rainout is ordered from the deepest layer toward the surface.
"""
import datetime
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

from condensate_chemistry_profile import audit_profile, HOT_JOIN_TEMPERATURE
from generate_nongrey_grid import composition
from import_nongrey_grid import read_text, source_inputs

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_composition(donor, provenance):
    """Verify the recorded composition against the actual atmosphere inputs."""
    x, y = provenance['XH'], provenance['X3']
    spec = provenance['specification']
    inputs = {key: read_text(donor / name) for key, name in [
        ('atmosphere_input', 'fort.5.gz'), ('parameters', 'tas.gz'),
        ('element_masses', 'ember-masses.dat.gz')]}
    if (donor / 'fort.8.gz').exists():
        inputs['initial_structure'] = read_text(donor / 'fort.8.gz')
    source_inputs(inputs, spec, x, y, provenance['teff_K'], provenance['log_g'],
                  read_text(donor / 'run.log.gz'))
    return x, y, composition(x, y, spec['metals'])


def main():
    donor, prepared_path, work = map(Path, sys.argv[1:4])
    work.mkdir(parents=True, exist_ok=False)
    archive = json.loads((donor / 'archive.json').read_text())
    for name, expected in archive['files_sha256'].items():
        assert digest(donor / name) == expected
    provenance = json.loads(read_text(donor / 'provenance.json.gz'))
    assert provenance['accepted_as_gas_source'] and provenance['mode'] == 'gas'
    x, y, (abundances, weights) = source_composition(donor, provenance)
    prepared = json.loads(prepared_path.read_text())
    assert digest(prepared['executable']) == prepared['executable_sha256']
    assert digest(ROOT / 'scripts/fastchem_probe.cpp') == prepared['adapter_sha256']
    for name, expected in prepared['chemistry_data_sha256'].items():
        assert digest(Path(prepared['source']) / 'input/logK' / name) == expected
    profile = []
    for line in read_text(donor / 'run.log.gz').rsplit('FINAL MODEL ATMOSPHERE', 1)[1].splitlines():
        words = line.replace('D', 'E').split()
        if len(words) == 11 and words[0].isdigit():
            profile.append(list(map(float, words)))
    profile.reverse()
    assert len(profile) == provenance['specification']['depths']
    assert all(b[1] < a[1] for a, b in zip(profile, profile[1:]))
    symbols = json.loads((ROOT / 'data/atmosphere/sources/synple-elements.json').read_text())['symbol']
    numbers, masses = {}, {}
    for symbol, abundance, weight in zip(symbols, abundances, weights, strict=True):
        if abundance > 1e-90:
            numbers[symbol.capitalize()] = abundance
            masses[symbol.capitalize()] = weight * 1.67333e-24 / 1.66053906660e-24
    abundance_path = work / 'abundances.dat'
    abundance_path.write_text('# GS98 baryonic element numbers\ne- 0\n' + ''.join(
        f'{s} {12 + math.log10(n):.17g}\n' for s, n in numbers.items()))
    records, hashes = [], {}

    def query(mode, rows, label):
        inputs = ''.join(f'{r[3]:.17g} {r[6]/1e6:.17g}\n' for r in rows)
        (work / (label + '.input')).write_text(inputs)
        result = subprocess.run([prepared['executable'], prepared['source'],
                                 str(abundance_path.resolve()), mode],
                                input=inputs, text=True, capture_output=True, timeout=600)
        (work / (label + '.stderr')).write_text(result.stderr)
        assert result.returncode == 0, (label, result.returncode)
        output = work / (label + '.json.gz')
        output.write_bytes(gzip.compress(result.stdout.encode(), mtime=0))
        hashes[str(output)] = digest(output)
        hashes[str(work / (label + '.input'))] = digest(work / (label + '.input'))
        return json.loads(result.stdout)

    cool = [r for r in profile if r[3] <= HOT_JOIN_TEMPERATURE]
    hot = []
    for row in profile:
        if row[3] > HOT_JOIN_TEMPERATURE:
            joined = row.copy()
            joined[3] = HOT_JOIN_TEMPERATURE
            hot.append(joined)
    hot_join = query('equilibrium', hot, 'hot_join') if hot else None
    for mode in ['equilibrium', 'rainout']:
        chemistry = query(mode, cool, mode)
        audit, condensed = audit_profile(chemistry, profile, numbers, masses, hot_join)
        records.append(dict(mode=mode, **audit,
                            layers_with_condensates=int(condensed.sum())))
        print(json.dumps(records[-1]), flush=True)
    for path in [donor / 'archive.json', prepared_path, Path(__file__),
                 ROOT / 'scripts/audit_condensates.py',
                 ROOT / 'scripts/condensate_chemistry_profile.py',
                 ROOT / 'scripts/generate_nongrey_grid.py',
                 ROOT / 'scripts/import_nongrey_grid.py', abundance_path]:
        hashes[str(path)] = digest(path)
    result = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  all_chemistry_closure_checks_pass=True, records=records,
                  atmosphere_source=str(donor), chemistry=prepared,
                  composition=dict(XH=x, XHe3=y, metals=provenance['specification']['metals']),
                  layers=len(profile), temperature_range_K=[min(r[3] for r in profile), max(r[3] for r in profile)],
                  input_sha256=hashes,
                  limitation='Fixed gas atmosphere: no grain absorption/scattering or gas-depletion feedback. Hotter layers use the declared fully vaporized continuation, checked at 6000 K and each actual pressure.')
    (work / 'review.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
