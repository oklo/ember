#!/usr/bin/env python3
"""Check a fresh stellar history and archive the data plotted in a daily paper."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess

from summarize_remnant_endpoints import summarize


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('track', type=Path)
    p.add_argument('paper', type=Path)
    p.add_argument('--mixing-probe', type=Path, required=True)
    a = p.parse_args()
    receipt_path = a.track.with_suffix('.receipt.json')
    receipt = json.loads(receipt_path.read_text())
    data = json.loads(a.track.read_text())
    assert sha(a.track) == receipt['output_sha256']
    assert not receipt['data_changed_during_run'] and not receipt['restart_changed_during_run']
    assert data['mass_Msun'] == .1 and data['points'] == 512
    assert data['history'][0][0] == 0 and len(data['profile']) == 512
    milestones = summarize(data)
    for row in data['history']:
        for key, value in zip(data['columns'], row, strict=True):
            assert (key == 'last_halfstep_gravothermal_Lsun' if value is None else math.isfinite(value))
    endpoint = dict(zip(data['columns'], data['history'][-1], strict=True))
    for row in data['profile']:
        assert all(math.isfinite(v) for v in row)
        assert abs(sum(row[5:8]) + .02 - 1) < 1e-12
    assert endpoint['central_X'] == data['profile'][0][5]
    assert endpoint['central_T_K'] == data['profile'][0][3]
    assert endpoint['central_rho'] == data['profile'][0][2]
    assert endpoint['surface_X'] == data['profile'][-1][5]
    log_path = a.track.with_suffix('.log')
    domain_errors = [line for line in log_path.read_text().splitlines()
                     if 'CompositionAtmosphereGrid:' in line and 'outside source grid' in line]
    if data['converged']:
        assert receipt['returncode'] == 0
        stop = 'Requested age completed.'
    else:
        assert receipt['returncode'] == 1 and domain_errors
        assert abs(endpoint['Teff_K'] - data['atmosphere_grid_support']['teff_K'][1]) < 1e-4
        stop = 'Accepted history ends at the atmosphere upper temperature; subsequent source-domain rejections and minimum-step failure are retained. Requested age was not reached.'
    profile_text = ''.join(' '.join(format(v, '.17g') for v in row)+'\n' for row in data['profile'])
    eos = data['eos_model'].split(':', 1)[1]
    result = subprocess.run([str(a.mixing_probe), eos, data['opacity_directory']],
                            input=profile_text, text=True, capture_output=True, check=True)
    regions = []
    for line in result.stdout.splitlines():
        region = json.loads(line)
        if regions and regions[-1]['convective'] == region['convective']:
            assert regions[-1]['end_exclusive'] == region['begin']
            regions[-1]['end_exclusive'] = region['end_exclusive']
            regions[-1]['outer_enclosed_fraction'] = region['outer_enclosed_fraction']
            regions[-1]['mass_fraction'] += region['mass_fraction']
        else:
            regions.append(region)
    assert abs(sum(r['mass_fraction'] for r in regions if r['convective'])
               - endpoint['convective_mass_fraction']) < 1e-12
    executable = Path(receipt['command'][0])
    assert sha(executable) == receipt['executable_sha256']
    a.paper.mkdir(parents=True, exist_ok=True)
    artifacts = a.paper / 'artifacts'
    artifacts.mkdir(exist_ok=True)
    manifest_path = a.paper / 'recovery_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    new_entries = []
    files = [(a.track, a.track.name), (receipt_path, receipt_path.name),
             (log_path, log_path.name), (executable, 'ember-evolve-'+sha(executable)[:12])]
    for name, expected in receipt.get('checkpoint_output_sha256', {}).items():
        path = Path(name)
        assert sha(path) == expected
        files.append((path, path.name))
    for path, name in files:
        raw = path.read_bytes()
        archive = artifacts / (name + '.gz')
        packed = gzip.compress(raw, mtime=0)
        if archive.exists():
            assert archive.read_bytes() == packed, 'archive already contains different data'
        else:
            archive.write_bytes(packed)
        entry = {'source': str(path.resolve()), 'archive': str(archive.relative_to(a.paper)),
                 'sha256': hashlib.sha256(raw).hexdigest(),
                 'gzip_sha256': hashlib.sha256(packed).hexdigest(), 'bytes': len(raw)}
        if not any(e['archive'] == entry['archive'] for e in manifest['entries']):
            new_entries.append(entry)
    manifest['entries'].extend(new_entries)
    manifest['scope'] = 'Local recovery files for the plotted fresh trajectory and earlier matched-mesh comparisons. Raw files remain local; physical tables require separate restoration with their recorded hashes.'
    manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
    csv_path = a.paper / 'evolution_latest.csv'
    with csv_path.open('w', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(data['columns'])
        writer.writerows(data['history'])
    profile_path = a.paper / 'evolution_latest_profile.csv'
    with profile_path.open('w', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(data['profile_columns'])
        writer.writerows(data['profile'])
    raw_entry = next(e for e in manifest['entries'] if e['archive'] == 'artifacts/'+a.track.name+'.gz')
    record = {'scope': 'Complete accepted history from one fresh calculation with fixed physical inputs; no history splicing.',
              'input': str(a.track), 'raw_sha256': sha(a.track),
              'raw_archive': raw_entry['archive'], 'raw_archive_sha256': raw_entry['gzip_sha256'],
              'history_csv': csv_path.name, 'history_csv_sha256': sha(csv_path),
              'profile_csv': profile_path.name, 'profile_csv_sha256': sha(profile_path),
              'states': len(data['history']), 'endpoint': endpoint,
              'requested_age_reached': data['converged'], 'stop': stop,
              'first_atmosphere_domain_rejection': domain_errors[0] if domain_errors else None,
              'last_atmosphere_domain_rejection': domain_errors[-1] if domain_errors else None,
              'rejected_attempts': data['rejected_steps'], 'milestones': milestones,
              'mixing_regions': regions, 'mixing_probe_sha256': sha(a.mixing_probe),
              'receipt_sha256': sha(receipt_path), 'executable_sha256': sha(executable),
              'EOS': data['eos_model'], 'atmosphere': data['atmosphere_model'],
              'opacity_directory': data['opacity_directory'],
              'step_error_tolerances': data['step_error_tolerances'],
              'data_changed_during_run': receipt['data_changed_during_run'],
              'source_changed_during_run': receipt['source_changed_during_run'],
              'source_note': 'Changes to the source tree do not alter the copied executable used for this calculation.',
              'cpu_seconds': receipt['child_user_seconds'] + receipt['child_system_seconds'],
              'awake_seconds': receipt['awake_elapsed_seconds'],
              'utc_seconds': receipt['utc_elapsed_seconds'], 'export_script_sha256': sha(__file__)}
    (a.paper / 'evolution_latest_provenance.json').write_text(json.dumps(record, indent=2)+'\n')
    print(f"Archived {len(data['history'])} states through {endpoint['age_yr']/1e12:.4g} trillion yr.")


if __name__ == '__main__':
    main()
