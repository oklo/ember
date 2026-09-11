#!/usr/bin/env python3
"""Export a checked checkpoint continuation for a daily working paper.

The plotted history is an assembly of actual accepted states, with each duplicate
restart state removed. Original run receipts remain attached to their separate
segments. The derived history is never assigned a fabricated run receipt.
"""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['joined', 'check', 'paper']:
        p.add_argument(name, type=Path)
    a = p.parse_args()
    check = json.loads(a.check.read_text())
    joined = json.loads(a.joined.read_text())
    files = check['input_sha256']
    if not check['physical_join_exact'] or files[str(a.joined.resolve())] != digest(a.joined):
        raise ValueError('history does not match the checked continuation')
    for path, expected in files.items():
        if digest(path) != expected:
            raise ValueError('checked input changed: '+path)
    assembly = joined['history_assembly']
    sources = [Path(p) for p in assembly['sources']]
    if len(sources) < 2 or not assembly['physical_join_exact']:
        raise ValueError('requires the checked original track and its continuations')
    if [digest(p) for p in sources] != assembly['source_sha256']:
        raise ValueError('assembly source checksum differs')
    tracks = [json.loads(p.read_text()) for p in sources]
    expected_history = tracks[0]['history'] + [row for track in tracks[1:] for row in track['history'][1:]]
    if (joined['history'] != expected_history
            or joined['profile'] != tracks[-1]['profile']
            or len(joined['history']) != check['joined_states']):
        raise ValueError('assembled values differ from accepted source states')
    receipts = [json.loads(p.with_suffix('.receipt.json').read_text()) for p in sources]
    for path, receipt in zip(sources, receipts, strict=True):
        if digest(path) != receipt['output_sha256']:
            raise ValueError('segment receipt does not describe its output')
    a.paper.mkdir(parents=True, exist_ok=True)
    artifacts = a.paper/'artifacts'
    artifacts.mkdir(exist_ok=True)
    manifest_path = a.paper/'recovery_manifest.json'
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {'entries': []}
    entries = []
    for i, path in enumerate([*files, str(a.check.resolve())]):
        path = Path(path)
        raw = path.read_bytes()
        packed = gzip.compress(raw, mtime=0)
        target = artifacts/(f'continuation-{digest(a.joined)[:12]}-{i:03d}-'+path.name+'.gz')
        if target.exists() and target.read_bytes() != packed:
            raise ValueError('archive path already contains different data')
        target.write_bytes(packed)
        entry = {'source': str(path.resolve()), 'archive': str(target.relative_to(a.paper)),
                 'sha256': hashlib.sha256(raw).hexdigest(),
                 'gzip_sha256': digest(target), 'bytes': len(raw)}
        existing = next((e for e in manifest['entries'] if e['archive'] == entry['archive']), None)
        if existing and existing != entry:
            raise ValueError('recovery record changed')
        if not existing:
            manifest['entries'].append(entry)
        entries.append(entry)
    manifest['scope'] = ('Local recovery of the plotted accepted history, its original run segments, '
                         'their receipts and checkpoints, the explicit assembly checks, and separate controls. '
                         'The joined JSON is derived data; physical tables require separate restoration.')
    manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
    for name, columns, rows in [('evolution_latest.csv', joined['columns'], joined['history']),
                                ('evolution_latest_profile.csv', joined['profile_columns'], joined['profile'])]:
        with (a.paper/name).open('w', newline='') as stream:
            writer = csv.writer(stream, lineterminator='\n')
            writer.writerow(columns)
            writer.writerows(rows)
    regions = []
    for r in check['mixing_regions']:
        r = r.copy()
        if regions and regions[-1]['convective'] == r['convective']:
            if regions[-1]['end_exclusive'] != r['begin']:
                raise ValueError('noncontiguous mixing regions')
            regions[-1]['end_exclusive'] = r['end_exclusive']
            regions[-1]['outer_enclosed_fraction'] = r['outer_enclosed_fraction']
            regions[-1]['mass_fraction'] += r['mass_fraction']
        else:
            regions.append(r)
    archive = next(e for e in entries if e['source'] == str(a.joined.resolve()))
    final = receipts[-1]
    record = {'scope': __doc__, 'input': str(a.joined), 'history_is_assembled': True,
              'raw_sha256': digest(a.joined), 'raw_archive': archive['archive'],
              'raw_archive_sha256': archive['gzip_sha256'],
              'serialized_history_note': 'The archived JSON is the explicitly derived joined history. Original segment outputs and receipts are archived separately.',
              'history_csv': 'evolution_latest.csv', 'history_csv_sha256': digest(a.paper/'evolution_latest.csv'),
              'profile_csv': 'evolution_latest_profile.csv', 'profile_csv_sha256': digest(a.paper/'evolution_latest_profile.csv'),
              'states': check['joined_states'], 'endpoint': check['endpoint'],
              'requested_age_reached': check['requested_age_reached'],
              'stop': check.get('stop', 'Atmosphere source temperature limit.'),
              'central_conduction_flux_fraction': check.get('central_conduction_flux_fraction'),
              'first_atmosphere_domain_rejection': check['first_atmosphere_domain_rejection'],
              'last_atmosphere_domain_rejection': check['last_atmosphere_domain_rejection'],
              'rejected_attempts': check['all_attempts_including_preceding_terminal_rejections'],
              'milestones': check['milestones'], 'mixing_regions': regions,
              'EOS': joined['eos_model'], 'atmosphere': joined['atmosphere_model'],
              'opacity_directory': joined['opacity_directory'], 'step_error_tolerances': joined['step_error_tolerances'],
              'executable_sha256': final['executable_sha256'],
              'segment_receipts': [{'source': str(p.with_suffix('.receipt.json')),
                                   'sha256': digest(p.with_suffix('.receipt.json'))} for p in sources],
              'continuation_check': str(a.check), 'continuation_check_sha256': digest(a.check),
              'cpu_seconds': sum(r['child_user_seconds']+r['child_system_seconds'] for r in receipts),
              'awake_seconds': sum(r['awake_elapsed_seconds'] for r in receipts),
              'utc_seconds': sum(r['utc_elapsed_seconds'] for r in receipts),
              'timing_note': 'Sums of the actual segment runs; no atmosphere or source-generation time included.',
              'export_script_sha256': digest(__file__)}
    (a.paper/'evolution_latest_provenance.json').write_text(json.dumps(record, indent=2)+'\n')
    print(f"Exported {record['states']} actual states through {record['endpoint']['age_yr']/1e12:.4g} trillion years.")


if __name__ == '__main__':
    main()
