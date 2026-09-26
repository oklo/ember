"""Archive the deeper member of a reviewed PMS atmosphere pair."""
import datetime
import gzip
import json
from pathlib import Path
import sys

from prepare_nongrey_sources import digest
from generate_nongrey_grid import input_fingerprint


def archive(review_path, destination):
    review_path, destination = Path(review_path), Path(destination)
    review = json.loads(review_path.read_text())
    assert review['all_checks_pass'] and review['accepted_as_gas_source']
    assert digest(Path(__file__).with_name('review_pms_atmosphere_pairs.py')) == review['reviewer_sha256']
    for name, expected in review['input_sha256'].items():
        assert digest(name) == expected, name
    selected = max(review['records'], key=lambda r: r['state']['optical_depth_range'][1])
    source = Path(selected['work'])
    plan = json.loads((source.parent / 'plan.json').read_text())
    case = next(c for c in plan['cases'] if c['name'] == source.name)
    spec = case['specification']
    assert spec['hydrogen'] == [.7] and spec['helium3'] == [0.]
    assert digest(plan['opacity']) == plan['opacity_sha256']
    provenance = dict(prepared=plan['prepared'], specification=spec, mode='gas',
                      XH=.7, X3=0., teff_K=selected['Teff_K'], log_g=selected['log_g'],
                      opacity_sha256=plan['opacity_sha256'], original_work=str(source),
                      plan_sha256=digest(source.parent / 'plan.json'),
                      independent_depth_review=str(review_path),
                      independent_depth_review_sha256=digest(review_path),
                      accepted_as_gas_source=True, selected_for_physical_PMS_track=False)
    destination.mkdir(parents=True, exist_ok=False)
    for name in ['fort.5', 'fort.7', 'fort.8', 'fort.9', 'run.log', 'tas',
                 'ember-masses.dat', 'physics.json', 'fort.15', 'opacity.sha256',
                 'completed.json', 'validated.json']:
        (destination / (name + '.gz')).write_bytes(gzip.compress((source / name).read_bytes(), mtime=0))
    (destination / 'provenance.json.gz').write_bytes(gzip.compress(
        (json.dumps(provenance, indent=2) + '\n').encode(), mtime=0))
    assert input_fingerprint(plan['prepared']['executables']['tlusty'], destination, archived=True) == json.loads(
        (source / 'completed.json').read_text())['input_sha256']
    receipt = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   files_sha256={p.name: digest(p) for p in destination.iterdir()},
                   opacity_source=plan['opacity'], opacity_sha256=plan['opacity_sha256'],
                   contains_opacity_binary=False, archiver_sha256=digest(__file__))
    (destination / 'archive.json').write_text(json.dumps(receipt, indent=2) + '\n')
    return dict(source=str(source), archive=str(destination), Teff_K=selected['Teff_K'],
                log_g=selected['log_g'], chemistry_checked=False)


if __name__ == '__main__':
    print(json.dumps(archive(*sys.argv[1:3]), indent=2))
