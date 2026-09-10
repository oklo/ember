#!/usr/bin/env python3
"""Compare recorded F77 milestones with the published 0.1-solar-mass example."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--paper', required=True, type=Path)
    parser.add_argument('--f77-source', required=True, type=Path)
    parser.add_argument('--history', type=Path,
                        default=root/'docs/reports/2026-09-10/f77_history.csv')
    parser.add_argument('--output', type=Path,
                        default=root/'docs/results/f77_lba97_comparison_v1.json')
    args = parser.parse_args()
    provenance_path = args.history.with_name('f77_history_provenance.json')
    provenance = json.loads(provenance_path.read_text())
    if sha(args.history) != provenance['csv_sha256']:
        raise ValueError('history differs from its extraction record')
    onset_path = root/'docs/results/f77_core_onset_comparison_v1.json'
    onset = json.loads(onset_path.read_text())
    expected = next(v for k, v in onset['input_files_sha256'].items()
                    if k.endswith('/baseline.f'))
    if sha(args.f77_source) != expected:
        raise ValueError('source differs from the recorded F77 calculation')
    with args.history.open() as stream:
        records = [{k: v if k == 'run' else int(v) if k in
                    ['model', 'central_convective_flag'] else float(v)
                    for k, v in row.items()} for row in csv.DictReader(stream)]
    published = {
        'reference': 'Laughlin, Bodenheimer and Adams (1997), ApJ 482, 420-432',
        'doi': 'https://doi.org/10.1086/304125',
        'source_pdf_sha256': sha(args.paper),
        'extraction': 'Numbers transcribed from section 3.2, pages 422-424, and visually checked Figure 1; no curve digitization.',
        'main_sequence_start': {'age_yr': 2e9, 'Teff_K': 2228,
                                'log10_L_Lsun': -3.38, 'page': 422},
        'peak_helium3': {'age_yr': 1.38e12, 'central_Y3': .0995, 'page': 423,
                        'note': 'Text gives 9.95 percent; Figure 1 labels 9.96 percent.'},
        'radiative_core': {'age_yr': 5.742e12, 'central_X': .16,
                           'Teff_K': 3450, 'log10_L_Lsun': -2.54,
                           'page': 423, 'age_source': 'Figure 1: 5742 Gyr',
                           'unit_note': 'Several prose ages on page 423 are printed as Gyr. Figure 1 and the surrounding discussion establish trillion-year evolution.'},
        'maximum_Teff_K': 5807,
        'cooling_endpoint': {'Teff_K': 1651, 'log10_L_Lsun': -5.287,
                             'time_since_core_yr': 540e9, 'page': 424},
        'nuclear_lifetime': 'Somewhat more than 6 trillion years; no exact threshold supplied in this comparison.',
    }
    runs = []
    for name in ['ajr-20000', 'ferguson-20000']:
        rows = [row for row in records if row['run'] == name]
        audit = next(run for run in onset['runs'] if run['name'] == name)
        if len(rows) != audit['recorded_models']:
            raise ValueError('accepted-state count differs')
        core = next(row for row in rows if row['central_convective_flag'] == 0)
        if core['model'] != audit['central_flag_transitions'][0]['after']['model']:
            raise ValueError('central transition differs from the independent audit')
        peak = max(rows, key=lambda row: row['central_Y3'])
        runs.append({
            'name': name,
            'early_luminosity_minimum': min((row for row in rows if row['age_yr'] < 1e10),
                                           key=lambda row: row['L_Lsun']),
            'early_minimum_note': 'A directly selected luminosity minimum, not the published definition that nuclear heating first supplies all surface luminosity.',
            'peak_helium3': peak,
            'first_central_nonconvective_flag': core,
            'maximum_Teff': max(rows, key=lambda row: row['Teff_K']),
            'last_accepted': rows[-1],
            'returncode': audit['receipt']['returncode'],
            'core_age_relative_difference': core['age_yr']/published['radiative_core']['age_yr']-1,
            'helium3_peak_age_relative_difference': peak['age_yr']/published['peak_helium3']['age_yr']-1,
        })
    report = {
        'scope': __doc__, 'published': published, 'current_runs': runs,
        'conclusion': 'The broad evolutionary sequence agrees, but the recorded reconstruction does not reproduce the published ages or composition at central radiative-core formation. A similar peak surface temperature is insufficient for quantitative validation.',
        'source_differences': 'The current OPACV uses the Cox-Tabor King IVa table with BFGH composition scaling, rather than the published Weiss-Keady-Magee interpolation between King IVa and Ross-Aller 2. Its nuclear rates use CF88 rather than the cited Bahcall rates. The Ferguson low-temperature option is a further change. These are documented differences, not a measured attribution of the discrepancy.',
        'limits': 'No new F77 runs, recovered historical executable, or mesh/timestep convergence study. F77 central flags have hysteresis; original brackets remain in the onset audit. Surface quantities retain the printed precision. The comparison uses the CSV solar-unit convention, without assuming the original paper used identical solar constants.',
        'input_sha256': {str(path): sha(path) for path in
                         [args.history, provenance_path, onset_path, args.f77_source]},
        'script_sha256': sha(Path(__file__)),
    }
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    for run in runs:
        core = run['first_central_nonconvective_flag']
        print(f"{run['name']}: core age {core['age_yr']/1e12:.4g} trillion yr, "
              f"XH={core['central_X']:.4g}")


if __name__ == '__main__':
    main()
