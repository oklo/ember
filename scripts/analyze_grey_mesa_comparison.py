#!/usr/bin/env python3
"""Freeze a grey/no-diffusion comparison at equal remaining hydrogen mass.

This reads running Ember output without changing its inputs or checkpoint.
It deliberately distinguishes Ember solver CPU from MESA process CPU.
"""
import argparse
import csv
import hashlib
import json
import gzip
import zlib
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import numpy as np


def check_convective_exchange(record, selected=False):
    """Allow only the reviewed neutral-face convection in a no-diffusion track."""
    fluxes = record['cn_boundary_fluxes']
    assert max(record['secular_diffusivity'], default=0) == 0
    if not fluxes:
        return
    assert selected, 'boundary exchange requires the reviewed convection treatment'
    face = record['neutral_face']
    assert len(fluxes) == 1 and fluxes[0][0] == face
    assert np.isfinite(fluxes[0][1]).all()
    assert record['mixing_parameter'] >= 0 and record['mixing_rate'] >= 0
    if record['mixing_parameter'] == 0:
        assert all(value == 0 for value in fluxes[0][1])
    else:
        assert 0 <= face < len(record['model'])-1
        assert abs(record['boundary_stability']) < 1e-9
        assert record['mlt_required_radiative_excess']/record['boundary_grad_ad'] < 1e-8
        assert any((lo == face+1 or hi == face+1) and hi-lo > 1
                   for lo, hi in record['mixing_regions'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ember', type=Path, required=True)
    parser.add_argument('--mesa-benchmark', type=Path, required=True)
    parser.add_argument('--mesa-track', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--segment', type=Path, action='append', default=[],
                        help='Retained continuation segments in chronological order')
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False, parents=True)
    receipts = {}

    def freeze(path, name):
        content = path.read_bytes()
        (args.output / name).write_bytes(content)
        receipts[str(path.resolve())] = dict(
            snapshot=name, sha256=hashlib.sha256(content).hexdigest())
        return content

    checkpoint = json.loads(freeze(args.ember/'checkpoint.json', 'checkpoint.json'))
    age = checkpoint['numerical_history_row']['age_years']
    rows = [json.loads(line) for line in (args.ember/'history.jsonl').read_text().splitlines()]
    rows = [row for row in rows if row['age_years'] <= age]
    assert rows[-1] == checkpoint['numerical_history_row']
    assert all(a['age_years'] < b['age_years'] for a, b in zip(rows, rows[1:]))
    assert all(row['error_norm'] <= 1 for row in rows)
    history_bytes = ''.join(json.dumps(row)+'\n' for row in rows).encode()
    (args.output/'ember_history.jsonl').write_bytes(history_bytes)
    receipts[str((args.ember/'history.jsonl').resolve())] = dict(
        snapshot='ember_history.jsonl', scope='Prefix through the frozen checkpoint',
        sha256=hashlib.sha256(history_bytes).hexdigest())
    plan = json.loads(freeze(args.ember/'plan.json', 'ember_plan.json'))
    mesa_receipt = json.loads(freeze(args.mesa_benchmark/'completion.json', 'mesa_completion.json'))
    raw = freeze(args.mesa_benchmark/'LOGS/history.data', 'mesa_benchmark_history.data')
    columns = raw.decode().splitlines()[5].split()
    history = np.loadtxt(args.output/'mesa_benchmark_history.data', skiprows=6)
    mesa = dict(zip(columns, history[-1]))
    assert mesa_receipt['returncode'] == 0
    log = freeze(args.mesa_benchmark/'run.log', 'mesa_run.log').decode()
    assert 'xa_central_lower_limit' in log
    assert np.all(history[:, columns.index('num_retries')] == 0)
    freeze(args.mesa_benchmark/'inlist_project', 'mesa_inlist_project')
    freeze(args.mesa_benchmark/'comparison.net', 'mesa_comparison.net')
    freeze(args.mesa_track, 'mesa_complete_track.csv')
    with (args.output/'mesa_complete_track.csv').open() as stream:
        mesa_track = [{key: float(value) for key, value in row.items() if value}
                      for row in csv.DictReader(stream)]

    # Each code's reported solar mass convention is retained. The tiny
    # difference in adopted M_sun is well below this diagnostic's accuracy.
    target = mesa['total_mass_h1']/.1
    index = next(i for i, row in enumerate(rows) if row['global_species'][0] <= target)
    assert index > 0
    left, right = rows[index-1:index+1]
    assert abs(left['mixed_mass_fraction']-1) < 1e-12
    assert abs(right['mixed_mass_fraction']-1) < 1e-12
    assert abs(mesa['mass_conv_core']/.1-1) < 1e-12
    weight = ((left['global_species'][0]-target) /
              (left['global_species'][0]-right['global_species'][0]))
    assert 0 <= weight <= 1
    ember = {key: (1-weight)*left[key]+weight*right[key] for key in
             ['age_years', 'Teff', 'luminosity', 'radius', 'central_T',
              'central_density', 'solver_cpu_seconds']}
    ember.update(hydrogen_mass_Msun=.1*target,
                 helium3_mass_Msun=.1*((1-weight)*left['global_species'][1]+
                                      weight*right['global_species'][1]),
                 luminosity_Lsun=ember['luminosity']/3.828e33,
                 radius_Rsun=ember['radius']/6.957e10)
    mesa_point = dict(age_years=mesa['star_age'], Teff=10**mesa['log_Teff'],
                      luminosity_Lsun=10**mesa['log_L'], radius_Rsun=10**mesa['log_R'],
                      hydrogen_mass_Msun=mesa['total_mass_h1'],
                      helium3_mass_Msun=mesa['total_mass_he3'],
                      central_T=10**mesa['log_cntr_T'], central_density=10**mesa['log_cntr_Rho'],
                      process_cpu_seconds=mesa_receipt['cpu_seconds'],
                      wall_seconds=mesa_receipt['wall_seconds'],
                      saved_models=len(history), zones=int(mesa['num_zones']))
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(),
        comparison='Grey Eddington atmosphere, no microscopic diffusion, 0.1 solar mass, X=0.7, Z=0.02',
        ember=ember, mesa=mesa_point, ember_snapshot_endpoint=rows[-1],
        ember_interpolation=dict(method='Linear in total hydrogen between two accepted states',
                                 weight=weight, bracketing_rows=[left, right]),
        ember_solver_to_mesa_process_cpu_ratio=ember['solver_cpu_seconds']/mesa_receipt['cpu_seconds'],
        ember_to_mesa_age_ratio=ember['age_years']/mesa['star_age'],
        ember_to_mesa_luminosity_ratio=ember['luminosity_Lsun']/mesa_point['luminosity_Lsun'],
        source_preparation_wall_seconds=plan['source_preparation_wall_seconds'],
        validation=dict(all_snapshot_time_errors_pass=True,
                        both_matched_states_fully_convective=True,
                        mesa_benchmark_retries=0),
        limitations=[
            'Ember has its own EOS, opacity, pp+CN network and static starting model. MESA uses supplied microphysics and a start obtained through contraction. This does not isolate atmosphere or diffusion as the cause of their difference.',
            'Equal remaining hydrogen is an operational comparison stage, not equal composition or equal released energy. The helium-3 inventories differ.',
            'Ember values, including CPU, are interpolated between the stated accepted models; no new model at exactly the target fuel is claimed.',
            'Ember timing covers native solver work, including required rejected attempts on the selected continuation. MESA timing is total native process CPU, including model loading and output. These scopes are close but not identical.',
            'Separate validation reruns, source-table construction, compilation, human/agent analysis, idle time and MESA pre-main-sequence construction are excluded from the efficiency comparison.',
            'The source preparation wall time is measured; its CPU time was not recorded. It must not be converted into an invented CPU count.',
            'These are early fully convective timings, not full-lifetime or cold-remnant benchmarks. The Ember history is a snapshot of a continuing diagnostic.',
            'Both codes retain their own solar-mass normalization; these differ by about 3e-5 relative.'
        ], input_receipts=receipts)
    (args.output/'report.json').write_text(json.dumps(report, indent=2)+'\n')

    if args.segment:
        limits = dict(region_species_error=2e-11, global_species_error=2e-11,
                      catalyst_number_error=2e-9, mass_energy_error=2e-6,
                      first_law_error=2e-8)
        maxima = dict.fromkeys(limits, 0.)
        attempts, profiles = [], []
        convection_selections = []
        for segment in args.segment:
            segment_plan = json.loads((segment/'plan.json').read_text())
            selection_path = segment_plan.get('boundary_selection')
            selected_convection = False
            if selection_path:
                selection_path = Path(selection_path)
                source = selection_path.read_bytes()
                selection = json.loads(source)
                digest = hashlib.sha256(source).hexdigest()
                assert segment_plan['inputs_sha256'][str(selection_path)] == digest
                assert selection['selected_for_grey_evolution']
                for path in [segment_plan['native_command'][0],
                             str(Path(__file__).with_name('grey_neutral_step.py').resolve()),
                             str(Path(__file__).with_name('evolve_grey_neutral_track.py').resolve())]:
                    assert segment_plan['inputs_sha256'][path] == selection['input_sha256'][path]
                    assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == selection['input_sha256'][path]
                selected_convection = True
                convection_selections.append(dict(path=str(selection_path), sha256=digest))
            end = min(age, json.loads((segment/'checkpoint.json').read_text())[
                'age_seconds']/(365.25*86400))
            local_attempts = [json.loads(line) for line in
                              (segment/'attempts.jsonl').read_text().splitlines()]
            selected = [a for a in local_attempts if a['accepted'] and
                        a['start_years']+a['interval_years'] <= end+1e-3]
            for attempt in selected:
                assert attempt['error_norm'] <= 1
                assert len(attempt['conservation']) == 3
                for audit in attempt['conservation']:
                    for key, limit in limits.items():
                        assert np.isfinite(audit[key]) and abs(audit[key]) <= limit
                        maxima[key] = max(maxima[key], abs(audit[key]))
            attempts.extend(selected)
            # A live gzip stream lacks its footer. Consume only complete
            # records, bounded by the already frozen accepted checkpoint.
            decoder = zlib.decompressobj(16+zlib.MAX_WBITS)
            data = decoder.decompress((segment/'accepted_models.jsonl.gz').read_bytes())
            for line in data[:data.rfind(b'\n')+1].splitlines():
                profile = json.loads(line)
                if profile['age_seconds']/(365.25*86400) > end+1e-3:
                    continue
                record = profile['model_record']
                model = np.asarray(record['model'])
                assert record['converged'] and np.isfinite(model).all()
                assert (model[:, 1:4] > 0).all() and (np.diff(model[:, 0]) > 0).all()
                check_convective_exchange(record, selected_convection)
                profiles.append(profile)
        assert len(profiles) == len(attempts) == len(rows)-1
        for profile, row in zip(profiles, rows[1:]):
            assert abs(profile['age_seconds']/(365.25*86400)-row['age_years']) < 1e-3
            assert profile['model_record']['Teff'] == row['Teff']
        payload = ''.join(json.dumps(profile)+'\n' for profile in profiles).encode()
        (args.output/'accepted_models.jsonl.gz').write_bytes(gzip.compress(payload, mtime=0))
        (args.output/'accepted_attempts.jsonl').write_text(
            ''.join(json.dumps(attempt)+'\n' for attempt in attempts))
        report['validation'].update(accepted_intervals_checked=len(attempts),
            retained_profiles_checked=len(profiles), conservation_maxima=maxima,
            conservation_limits=limits, all_retained_profiles_finite_positive_and_no_diffusion=True)
        report['profile_segments'] = [str(path.resolve()) for path in args.segment]
        report['convective_boundary_selections'] = convection_selections
        (args.output/'report.json').write_text(json.dumps(report, indent=2)+'\n')

    timing_segments = []
    for segment in args.segment:
        completion_path = segment/'completion.json'
        if not completion_path.is_file():
            break
        completion = json.loads(completion_path.read_text())
        if 'total_native_cpu_seconds' not in completion:
            break
        timing_segments.append(dict(path=str(completion_path.resolve()),
            sha256=hashlib.sha256(completion_path.read_bytes()).hexdigest(),
            process_cpu_seconds=completion['total_native_cpu_seconds']))
    if timing_segments and len(timing_segments) == len(args.segment):
        report['ember_total_computation'] = dict(
            process_cpu_seconds=sum(item['process_cpu_seconds'] for item in timing_segments),
            segments=timing_segments,
            scope='Total native process CPU including startup, rejected and terminal failed steps; table preparation and independent controls excluded.')
        endpoint_time = f"{report['ember_total_computation']['process_cpu_seconds']/60:.4g} CPU-min total"
    else:
        endpoint_time = f"{rows[-1]['solver_cpu_seconds']/60:.4g} CPU-min (solver)"
    (args.output/'report.json').write_text(json.dumps(report, indent=2)+'\n')

    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42})
    fig, (hr, fuel) = plt.subplots(1, 2, figsize=(7.6, 4.2))
    fig.subplots_adjust(left=.105, right=.985, bottom=.16, top=.83, wspace=.36)
    orange, blue = '#a94d2e', '#326b91'
    hr.plot([10**r['log_Teff'] for r in mesa_track], [r['log_L'] for r in mesa_track],
            color=blue, lw=1.4, alpha=.75)
    hr.plot([r['Teff'] for r in rows], np.log10([r['luminosity']/3.828e33 for r in rows]),
            color=orange, lw=1.8)
    fuel.plot([r['star_age']/1e12 for r in mesa_track], [r['total_mass_h1'] for r in mesa_track],
              color=blue, lw=1.4, alpha=.75)
    fuel.plot([r['age_years']/1e12 for r in rows], [.1*r['global_species'][0] for r in rows],
              color=orange, lw=1.8)
    for point, color in [(ember, orange), (mesa_point, blue)]:
        hr.plot(point['Teff'], np.log10(point['luminosity_Lsun']), 'o', ms=5,
                markerfacecolor='white', markeredgecolor=color)
        fuel.plot(point['age_years']/1e12, point['hydrogen_mass_Msun'], 'o', ms=5,
                  markerfacecolor='white', markeredgecolor=color)
    hr.plot(rows[-1]['Teff'], np.log10(rows[-1]['luminosity']/3.828e33), '.', color=orange, ms=6)
    fuel.plot(rows[-1]['age_years']/1e12, .1*rows[-1]['global_species'][0], '.', color=orange, ms=6)
    fuel.annotate(f"{rows[-1]['age_years']/1e12:.4g} Tyr\n{endpoint_time}", xy=(
        rows[-1]['age_years']/1e12, .1*rows[-1]['global_species'][0]),
        xytext=(.97, .68), textcoords='axes fraction', ha='right', va='bottom',
        arrowprops=dict(arrowstyle='-', color=orange, alpha=.5, lw=.7),
        fontsize=8, color=orange)
    fuel.text(.05, .12, 'At the open circles:\n'
              f"Ember {ember['solver_cpu_seconds']/60:.4g} CPU-min\n"
              f"MESA {mesa_receipt['cpu_seconds']/60:.4g} CPU-min",
              transform=fuel.transAxes, fontsize=8)
    hr.set(xlim=(7200, 2500), ylim=(-3.15, -1.85),
           xlabel='Effective surface temperature (K)', ylabel=r'$\log_{10}(L/L_\odot)$')
    fuel.set(xlim=(0, max(2.7, rows[-1]['age_years']/1e12*1.04)), ylim=(0, .072),
             xlabel='Age (Tyr)', ylabel=r'Remaining hydrogen mass ($M_\odot$)')
    for ax in (hr, fuel):
        ax.grid(alpha=.13)
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_formatter(FuncFormatter(lambda value, _: f'{value:.4g}'))
    fig.legend(handles=[Line2D([], [], color=orange, lw=1.8, label='Ember: grey, no diffusion'),
                        Line2D([], [], color=blue, lw=1.4, alpha=.75, label='MESA: grey, no diffusion')],
               loc='upper center', ncol=2, frameon=False, fontsize=9)
    fig.savefig(args.output/'grey_mesa_comparison.pdf', metadata={'CreationDate': None})
    fig.savefig(args.output/'grey_mesa_comparison.png', dpi=180)
    plt.close(fig)
    print(json.dumps({key: report[key] for key in ['ember', 'mesa',
          'ember_solver_to_mesa_process_cpu_ratio']}, indent=2))


if __name__ == '__main__':
    main()
