#!/usr/bin/env python3
"""Check actual atmosphere-extended run segments through an opacity density stop.

Each segment retains its own receipt. Only duplicate restart rows are removed;
no ages or physical values are interpolated. This checks source identity and
the implemented physics, not mesh convergence or physical model accuracy.
"""
import argparse
import datetime
import hashlib
import json
import math
from pathlib import Path
import subprocess
from summarize_remnant_endpoints import summarize

STEP_ONLY = {'step_yr', 'step_error_norm', 'discrete_luminosity_balance',
             'nuclear_rest_mass_balance', 'last_halfstep_coupling_iterations',
             'last_halfstep_gravothermal_Lsun'}
FIXED = ['mass_basis', 'mass_Msun', 'points', 'nuclear_model', 'nuclear_physics',
         'thermal_neutrino_model', 'transport_model', 'eos_model', 'metal_inventory',
         'convection_criterion', 'secular_mixing', 'step_error_tolerances',
         'opacity_directory', 'conduction', 'atmosphere_tau_match']


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def exact_join(old, new, accepted_before):
    if new['columns'] != old['columns'] or new['profile_columns'] != old['profile_columns']:
        raise ValueError('segment columns differ')
    endpoint = dict(zip(old['columns'], old['history'][-1], strict=True))
    start = dict(zip(new['columns'], new['history'][0], strict=True))
    for key, value in start.items():
        expected = (None if key == 'last_halfstep_gravothermal_Lsun' else 0) if key in STEP_ONLY else endpoint[key]
        if value != expected:
            raise ValueError('checkpoint join differs: '+key)
    if (new['restart']['history_start_age_yr'] != endpoint['age_yr'] or
            new['restart']['accepted_steps_before_restart'] != accepted_before):
        raise ValueError('checkpoint age or cumulative accepted count differs')
    for key in FIXED:
        if old[key] != new[key]:
            raise ValueError('physical selection differs: '+key)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('segments', nargs='+', type=Path)
    for name in ['output', 'joined-output', 'mixing-probe', 'physics-probe', 'profile-eos-check']:
        parser.add_argument('--'+name, required=True, type=Path)
    parser.add_argument('--atmosphere-acceptances', required=True, nargs='+', type=Path)
    parser.add_argument('--evidence', required=True, nargs='+', type=Path)
    args = parser.parse_args()
    if len(args.segments) < 2 or len(args.atmosphere_acceptances) != len(args.segments)-1:
        raise ValueError('provide one acceptance for each atmosphere continuation')
    if args.output.exists() or args.joined_output.exists():
        raise FileExistsError('use new check and joined-history paths')
    files = {}

    def remember(path):
        path = Path(path).resolve()
        checksum = digest(path)
        files[str(path)] = checksum
        return checksum

    remember(__file__)
    remember(Path(__file__).with_name('summarize_remnant_endpoints.py'))
    tracks, receipts, physics_files = [], [], {}
    for path in args.segments:
        track = json.loads(path.read_text())
        receipt_path = path.with_suffix('.receipt.json')
        receipt = json.loads(receipt_path.read_text())
        if remember(path) != receipt['output_sha256']:
            raise ValueError('output differs from its own receipt')
        remember(receipt_path)
        remember(path.with_suffix('.log'))
        if receipt['data_changed_during_run'] or receipt['restart_changed_during_run']:
            raise ValueError('physical input changed during evolution')
        if remember(receipt['command'][0]) != receipt['executable_sha256']:
            raise ValueError('executable differs from its receipt')
        for filename, checksum in receipt['checkpoint_output_sha256'].items():
            if remember(filename) != checksum:
                raise ValueError('saved checkpoint changed')
        if track['mass_Msun'] != .1 or track['points'] != 512:
            raise ValueError('expected 512-point 0.1-solar-mass track')
        summarize(track)
        for row in track['history']:
            for key, value in zip(track['columns'], row, strict=True):
                if value is None:
                    if key != 'last_halfstep_gravothermal_Lsun':
                        raise ValueError('missing diagnostic')
                elif not math.isfinite(value):
                    raise ValueError('nonfinite diagnostic')
        if len(track['profile']) != 512:
            raise ValueError('incomplete final structure')
        for row in track['profile']:
            if (len(row) != 8 or not all(math.isfinite(v) for v in row) or
                    min(row[:4]) <= 0 or min(row[5:]) < 0 or abs(sum(row[5:])+.02-1) > 1e-12):
                raise ValueError('invalid profile or composition')
        end = dict(zip(track['columns'], track['history'][-1], strict=True))
        for key, i, j in [('central_X', 0, 5), ('central_Y3', 0, 6),
                          ('central_T_K', 0, 3), ('central_rho', 0, 2), ('surface_X', -1, 5)]:
            if end[key] != track['profile'][i][j]:
                raise ValueError('profile/history mismatch: '+key)
        # Bulk tables are checked but remain outside the paper recovery archive.
        eos = Path(track['eos_model'].split(':', 1)[1]).resolve()
        atmo = Path(track['atmosphere_model'].split(':', 1)[1]).resolve()
        opacity = Path(track['opacity_directory']).resolve()
        selected = {}
        for filename, checksum in receipt['data_sha256'].items():
            f = Path(filename).resolve()
            if f == atmo or f.parent in [eos.parent, opacity, Path('data/conduction').resolve()]:
                if str(f) in physics_files and physics_files[str(f)] != checksum:
                    raise ValueError('shared physical table changed between runs')
                selected[str(f)] = checksum
                physics_files[str(f)] = checksum
        if str(eos) not in selected or str(atmo) not in selected or not any(Path(f).parent == opacity for f in selected):
            raise ValueError('receipt does not cover selected physical data')
        tracks.append(track)
        receipts.append(receipt)
    for path, checksum in physics_files.items():
        if digest(path) != checksum:
            raise ValueError('selected physical data changed: '+path)
    if tracks[0]['history'][0][0] != 0 or 'restart' in tracks[0]:
        raise ValueError('chain must start with a fresh age-zero calculation')
    joined_rows = tracks[0]['history'].copy()
    joins = []
    for i, (old, new) in enumerate(zip(tracks, tracks[1:])):
        exact_join(old, new, len(joined_rows)-1)
        ro, rn = receipts[i:i+2]
        checkpoint = Path(new['restart']['source'])
        if not checkpoint.is_absolute():
            checkpoint = Path(rn['working_directory'])/checkpoint
        checkpoint = checkpoint.resolve()
        old_outputs = {str(Path(k).resolve()): v for k, v in ro['checkpoint_output_sha256'].items()}
        new_inputs = {str(Path(k).resolve()): v for k, v in rn['restart_input_sha256'].items()}
        if old_outputs.get(str(checkpoint)) != digest(checkpoint) or new_inputs.get(str(checkpoint)) != digest(checkpoint):
            raise ValueError('restart is not the preceding saved checkpoint')
        if any('extension' in key and key != 'atmosphere_extension' for key in new):
            raise ValueError('only atmosphere extensions are supported by this check')
        extension = new['atmosphere_extension']
        if (digest(extension['source_executable']) != ro['executable_sha256'] or
                'nongrey:'+extension['source_atmosphere'] != old['atmosphere_model']):
            raise ValueError('extension source differs from preceding calculation')
        acceptance = json.loads(args.atmosphere_acceptances[i].read_text())
        remember(args.atmosphere_acceptances[i])
        if (not acceptance['accepted_for_gas_trajectory'] or
                acceptance['table_sha256'] != digest(new['atmosphere_model'].split(':', 1)[1])):
            raise ValueError('atmosphere extension lacks matching acceptance')
        joins.append({'age_yr': new['restart']['history_start_age_yr'],
                      'accepted_steps_before_restart': len(joined_rows)-1,
                      'checkpoint_sha256': digest(checkpoint), 'physical_join_exact': True})
        joined_rows.extend(new['history'][1:])
    new, receipt = tracks[-1], receipts[-1]
    final = dict(zip(new['columns'], new['history'][-1], strict=True))
    profile = ''.join(' '.join(format(v, '.17g') for v in row)+'\n' for row in new['profile'])
    eos_path = new['eos_model'].split(':', 1)[1]
    mixing = subprocess.run([str(args.mixing_probe), eos_path, new['opacity_directory']],
                            input=profile, text=True, capture_output=True, check=True)
    regions = [json.loads(line) for line in mixing.stdout.splitlines()]
    if abs(sum(r['mass_fraction'] for r in regions if r['convective'])-final['convective_mass_fraction']) > 1e-12:
        raise ValueError('independent mixing classification differs')
    command = [str(args.physics_probe), '--gs98', '--eos-family', eos_path,
               '--opacity-directory', new['opacity_directory']]
    response = subprocess.run(command, input=profile, text=True, capture_output=True, check=True)
    local = [json.loads(line) for line in response.stdout.splitlines()]
    if len(local) != 512 or any(not all(math.isfinite(v) for v in row) for row in local):
        raise ValueError('incomplete independent local-physics response')
    log = args.segments[-1].with_suffix('.log').read_text()
    opacity_errors = [line for line in log.splitlines() if 'TOPS_ATOMIC_GS98_high:' in line and 'outside table' in line]
    atmosphere_errors = [line for line in log.splitlines() if 'CompositionAtmosphereGrid:' in line and 'outside source grid' in line]
    bound = local[0][11]
    if (new['converged'] or receipt['returncode'] != 1 or len(opacity_errors) < 3 or
            atmosphere_errors or not .99*bound < final['central_rho'] < bound or
            final['Teff_K'] >= new['atmosphere_grid_support']['teff_K'][1]):
        raise ValueError('endpoint does not match the expected interior-density stop')
    trial = new['profile'][0].copy()
    trial[2] = bound*(1+1e-6)
    beyond = subprocess.run(command, input=' '.join(format(v, '.17g') for v in trial)+'\n',
                            text=True, capture_output=True)
    if beyond.returncode == 0 or 'TOPS_ATOMIC_GS98_high:' not in beyond.stderr or 'outside table' not in beyond.stderr:
        raise ValueError('native opacity did not reject a query beyond the density bound')
    profile_check = json.loads(args.profile_eos_check.read_text())
    if (not profile_check['profile_passed'] or profile_check['profile_zones'] != 512 or
            profile_check['input_sha256'].get(str(args.segments[-1].resolve())) != digest(args.segments[-1]) or
            profile_check['input_sha256'].get(str(Path(eos_path).resolve())) != digest(eos_path)):
        raise ValueError('EOS profile check does not describe this selected model')
    for path, checksum in profile_check['input_sha256'].items():
        if digest(path) != checksum:
            raise ValueError('EOS check input changed: '+path)
    for path in [args.mixing_probe, args.physics_probe, args.profile_eos_check, *args.evidence]:
        remember(path)
    joined = {**new, 'history': joined_rows, 'history_assembly': {
        'scope': __doc__, 'sources': [str(p.resolve()) for p in args.segments],
        'source_sha256': [digest(p) for p in args.segments], 'physical_join_exact': True}}
    milestones = summarize(joined)
    args.joined_output.write_text(json.dumps(joined, separators=(',', ':'), allow_nan=False)+'\n')
    remember(args.joined_output)
    rejected = tracks[0]['rejected_steps']+sum(t['rejected_steps']-t['restart']['rejected_steps_before_restart'] for t in tracks[1:])
    for path, checksum in files.items():
        if digest(path) != checksum:
            raise ValueError('check input changed')
    report = {'scope': __doc__, 'checked_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'endpoint': final, 'segments': len(tracks), 'joined_states': len(joined_rows),
              'physical_join_exact': True, 'joins': joins,
              'duplicate_initial_diagnostics_excluded': sorted(STEP_ONLY),
              'mixing_regions': regions, 'milestones': milestones, 'requested_age_reached': False,
              'stop': 'Repeated interior opacity density-domain rejections followed by decreasing-step line-search failure.',
              'density_boundary_g_cm3': bound, 'native_beyond_density_rejection': beyond.stderr.strip(),
              'first_opacity_domain_rejection': opacity_errors[0], 'last_opacity_domain_rejection': opacity_errors[-1],
              'opacity_domain_rejection_count': len(opacity_errors),
              'first_atmosphere_domain_rejection': None, 'last_atmosphere_domain_rejection': None,
              'central_conduction_flux_fraction': local[0][9],
              'central_radiative_opacity': local[0][4], 'central_conductive_opacity': local[0][5],
              'all_attempts_including_preceding_terminal_rejections': rejected,
              'profile_EOS_check': str(args.profile_eos_check),
              'CPU_seconds': sum(r['child_user_seconds']+r['child_system_seconds'] for r in receipts),
              'awake_seconds': sum(r['awake_elapsed_seconds'] for r in receipts),
              'physical_inputs_sha256': physics_files, 'input_sha256': files}
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: report[k] for k in ['joined_states', 'endpoint', 'central_conduction_flux_fraction', 'density_boundary_g_cm3']}))


if __name__ == '__main__':
    main()
