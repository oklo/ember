#!/usr/bin/env python3
"""Independent conservation and timestep comparisons for coupled CN stellar steps."""
import argparse
import hashlib
import json
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

YEAR = 365.25 * 86400
C0, N0 = .02 * .171836 / 12, .02 * .050335 / 14
ATOMIC = np.array([1.00782503, 3.01602932, 4.00260325, 12., 13.00335484, 14.003074])
MASS = np.array([1., 3., 4., 12., 13., 14.])
LIGHT = 2.99792458e10


def physical(model):
    cn = model[:, 7:10] * MASS[3:]
    extra = cn.sum(axis=1) - (12*C0 + 14*N0)
    return np.column_stack((model[:, 5:7], .98-model[:, 5:7].sum(axis=1)-extra, cn))


def weights(m):
    w = np.zeros(len(m)); w[0] = m[0]
    dm = np.diff(m); w[:-1] += dm/2; w[1:] += dm/2
    return w


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--native', type=Path, required=True)
    ap.add_argument('--scratch', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--secular', action='store_true', help='enable alpha_sc=0.1 and alpha_th=1')
    ap.add_argument('--microscopic', action='store_true', help='enable H/He diffusion with explicitly selected trace CN advection')
    ap.add_argument('--abundance-tolerance', type=float, default=1e-13)
    ap.add_argument('--case', action='append', help='run only named controls; repeat for several')
    a = ap.parse_args(); a.scratch.mkdir()
    if a.output.exists(): raise ValueError('preserve previous results')
    start = time.monotonic(); identities = {}
    def pin(path):
        path = Path(path); digest = hashlib.sha256()
        with path.open('rb') as f:
            for chunk in iter(lambda: f.read(8*1024*1024), b''): digest.update(chunk)
        identities[str(path)] = digest.hexdigest()
        return path
    tables = [
        '/tmp/ember-smooth-eos-refined-family-v1/freeeos300_gs98_z020.dat',
        '/tmp/ember-refractive-hot-family-refined-v1/aesopus21_gs98_mixture.dat',
        '/tmp/ember-op-warm-family-v4/density0025/op_gs98_ordinary.dat',
        '/tmp/ember-tops-bridge-v2/tables/tops_bridge.dat',
        '/tmp/ember-refractive-hot-family-refined-v1/tops_gs98_mixture_high.dat',
        'data/conduction/condtab21wd_metals.dat',
        'data/atmosphere/nongrey_gs98_z020_exhaustion_t6400_g600_v1.dat',
        '/tmp/ember-refractive-hot-family-refined-v1/tops_gs98_mixture_low.dat']
    if a.microscopic:
        tables[6] = '/tmp/ember-fable-nongrey-x030-candidate-v1.dat'
        tables += ['.1', '1', '/tmp/ember-collision-transport-table-v1.dat']
    for p in [a.native, __file__, 'scripts/cn_thermal_probe.cpp', *[t for t in tables if '/' in t]]: pin(p)
    cases = [
        ('2500gyr', 'out/evolution-metal-512-2500gyr-gas-checkpoint.json', 0., 1e8, (1., 1.), 1),
        ('3400gyr', 'out/evolution-cold-remnant-forward-512-3400gyr-v1.json', .5, 1e7, (1., 1.), 1),
        ('800myr', '/tmp/ember-adaptive-diffusion-v11/checkpoint.json', 1., 2.5e5, (1., 1.), 0),
        # Controlled paired perturbation from the verified half-conversion 3200 K atmosphere.
        ('3400gyr_boundary', 'out/evolution-cold-remnant-forward-512-3400gyr-v1.json', .5, 1e7,
         (4992.034853/4913.3726575, 10255895.29/10873599.27), 1),
        ('800myr_boundary', '/tmp/ember-adaptive-diffusion-v11/checkpoint.json', 1., 2.5e5,
         (9938.914075872983/9939.47464515195, 1487364.0713018973/1475825.603305754), 0),
    ]
    # The maximum X=.25 held-out error, applied as a paired sensitivity to
    # this nearby stellar structure. It is not an evolved X=.25 model.
    spacing = json.loads(pin('docs/research/fable/transport/atmosphere/x025_heldout_interpolation_v1.json').read_text())
    sample = spacing['rows'][0]
    cases.append(('800myr_xspacing', '/tmp/ember-adaptive-diffusion-v11/checkpoint.json', 1., 2.5e5,
                  (1/(1+sample['T']['loglinear_error']), 1/(1+sample['Pgas']['loglinear_error'])), 0))
    if a.microscopic:
        cases.append(('800myr_zero_cn_drift', '/tmp/ember-adaptive-diffusion-v11/checkpoint.json', 1., 2.5e5, (1., 1.), 0))
    queries = (a.scratch/'queries.txt').open('w')
    raw = (a.scratch/'responses.jsonl').open('w')
    stderr = (a.scratch/'native.stderr').open('w')
    if a.case:
        if set(a.case)-{c[0] for c in cases}: raise ValueError('unknown case')
        cases = [c for c in cases if c[0] in a.case]
    process = subprocess.Popen([str(a.native), *tables, *(['.1', '1'] if a.secular and not a.microscopic else [])], stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE, stderr=stderr, text=True)
    (a.scratch/'native_pid.txt').write_text(str(process.pid)+'\n')
    summaries, failures = [], []
    def step(old, dt, factors, warm_choice, label):
        values = [warm_choice, *factors, dt, a.abundance_tolerance, len(old),
                  *([int('zero_cn_drift' not in label), 1e-6] if a.microscopic else []), *old.flat]
        query = ' '.join(format(v, '.17g') for v in values)+'\n'
        queries.write(query); queries.flush(); process.stdin.write(query); process.stdin.flush()
        line = process.stdout.readline(); raw.write(line); raw.flush()
        if not line: raise RuntimeError('native process ended without reply')
        result = json.loads(line)
        (a.scratch/(label+'.json')).write_text(json.dumps(result)+'\n')
        if not result.get('converged'): raise RuntimeError(result.get('message', result.get('error')))
        new = np.array(result['model']); source = np.array(result['physical_sources'])
        energy = np.array(result['energy_cells']); w = weights(old[:, 0]); mass = old[-1, 0]
        assert np.isfinite(new).all() and np.min(physical(new)) >= 0
        assert np.array_equal(new[:, 0], old[:, 0])
        assert result['input_preserved'] and result['elapsed_age_seconds'] == dt
        catalyst_error = float(abs(w@(new[:, 7:10].sum(axis=1)-old[:, 7:10].sum(axis=1))/mass/(C0+N0)))
        assert catalyst_error < 2e-9, ('total catalyst number', catalyst_error)
        if not a.microscopic: assert abs(new[:, 7:10].sum(axis=1)/(C0+N0)-1).max() < 2e-12
        # Uniform inert material is not explicitly represented in this six-species audit.
        normalization = physical(new).sum(axis=1)+.02-(12*C0+14*N0)
        assert abs(normalization-1).max() < 2e-14
        delta = physical(new)-physical(old)
        residual = delta-dt*source
        diffusivity = np.array(result.get('secular_diffusivity', np.zeros(len(old)-1)))
        radius = (new[:-1, 1]+new[1:, 1])/2
        density = (new[:-1, 2]+new[1:, 2])/2
        g = dt/mass*(4*np.pi*radius**2*density)**2*diffusivity/np.diff(new[:, 0])
        flux = g[:, None]*(physical(new)[:-1]-physical(new)[1:])
        if a.microscopic:
            flux[:] = 0
            for face, rate in result['cn_boundary_fluxes']:
                flux[face] = dt/mass*np.array([rate[0], rate[1], -sum(rate), *rate[2:]])
        region_error = 0.
        for begin, end in result['mixing_regions']:
            balance = w[begin:end]@residual[begin:end]/mass
            if begin: balance -= flux[begin-1]
            if end < len(new): balance += flux[end-1]
            region_error = max(region_error, float(abs(balance).max()))
            if end > begin+1: assert np.all(new[begin:end, 5:10] == new[begin, 5:10])
        assert region_error < 2e-12, ('region balance', region_error)
        global_error = float(abs(w@residual/mass).max())
        assert global_error < 1e-13, ('global balance', global_error)
        # Direct physical isotope masses, without the native lookup binding correction.
        rest_release = -float(w@delta@(ATOMIC/MASS-1))*LIGHT**2/dt
        nuclear = float(w@(energy[:, 3]+energy[:, 4]))
        rest_error = rest_release/nuclear-1
        assert abs(rest_error) < 2e-6, ('nuclear mass-energy', rest_error)
        grav = -(energy[:, 0]-energy[:, 2]+energy[:, 1]*(1/new[:, 2]-1/old[:, 2]))/dt
        luminosity_error = float(w@(energy[:, 3]-energy[:, 5]+grav))/new[-1, 4]-1
        assert abs(luminosity_error) < 2e-8, ('first law', luminosity_error)
        summary = {k: v for k, v in result.items() if k not in ('model', 'physical_sources', 'energy_cells', 'mixing_regions', 'secular_diffusivity', 'cn_boundary_fluxes', 'total_species_rates')}
        summary.update(maximum_region_species_residual=region_error, global_species_residual=global_error,
                       relative_total_catalyst_number_error=catalyst_error,
                       nonzero_secular_faces=int(np.count_nonzero(diffusivity)), maximum_secular_diffusivity=float(diffusivity.max()),
                       independently_recomputed_nuclear_mass_balance=rest_error,
                       independently_recomputed_luminosity_balance=luminosity_error)
        return new, summary
    try:
        for name, path, fraction, years, factors, warm_choice in cases:
            try:
                data = json.loads(pin(path).read_text())
                m = np.array(data['profile'])[:, :7] if 'profile' in data else np.array(data['model_record']['model'])
                old = np.column_stack((m, np.tile([C0*(1-fraction), 0, N0+C0*fraction], (len(m), 1))))
                coarse, c1 = step(old, years*YEAR, factors, warm_choice, name+'-full')
                half, c2 = step(old, years*YEAR/2, factors, warm_choice, name+'-half1')
                fine, c3 = step(half, years*YEAR/2, factors, warm_choice, name+'-half2')
                differences = dict(log_structure=float(abs(np.log(coarse[:, 1:4]/fine[:, 1:4])).max()),
                                   physical_abundances=float(abs(physical(coarse)-physical(fine)).max()),
                                   surface_luminosity=float(abs(coarse[-1, 4]/fine[-1, 4]-1)))
                summaries.append(dict(case=name, interval_years=years, imposed_initial_carbon_conversion=fraction,
                                      warm_opacity='retained_TOPS_low' if warm_choice else 'OP',
                                      boundary_temperature_factor=factors[0], boundary_gas_pressure_factor=factors[1],
                                      steps=[c1, c2, c3], full_vs_two_halves=differences,
                                      fine_surface_R_L_Teff=[fine[-1, 1], fine[-1, 4], c3['Teff']],
                                      fine_central_T_H_He3_CN=fine[0, [3, 5, 6, 7, 8, 9]].tolist()))
                print(name, 'three coupled steps pass; full/half', differences, flush=True)
            except Exception as e:
                failures.append({'case': name, 'error': str(e)})
                print(name, 'FAILED', str(e), flush=True)
    finally:
        process.stdin.close(); code = process.wait(timeout=120)
        queries.close(); raw.close(); stderr.close()
    for p in a.scratch.iterdir():
        if p.is_file(): pin(p)
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  outcome='passed_coupled_cn_controls' if not failures and code == 0 else 'failed',
                  elapsed_seconds=time.monotonic()-start, native_exit_code=code, cases=summaries,
                  failures=failures, input_sha256=identities, limitations=[
                      'Imposed catalyst abundances on saved structures; no physical restart history or trajectory advance is inferred.',
                      'Fixed-GS98 EOS/opacity and bulk screening/pp He4 lookup; C/N material response is not incorporated.',
                      ('H/He microscopic diffusion and secular mixing enabled; CN shares the helium-group velocity, except the explicitly named zero-CN-drift control. Separate catalyst settling, CN material feedback and oxygen branches are omitted.' if a.microscopic else
                       'Secular mixing uses alpha_sc=0.1, alpha_th=1; microscopic drift and oxygen branches absent.' if a.secular
                       else 'No microscopic or secular species transport in these controls; no oxygen branches or catalyst settling.'),
                      'Earlier controls explicitly select retained low-temperature TOPS where the OP table lacks dense-envelope coverage; no source fallback or extrapolation.',
                      'Paired atmosphere perturbation is constant; one measured source node does not establish a global correction.',
                      'Full-versus-half differences are reported; conservation alone does not establish production timestep accuracy.'])
    a.output.write_text(json.dumps(report, indent=2)+'\n')
    return 1 if failures or code else 0


if __name__ == '__main__': raise SystemExit(main())
