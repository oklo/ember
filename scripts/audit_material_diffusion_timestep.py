#!/usr/bin/env python3
"""Frozen-coefficient linear diffusion time scales in the saved hot interior.

This estimates an explicit method's linear stability limit at fixed T/rho.
It is neither an accuracy criterion nor a bound for changing coefficients,
burning, convection, or the cooler envelope.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import eigh
from native_material_transport import RGAS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    assert not args.output.exists()
    root = Path(__file__).resolve().parents[1];inputs = {}

    def source(p):
        p = Path(p);raw = p.read_bytes();inputs[str(p)] = hashlib.sha256(raw).hexdigest();return raw

    d = json.loads(source(root/'docs/results/stellar_diffusion_pair_transport_v1.json'))
    retained = json.loads(source(root/'docs/results/smooth_eos_refined_retained_sources_v1.json'))
    inp = Path('/tmp/ember-smooth-eos-refined-physical-run-v1/retained/retained.input')
    out = inp.with_suffix('.output')
    for p in (inp, out):
        source(p);assert inputs[str(p)] == retained['artifacts_sha256'][str(p)]
    cache = {}
    for row, reply in zip(inp.read_text().splitlines(), out.read_text().splitlines(), strict=True):
        x, y, T, rho, chemical = map(float, row.split());value = json.loads(reply)
        if chemical == 1 and value['ok']:
            cache[x, y, T, rho] = np.asarray(value['values'][24:28]).reshape(2, 2)/RGAS
    pp = root/'docs/reports/2026-09-11/evolution_latest_profile.csv';source(pp)
    profile = list(csv.DictReader(pp.open()))
    hp = root/'docs/reports/2026-09-11/evolution_latest.csv';source(hp)
    history = list(csv.DictReader(hp.open()))
    source(Path(__file__))
    zones = [r['zone'] for r in d['records'] if not r['convective']]
    assert zones == list(range(len(zones)))
    n = len(zones);mass_grid = np.array([float(r['mass_g']) for r in profile])
    weights = np.zeros(len(profile));weights[0] = mass_grid[0]
    weights[:-1] += .5*np.diff(mass_grid);weights[1:] += .5*np.diff(mass_grid)
    h = np.array([cache[tuple(float(profile[i][k]) for k in ['X','Y3','temperature_K','density_g_cm3'])] for i in zones])
    lam, v = np.linalg.eigh(h/weights[:n, None, None]);assert np.all(lam > 0)
    b = (v*np.sqrt(lam)[:, None, :])@v.swapaxes(-1, -2)
    # The similar symmetric rate matrix is B (face Laplacian K) B,
    # where B_i=sqrt(H_i/m_i). Its eigenvalues are the physical decay rates.
    start = time.process_time();records = [];checks = []
    rng = np.random.default_rng(219)
    for screening in [c['screening'] for c in d['records'][0]['cases']]:
        conductance = []
        for i in range(n):
            matrices = [np.asarray(next(c for c in d['records'][j]['cases'] if c['screening'] == screening)['mobility_scaled_heat'])[:2, :2] for j in (i, i+1)]
            k = RGAS*.5*(matrices[0]+matrices[1]);k = (k+k.T)/2
            rho = .5*(float(profile[i]['density_g_cm3'])+float(profile[i+1]['density_g_cm3']))
            radius = .5*(float(profile[i]['radius_cm'])+float(profile[i+1]['radius_cm']))
            area = 4*np.pi*radius**2
            conductance.append(area**2*rho/(mass_grid[i+1]-mass_grid[i])*k)
        for outer in ('sealed', 'fixed_composition_reservoir'):
            a = np.zeros((2*n, 2*n))
            for i in range(n-1):
                left = slice(2*i, 2*i+2);right = slice(2*i+2, 2*i+4);k = conductance[i]
                a[left, left] += b[i]@k@b[i];a[right, right] += b[i+1]@k@b[i+1]
                a[left, right] -= b[i]@k@b[i+1];a[right, left] -= b[i+1]@k@b[i]
            if outer != 'sealed':
                a[-2:, -2:] += b[-1]@conductance[-1]@b[-1]
            a = (a+a.T)/2
            z = rng.normal(size=(n, 2));potential = np.einsum('nij,nj->ni', b, z)
            face = np.einsum('nij,nj->ni', np.asarray(conductance[:-1]), potential[1:]-potential[:-1])
            divergence = np.zeros((n, 2));divergence[:-1] -= face;divergence[1:] += face
            if outer != 'sealed':divergence[-1] += conductance[-1]@potential[-1]
            action = np.einsum('nij,nj->ni', b, divergence).ravel()
            error = float(np.max(abs(a@z.ravel()-action))/np.max(abs(action)))
            checks.append(dict(name='independent face-divergence matrix action', screening=screening, outer=outer,
                               error=error, tolerance=2e-12, passed=error<2e-12))
            low = float(eigh(a, subset_by_index=[0,0], eigvals_only=True)[0])
            high = float(eigh(a, subset_by_index=[2*n-1,2*n-1], eigvals_only=True)[0])
            checks.append(dict(name='nonnegative diffusion decay rates', screening=screening, outer=outer,
                               error=min(low/high, 0.), tolerance=2e-12, passed=low>=-2e-12*high))
            records.append(dict(screening=screening, outer=outer, maximum_decay_rate_per_second=high,
                                forward_euler_linear_stability_limit_years=2/high/(365.25*86400),
                                minimum_decay_rate_over_maximum=low/high))
    steps = [float(r['step_yr']) for r in history[-500:]]
    large = [v for v in steps if v>=1e6]
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(), outcome='passed' if all(c['passed'] for c in checks) else 'failed',
                  scope=__doc__, zones=zones, records=records, checks=checks, input_sha256=inputs,
                  saved_final_step_years=steps[-1], largest_step_in_last_500_years=max(steps),
                  median_step_of_at_least_one_million_years_in_last_500=float(np.median(large)) if large else None,
                  calculation_cpu_seconds=time.process_time()-start, new_EOS_queries=0, new_collision_integrals=0,
                  accepted_for_timestep_selection=False,
                  limitations=['Temperatures, densities and mobility matrices are held fixed; mobility derivatives are omitted.',
                               'Only the hot nonconvective interior is included. Two outer conditions test sensitivity, not full convective coupling.',
                               'These eigenvalues do not ensure abundance positivity or integration accuracy and omit nuclear reaction stiffness.'])
    with args.output.open('x') as stream:json.dump(result, stream, indent=2, allow_nan=False);stream.write('\n')
    print(json.dumps({k:result[k] for k in ['outcome','records','saved_final_step_years','largest_step_in_last_500_years','median_step_of_at_least_one_million_years_in_last_500']}, indent=2))
    if result['outcome'] != 'passed':raise SystemExit(1)


if __name__ == '__main__':main()
