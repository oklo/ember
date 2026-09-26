"""Audit an initial, homogeneous, D-dominated PMS contraction control.

The binding-energy audit retains its original relative threshold and reports
it separately. A second criterion adds only the calculable representation
uncertainty of stored binary64 abundances, and also caps the error against
surface power. Independent D consumption/heat and stoichiometry are required.
This does not select the atmosphere or accept an evolutionary track.
"""
import datetime
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

C = 2.99792458e10
LSUN = 3.828e33
YEAR = 31557600.
ATOMIC = np.array([1.00782503, 3.01602932, 4.00260325, 12.,
                   13.00335484, 14.0030740, 15.9949146, 20., 2.01410177812])
BARYON = np.array([1., 3., 4., 12., 13., 14., 16., 20., 2.])
BINDING = (ATOMIC / BARYON - 1) * C**2
D_HEAT = float((ATOMIC[0] + ATOMIC[8] - ATOMIC[1]) * C**2 / 2)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def checkpoint(path):
    lines = Path(path).read_text().splitlines()
    assert lines[0] == 'EMBER_EVOLUTION_CHECKPOINT 4'
    count = int(lines[8]); offset = 9 + count
    n, mass, age, step, accepted, rejected = lines[offset].split()
    n = int(n)
    data = np.loadtxt(lines[offset+1:offset+1+n])
    assert data.shape == (n, 16) and np.isfinite(data).all()
    assert np.all(data[:, 5:7] == 1), 'baryonic GS98 composition required'
    assert np.all(data[:, 7:] >= 0)
    assert np.max(np.abs(data[:, 7:].sum(axis=1) - 1)) < 2e-12
    return dict(mass=float(mass), age=float(age), data=data)


def review(work):
    work = Path(work).resolve()
    plan = json.loads((work / 'plan.json').read_text())
    receipt = json.loads((work / 'result.json').read_text())
    assert digest(work / 'plan.json') == receipt['plan_sha256']
    for name, expected in plan['input_sha256'].items():
        assert digest(name) == expected, name
    result = json.loads((work / 'stdout.json').read_text())
    assert result == receipt['native_result']
    models = {name: checkpoint(work / (name + '.checkpoint'))
              for name in ['seed', 'full', 'half1', 'final']}
    seed = models['seed']; mass = seed['mass']; m = seed['data'][:, 0]
    assert np.all(np.diff(m) > 0) and m[-1] == mass
    weights = np.diff(np.r_[0., (m[:-1]+m[1:])/2, mass])
    birth = np.loadtxt(plan['command'][10])
    assert np.all(seed['data'][:, 7:] == birth)
    for state in models.values():
        assert state['mass'] == mass and np.array_equal(state['data'][:, 0], m)
        assert np.all(state['data'][:, 10:15] == seed['data'][:, 10:15]), 'inert metals changed'
        # This audit is restricted to a homogeneous initial contraction.
        assert np.all(state['data'][:, 7:] == state['data'][0, 7:])
    balances = []
    pairs = [('seed', 'full', result['years']), ('seed', 'half1', result['years']/2),
             ('half1', 'final', result['years']/2)]
    for j, (start, end, years) in enumerate(pairs):
        dt = years * YEAR
        assert models[end]['age'] - models[start]['age'] == dt
        x0 = models[start]['data'][:, 7:]; x1 = models[end]['data'][:, 7:]
        delta = x1 - x0
        consumed = math.fsum(float(w)*float(-d) for w, d in zip(weights, delta[:, 8]))
        heat = result[end]['deuterium_luminosity'] * dt
        relative_heat = heat / (consumed * D_HEAT) - 1
        # Each rounded endpoint contributes half one unit in the last place.
        # This bounds representation uncertainty, not arbitrary solver error.
        ulps = np.array([[.5*(math.ulp(float(a))+math.ulp(float(b)))
                          for a, b in zip(row0, row1)] for row0, row1 in zip(x0, x1)])
        roundoff_power = math.fsum(float(weights[i]) * math.fsum(
            float(ulps[i, k])*abs(float(BINDING[k])) for k in range(9))
            for i in range(len(m))) / dt
        mass_power = -math.fsum(float(weights[i])*math.fsum(
            float(delta[i, k])*float(BINDING[k]) for k in range(9))
            for i in range(len(m))) / dt
        nuclear_power = result[end]['nuclear_luminosity']
        surface_power = result[end]['luminosity_Lsun'] * LSUN
        # This early control has negligible pp neutrino power; require that
        # all non-D deposited power is below the independent fuel tolerance.
        other_fraction = abs(nuclear_power-result[end]['deuterium_luminosity'])/nuclear_power
        error = mass_power - nuclear_power
        original_relative = error / nuclear_power
        stoich_H = np.abs(delta[:, 0] - .5*delta[:, 8])
        stoich_He3 = np.abs(delta[:, 1] + 1.5*delta[:, 8])
        checks = dict(
            non_D_power_negligible=other_fraction < 2e-6,
            fuel_heat=abs(relative_heat) < 2e-6,
            hydrogen_stoichiometry=bool(np.all(stoich_H <= ulps[:, 0]+.5*ulps[:, 8])),
            helium3_stoichiometry=bool(np.all(stoich_He3 <= ulps[:, 1]+1.5*ulps[:, 8])),
            mass_energy_with_representation=abs(error) <= 2e-6*nuclear_power+roundoff_power,
            mass_energy_surface_budget=abs(error) < 2e-7*surface_power,
            first_law=abs(result['steps'][j]['first_law_error']) < 2e-7,
        )
        balances.append(dict(start=start, end=end, checks=checks,
            D_consumed_g=consumed, D_fuel_heat_relative_error=relative_heat,
            independent_mass_power=mass_power, nuclear_power=nuclear_power,
            original_relative_mass_error=original_relative,
            original_relative_check_pass=abs(original_relative) < 2e-6,
            representation_uncertainty_power=roundoff_power,
            representation_uncertainty_relative_to_nuclear=roundoff_power/nuclear_power,
            error_relative_to_surface=error/surface_power,
            maximum_hydrogen_stoichiometry_error=float(stoich_H.max()),
            maximum_helium3_stoichiometry_error=float(stoich_He3.max()),
            non_D_deposited_fraction=other_fraction))
    checks = dict(completed=receipt['returncode'] == 0,
        coupled_steps=all(s['converged'] for s in result['steps']),
        interval_audits=all(all(b['checks'].values()) for b in balances),
        time_structure=result['maximum_log_structure_difference'] < 1e-4,
        time_composition=result['maximum_abundance_difference'] < 1e-8,
        contracting=result['final']['radius_Rsun'] < result['seed']['radius_Rsun'],
        entropy_decreasing=result['final']['mean_entropy'] < result['seed']['mean_entropy'],
        fully_convective=all(s['convective_mass_fraction'] > .99999 for s in result['steps']))
    files = [work/name for name in ['plan.json', 'result.json', 'stdout.json',
             'seed.checkpoint', 'full.checkpoint', 'half1.checkpoint', 'final.checkpoint']]
    return dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        all_initial_contraction_checks_pass=all(checks.values()), checks=checks,
        balances=balances, seed=result['seed'], final=result['final'],
        original_relative_mass_check_pass=all(b['original_relative_check_pass'] for b in balances),
        input_sha256={str(p):digest(p) for p in files+[Path(__file__).resolve()]},
        selected_for_physical_PMS_track=False,
        limitations=['The atmosphere, source resolution and interpolation require separate acceptance.',
                     'The representation allowance applies only with independent fuel/heat, isotope and total thermal checks; it does not justify inaccurate reaction integration.',
                     'This reviewer requires homogeneous composition, unchanged inert metals and D-dominated nuclear power; it is not a general later-evolution audit.'])


if __name__ == '__main__':
    output = review(sys.argv[1]); path = Path(sys.argv[2])
    serialized = json.dumps(output, indent=2) + '\n'
    with path.open('x') as stream:
        stream.write(serialized)
    print(json.dumps({k:output[k] for k in ['all_initial_contraction_checks_pass', 'checks', 'balances']}, indent=2))
    raise SystemExit(0 if output['all_initial_contraction_checks_pass'] else 1)
