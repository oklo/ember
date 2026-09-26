"""Check the bounded stellar coupling control, independently of its native audits.

The artificial D injection tests integration; it is not a physical PMS start.
Usage: python review_pms_deuterium_coupling.py WORK_DIR OUTPUT_JSON
"""
import datetime
import hashlib
import json
from pathlib import Path
import sys


def review(work):
    work = Path(work)
    plan = json.loads((work / 'plan.json').read_text())
    result = json.loads((work / 'result.json').read_text())
    report = json.loads((work / 'report.json').read_text())
    digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    assert digest(work / 'native') == plan['native_sha256']
    # Atomic mass differences, conserved baryonic D mass, and c in CGS.
    # The factor two is D's baryon count, not its measured atomic mass.
    energy_per_D_gram = (1.00782503 + 2.01410177812 - 3.01602932) * 2.99792458e10**2 / 2
    duration = result['years'] * 31557600.
    balances = []
    for start, end, dt in [('seed', 'full', duration), ('seed', 'half1', duration / 2),
                           ('half1', 'final', duration / 2)]:
        fuel_energy = (result[start]['deuterium_g'] - result[end]['deuterium_g']) * energy_per_D_gram
        deposited = dt * result[end]['deuterium_luminosity']
        balances.append(dict(start=start, end=end, fuel_energy_erg=fuel_energy,
                             deposited_energy_erg=deposited, relative_difference=deposited/fuel_energy - 1))
    checks = dict(
        completed=report.get('returncode') == 0,
        coupled_steps_converged=all(s['converged'] for s in result['steps']),
        first_law=all(abs(s['first_law_error']) < 2e-7 for s in result['steps']),
        nuclear_mass=all(abs(s['nuclear_mass_error']) < 2e-6 for s in result['steps']),
        time_structure=result['maximum_log_structure_difference'] < 1e-4,
        time_composition=result['maximum_abundance_difference'] < 1e-8,
        deuterium_heat=all(abs(b['relative_difference']) < 2e-6 for b in balances),
        contracting=result['final']['radius_Rsun'] < result['seed']['radius_Rsun'],
        fully_convective=all(s['convective_mass_fraction'] > .99999 for s in result['steps']),
        positive_remaining_deuterium=0 < result['final']['deuterium_g'] < result['seed']['deuterium_g'],
    )
    return dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                all_checks_pass=all(checks.values()), checks=checks, deuterium_energy_balances=balances,
                scope=plan['scope'], selected_for_physical_track=False, result=result, completion=report,
                input_sha256={str(work/name):digest(work/name) for name in
                              ['plan.json', 'result.json', 'report.json', 'native']},
                reviewer_sha256=digest(__file__), limitations=[
                    'A small D inventory is injected into a relaxed late-PMS numerical control. This is not a Hayashi-start model.',
                    'The implicit million-year interval checks fuel/heat coupling, not the physical duration or luminosity profile of the artificially introduced burning transient.',
                    'The physical initial D/He3 inventory, low-gravity boundary grid, finite transport and continuous late atmosphere remain separate requirements.',
                ])


if __name__ == '__main__':
    output = review(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps(dict(all_checks_pass=output['all_checks_pass'], checks=output['checks'],
                         deuterium_energy_balances=output['deuterium_energy_balances']), indent=2))
    raise SystemExit(0 if output['all_checks_pass'] else 1)
