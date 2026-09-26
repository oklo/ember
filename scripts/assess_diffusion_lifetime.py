#!/usr/bin/env python3
"""Fuel and timescale estimates from the completed matched diffusion experiment."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
from audit_cn_thermal import weights, YEAR


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--comparison', type=Path, required=True)
    ap.add_argument('--frontier', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    msun, lsun = 1.988409870698051e33, 3.828e33
    q, kb, grav, mp, mu = 6e18, 1.380649e-16, 6.67430e-8, 1.6726219e-24, 1.6605391e-24
    hashes = {}

    def state(path):
        raw = path.read_bytes()
        hashes[str(path)] = hashlib.sha256(raw).hexdigest()
        d = json.loads(raw)
        r = d['model_record']
        m = np.array(r['model'])
        return dict(age_years=d['age_seconds']/YEAR,
                    hydrogen_mass_g=float(weights(m[:, 0]) @ m[:, 5]),
                    luminosity=float(m[-1, 4]), radius=float(m[-1, 1]),
                    nuclear_luminosity=r['nuclear_luminosity'], Teff=r['Teff'],
                    central_hydrogen=float(m[0, 5]), surface_hydrogen=float(m[-1, 5]))

    states = {k: state(a.comparison/p) for k, p in [
        ('initial', 'initial_checkpoint.json'), ('on', 'on/checkpoint.json'),
        ('off', 'off/checkpoint.json')]}
    frontier = state(a.frontier)
    on, off, initial = [states[k] for k in ['on', 'off', 'initial']]
    assert on['age_years'] == off['age_years']
    assert np.isclose(on['age_years']-initial['age_years'], 5e10)
    relative = {k: on[k]/off[k]-1 for k in on if k != 'age_years'}
    burned = {k: initial['hydrogen_mass_g']-states[k]['hydrogen_mass_g'] for k in ['on', 'off']}
    # Representative early radiative center: retained physical diagnostics.
    tc, rho, diffusion, xc = 8.587e6, 325., .6635, .1459
    scales = dict(
        homogeneous_center_separation_years_unit_force=kb*tc/(4*np.pi*grav*rho*mp*diffusion)/YEAR,
        central_burning_years=xc/1.723e-12,
        central_burning_plus_measured_drift_years=xc/(1.723e-12+2.626e-12),
        frontier_hydrogen_Msun=frontier['hydrogen_mass_g']/msun,
        frontier_remaining_fraction_initial_hydrogen=frontier['hydrogen_mass_g']/(.07*msun),
        frontier_constant_nuclear_power_fuel_years=q*frontier['hydrogen_mass_g']/frontier['nuclear_luminosity']/YEAR,
        burning_1em3_Msun_at_1em4_Lsun_years=q*.001*msun/(1e-4*lsun)/YEAR,
        helium_ideal_ion_thermal_energy_at_10MK_erg=1.5*.1*msun/(4*mu)*kb*1e7,
        burning_1em3_Msun_energy_erg=q*.001*msun,
        illustrative_100Gyr_shift_fraction_total_age=1e11/3.6e12)
    report = dict(utc=datetime.now(timezone.utc).isoformat(),
                  outcome='order_of_magnitude_assessment', states=states,
                  on_relative_to_off=relative, consumed_hydrogen_g=burned,
                  relative_change_consumed_hydrogen=burned['on']/burned['off']-1,
                  frontier=frontier, scales=scales, input_sha256=hashes,
                  assumptions=[
                      'The separation estimate is a homogeneous isothermal center with spherical gravity and a force coefficient of order unity. The coefficient D comes from the current collision model; this is not independent rate validation.',
                      'Local depletion times hold rates fixed and are not predicted crossing ages.',
                      'Fuel clocks assume complete H-to-He conversion and constant stated power. They do not show that remaining hydrogen can reach a hot enough burning layer.',
                      'The ideal helium ion thermal energy is a dimensional comparison, not a white dwarf cooling calculation.',
                      'Matched branches share an already radiative initial structure. They isolate the implemented microscopic drift and its associated heat over 50 Gyr, not a full no-diffusion history.'])
    with a.output.open('x') as f:
        json.dump(report, f, indent=2)
        f.write('\n')
    print('Remaining hydrogen difference: %.4g percent; consumed difference: %.4g percent' %
          (relative['hydrogen_mass_g']*100, report['relative_change_consumed_hydrogen']*100))


if __name__ == '__main__':
    main()
