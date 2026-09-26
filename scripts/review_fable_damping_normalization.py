#!/usr/bin/env python3
"""Check retained mode damping by dividing oscillator energy by radiated power.

No mode eigenfunction, stellar model or encounter is recalculated. The tensor
power is averaged over phase independently of the worker's damping formula.
"""
import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import shutil


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--scratch', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.scratch.exists():
        raise FileExistsError('Review outputs must be new.')
    root = Path(__file__).resolve().parents[1]
    batch = root / 'docs/research/fable/results/fable-mode-damping-v1'
    ready = json.loads((batch / 'READY.json').read_text())
    for relative, expected in ready['files_sha256'].items():
        assert sha(batch / relative) == expected, relative
    data = json.loads((batch / 'mode_damping_he010_v1.json').read_text())
    benchmark = root / 'docs/research/fable/results/mode_damping_benchmark_v1.json'
    inputs = [Path(__file__), batch / 'READY.json', batch / 'mode_damping.py',
              batch / 'mode_damping_he010_v1.json',
              batch / 'tidal_linear_n15_equal_mass_v1.json', benchmark,
              root / 'docs/research/fable/sph/tidal_linear.py',
              Path('/tmp/ember-fable-lit/rathore2005.pdf')]
    identities = {str(p): sha(p) for p in inputs}
    args.scratch.mkdir()
    for i, path in enumerate(inputs):
        shutil.copyfile(path, args.scratch / f'{i:02d}-{path.name}')

    G, c, yr = 6.67430e-8, 2.99792458e10, 3.15576e7
    M, R = data['constants']['M_g'], data['constants']['R_cm']
    angular = 4 * math.sqrt(math.pi / 5) / 3
    # Reduced quadrupole I = diag(-Izz/2, -Izz/2, Izz).
    # <sum (I'''_ij)^2> = (3/2)*(1/2)*omega^6*Izz_amplitude^2.
    power_prefactor = Fraction(1, 5) * Fraction(3, 2) * Fraction(1, 2)
    lifetime_prefactor = Fraction(1, 2) / power_prefactor
    assert lifetime_prefactor == Fraction(10, 3)
    controls, modes = [], []
    for mode in data['modes']:
        omega, overlap = mode['omega_rad_s'], mode['Q']
        tau = float(lifetime_prefactor) * c**5 / (G * omega**4 * (angular * overlap)**2 * M * R**2)
        for amplitude in [1e-6, 1e-3]:
            for phase in [0.0, 0.37]:
                zz = angular * overlap * amplitude * M * R**2
                energies, powers = [], []
                for j in range(128):
                    theta = 2 * math.pi * j / 128 + phase
                    displacement = amplitude * math.cos(theta)
                    velocity = -omega * amplitude * math.sin(theta)
                    energies.append(0.5 * M * R**2 * (velocity**2 + omega**2 * displacement**2))
                    third = omega**3 * zz * math.sin(theta)
                    powers.append(G / (5 * c**5) * ((-third/2)**2 + (-third/2)**2 + third**2))
                energy = math.fsum(energies) / len(energies)
                power = math.fsum(powers) / len(powers)
                classical_power = G / (60 * c**5) * omega**6 * (3 * zz)**2
                # Same calculation in SI, preserving dimensionless mode amplitude.
                si_energy = 0.5 * (M/1000) * (R/100)**2 * amplitude**2 * omega**2
                si_zz = zz / 1e7
                si_power = 6.67430e-11 / (60 * (c/100)**5) * omega**6 * (3 * si_zz)**2
                differences = dict(phase_average=(energy/power)/tau-1,
                                   classical_quadrupole=power/classical_power-1,
                                   si_units=(si_energy/si_power)/tau-1)
                assert max(map(abs, differences.values())) < 1e-12, differences
                controls.append(dict(order=mode['order'], amplitude=amplitude,
                                     phase=phase, relative_differences=differences))
        modes.append(dict(order=mode['order'], corrected_energy_time_yr=tau/yr,
                          amplitude_time_yr=2*tau/yr,
                          worker_energy_time_yr=mode['tau_GW_energy_yr'],
                          worker_to_energy_time_ratio=mode['tau_GW_energy_s']/tau,
                          worker_matches_energy_time=math.isclose(mode['tau_GW_energy_s'],tau,rel_tol=1e-12)))
    viscous = data['electron_viscosity']['tau_viscous_energy_yr']
    f_time = modes[0]['corrected_energy_time_yr']
    scaled_reference = json.loads(benchmark.read_text())['scaled_to_0p1Msun_yr']
    assert all(sha(p) == value for p, value in identities.items())
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  outcome='worker_energy_damping_time_is_twice_E_over_power',
                  accepted_new_stellar_physics=False,
                  source_url='https://arxiv.org/abs/astro-ph/0406102',
                  source_location='Appendix C, equation C1 and total isolated-mode energy, PDF page 14',
                  oscillator_energy='E = (1/2) a^2 omega^2 M R^2',
                  reduced_tensor_power='P = (3G/20c^5) omega^6 (C Q a M R^2)^2',
                  energy_time_prefactor=str(lifetime_prefactor), worker_prefactor='20/3',
                  controls=controls, controls_pass=True, modes=modes,
                  f_mode_viscous_only_heat_fraction=f_time/(f_time+viscous),
                  corrected_to_worker_scaled_Rathore_reference=f_time/scaled_reference,
                  limitations=['Retained mode eigenfunctions and overlap normalization are not independently recomputed.',
                               'Viscosity is retained as an assumed scenario; this check does not validate it.',
                               'An assigned crystal quality factor does not bound real mode thermalization.',
                               'The source formula identity alone did not validate the worker E/P division.'],
                  input_sha256=identities,
                  preserved_inputs_sha256={str(p):sha(p) for p in args.scratch.iterdir()})
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ['outcome','controls_pass','modes','f_mode_viscous_only_heat_fraction']},indent=2))


if __name__ == '__main__':
    main()
