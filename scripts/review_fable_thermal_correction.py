#!/usr/bin/env python3
"""Review the sealed Fable correction without changing or rerunning its survey."""
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
BATCH = ROOT/'docs/research/fable/results/fable-batch-3-thermal-correction-v1'
COMMON = ROOT/'docs/research/fable/sph/fable_common.py'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ready = json.loads((BATCH/'READY.json').read_text())
    inputs = {str(BATCH/'READY.json'): sha(BATCH/'READY.json'),
              str(COMMON): sha(COMMON), str(Path(__file__)): sha(Path(__file__))}
    for name, expected in ready['files_sha256'].items():
        path = BATCH/name
        assert sha(path) == expected, name
        inputs[str(path)] = expected
    sys.path.insert(0, str(COMMON.parent))
    spec = importlib.util.spec_from_file_location('fable_sealed_thermal', BATCH/'thermal_response_v2.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    r = module.Remnant()
    # Independent SI expression and analytic low-temperature Debye/Sommerfeld energy.
    kb, hbar, mu, me = 1.380649e-23, 1.054571817e-34, 1.66053906660e-27, 9.1093837015e-31
    charge, eps0, mass = 1.602176634e-19, 8.8541878128e-12, .1*1.98847e30
    rho, A, Z, T = 6782.7e3, 4., 2., 100.
    ni, N = rho/(A*mu), mass/(A*mu)
    tp = hbar/kb*math.sqrt(ni*(Z*charge)**2/(eps0*A*mu))
    ef = hbar*hbar/(2*me)*(3*math.pi**2*Z*ni)**(2/3)
    uion = N*3*math.pi**4/5*kb*T**4/(.45*tp)**3
    ue = Z*N*math.pi**2/4*(kb*T)**2/ef
    relative = abs(r.U(T)/(uion+ue)-1)
    assert relative < 1e-12 and abs(r.T_p/tp-1) < 1e-12
    solid = r.U(r.T_m)
    latent = r.N_ion*r.L_ion
    gaps = []
    for fraction in [.1, .5, .9]:
        target = solid + fraction*latent
        returned = r.T_of_U(target)
        recovered = r.U(returned)
        gaps.append(dict(liquid_fraction_required_at_fixed_density=fraction,
                         target_energy_J=target, returned_temperature_K=returned,
                         returned_energy_J=recovered,
                         energy_residual_over_latent_heat=(recovered-target)/latent,
                         relative_energy_residual=(recovered-target)/target))
    # These are findings of a failing inversion, not checks declared passed.
    assert all(abs(v['energy_residual_over_latent_heat']) > .09 for v in gaps)
    output = dict(
        scope=__doc__, outcome='correction_partly_accepted_with_new_inversion_defect',
        input_sha256=inputs,
        accepted_under_stated_uniform_density_model=dict(
            plasma_temperature_K=tp, total_energy_100K_J=uion+ue,
            ion_energy_100K_J=uion, electron_energy_100K_J=ue,
            ion_fraction_100K=uion/(uion+ue), relative_energy_difference=relative),
        melting_model=dict(classical_melting_temperature_K=r.T_m,
                           plasma_to_melting_temperature_ratio=r.T_p/r.T_m,
                           solid_energy_at_melting_J=solid, latent_heat_total_J=latent,
                           inversion_gap_cases=gaps),
        review_decisions=[
            'Accept the corrected plasma temperature and low-temperature energy as conditional analytic benchmarks, not a cooled Ember structure.',
            'Require a phase-fraction state for energies between solid and liquid energies at melting; a root finder across an energy jump silently violates energy conservation.',
            'Classical Gamma=175 melting and constant liquid ionic heat capacity need quantum and liquid-correlation comparisons before assigning helium melting or pulse temperatures.',
            'Keep effective temperature, interior temperature and cooling age separate; retain the conditional conversion of mode energy to heat.',
            'No sealed SPH orbit needs rerunning for this correction.'],
        primary_reference='https://arxiv.org/abs/1910.06771',
        accepted_for_stellar_evolution=False)
    path = ROOT/'docs/results/fable_thermal_correction_review_v2.json'
    with path.open('x') as stream:
        json.dump(output, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({'outcome': output['outcome'], 'U100_J': uion+ue,
                      'Tp_over_Tm': r.T_p/r.T_m,
                      'gap_relative_energy_residuals': [v['relative_energy_residual'] for v in gaps]}))


if __name__ == '__main__':
    main()
