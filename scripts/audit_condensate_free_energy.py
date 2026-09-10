#!/usr/bin/env python3
"""Audit an equilibrium-condensation Gibbs correction at fixed bulk composition.

Delta(g/T) = k_B sum_j (N_j/M) ln(n_j,eq/n_j,gas), using neutral atomic
chemical potentials. Atomic reference partition functions cancel between
the two equilibria. The same potential gives volume, entropy and enthalpy
corrections; gas-species rearrangement is included along with condensates.

This is an offline ideal-chemistry correction, not an installed atmosphere
EOS. Condensates retain mass, occupy negligible volume and remain in local
equilibrium. Rainout, grain opacity, transport and nonideal chemistry are
separate physics. A native gas EOS would need this correction applied to its
bulk-composition potential and all its derivatives, not merely its cp.
"""
import argparse
import gzip
import json
import math
from pathlib import Path
import subprocess

from estimate_condensate_enthalpy import coefficients, formation_enthalpy
from generate_nongrey_grid import composition, SOURCE_HMASS
from prepare_nongrey_sources import SOURCES, digest

KB = 1.380649e-16
OFFSETS = [(0, 0), (1, 0), (-1, 0), (2, 0), (-2, 0),
           (0, 1), (0, -1), (0, 2), (0, -2),
           (1, 1), (1, -1), (-1, 1), (-1, -1)]


def gas_coefficients(path):
    lines = path.read_text().splitlines()
    result = {}
    for i, line in enumerate(lines):
        if not line.strip() or line.startswith('#') or ':' not in line:
            continue
        fit = list(map(float, lines[i+1].split()))
        if len(fit) != 5 or not all(math.isfinite(v) for v in fit):
            raise ValueError('invalid gas thermochemical fit')
        result[line.split()[0]] = fit
    return result


def state(chemistry, row, masses, gas_fits, solid_fits):
    if row['flag'] or not all(row['element_conserved']):
        raise ValueError('chemical equilibrium/conservation failure')
    symbols, abundance = chemistry['elements'], chemistry['abundances']
    total_mass = sum(a*masses[s] for s, a in zip(symbols, abundance, strict=True) if s != 'e-')
    counts = {s: a/total_mass for s, a in zip(symbols, abundance, strict=True) if s != 'e-'}
    gas = {s['symbol']: n for s, n in zip(chemistry['gas_species'], row['gas'], strict=True)}
    atoms = {s: gas[s] for s in counts}
    if any(not math.isfinite(v) or v <= 0 for v in atoms.values()):
        raise ValueError('neutral-atom chemical potential underflow')
    rho = 0.; condensed_mass = 0.; formation = 0.; solid_mass = {}
    for category, names, densities in [('gas', chemistry['gas_species'], row['gas']),
                                      ('solid', chemistry['condensates'], row['condensed'])]:
        for species, number in zip(names, densities, strict=True):
            if not math.isfinite(number) or number < 0:
                raise ValueError('invalid chemical species density')
            particle_mass = sum(n*masses[s] for s, n in zip(symbols, species['stoichiometry'], strict=True) if s != 'e-')
            rho += number*particle_mass
            t = row['T_K']; symbol = species['symbol']
            if category == 'solid':
                condensed_mass += number*particle_mass
                solid_mass[symbol] = number*particle_mass
                h = formation_enthalpy(t, solid_fits[symbol]) if number else 0.
            elif symbol in gas_fits:
                a = gas_fits[symbol]
                h = KB*(-a[0]+a[1]*t+a[3]*t*t+2*a[4]*t**3)
            elif symbol in symbols:
                h = 0.  # Reference atoms/electron; their bulk terms cancel.
            else:
                raise ValueError('gas species lacks its thermochemical fit')
            formation += number*h
    pressure = sum(row['gas'])*KB*row['T_K']
    if rho <= 0 or abs(pressure/(row['P_bar']*1e6)-1) > 2e-7:
        raise ValueError('chemical pressure closure failed')
    return {'atom_number_density': atoms, 'bulk_counts_per_gram': counts,
            'rho_total': rho, 'formation_enthalpy': formation/rho,
            'condensed_mass_fraction': condensed_mass/rho,
            'active_condensates': sorted(s for s, v in solid_mass.items() if v/rho > 1e-10)}


def correction(gas, equilibrium):
    if gas['bulk_counts_per_gram'] != equilibrium['bulk_counts_per_gram']:
        raise ValueError('gas and condensate calculations have different bulk mixtures')
    phi = KB*math.fsum(n*(math.log(equilibrium['atom_number_density'][s])-math.log(gas['atom_number_density'][s]))
                      for s, n in gas['bulk_counts_per_gram'].items())
    return {'delta_g_over_T': phi,
            'delta_specific_volume': 1/equilibrium['rho_total']-1/gas['rho_total'],
            'delta_formation_enthalpy': equilibrium['formation_enthalpy']-gas['formation_enthalpy'],
            'gas_rho': gas['rho_total'], 'equilibrium_rho': equilibrium['rho_total'],
            'active_condensates': equilibrium['active_condensates'],
            'condensed_mass_fraction': equilibrium['condensed_mass_fraction']}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('model', type=Path, help='source model supplies only pinned chemistry and bulk abundances')
    p.add_argument('work', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--states', type=Path, required=True, help='JSON list of temperature_K, pressure_bar records')
    p.add_argument('--steps', type=float, nargs='+', default=[.008, .004, .002, .001],
                   help='decreasing positive logarithmic finite-difference steps')
    a = p.parse_args()
    if not a.steps or not all(math.isfinite(h) and 0 < h < .1 for h in a.steps) or any(
            a.steps[i] >= a.steps[i-1] for i in range(1, len(a.steps))):
        raise ValueError('invalid derivative refinement sequence')
    if a.work.exists():
        raise FileExistsError('use a new thermodynamic audit directory')
    provenance = json.loads((a.model/'provenance.json').read_text())
    fc = provenance['prepared']['depletion']['fastchem']; source = Path(fc['source'])
    if digest(fc['executable']) != fc['executable_sha256']:
        raise ValueError('chemistry executable changed')
    dependencies = {str((a.model/'provenance.json').resolve()): digest(a.model/'provenance.json'),
                    str((a.model/'condensate-abundances.dat').resolve()): digest(a.model/'condensate-abundances.dat'),
                    str(a.states.resolve()): digest(a.states), str(Path(fc['executable']).resolve()): fc['executable_sha256']}
    for name, expected in fc['chemistry_data_sha256'].items():
        path = source/'input/logK'/name
        if digest(path) != expected:
            raise ValueError('chemistry data changed')
        dependencies[str(path.resolve())] = expected
    _, weights = composition(provenance['XH'], provenance['X3'], provenance['specification']['metals'])
    symbols = json.loads((SOURCES/'synple-elements.json').read_text())['symbol']
    masses = {s.capitalize(): m*SOURCE_HMASS for s, m in zip(symbols, weights, strict=True)}
    gas_fits = gas_coefficients(source/'input/logK/logK.dat')
    solid_fits = coefficients(source/'input/logK/logK_condensates.dat')
    centers = json.loads(a.states.read_text()); steps = a.steps
    points = [(r['temperature_K']*math.exp(i*h), r['pressure_bar']*math.exp(j*h))
              for r in centers for h in steps for i, j in OFFSETS]
    inputs = ''.join(f'{t:.17g} {p:.17g}\n' for t, p in points)
    a.work.mkdir(parents=True); (a.work/'input.dat').write_text(inputs)
    (a.work/'audit_source.py').write_bytes(Path(__file__).read_bytes())
    sources = {}
    for mode in ['gas', 'equilibrium']:
        result = subprocess.run([fc['executable'], str(source), str((a.model/'condensate-abundances.dat').resolve()), mode],
                                input=inputs, capture_output=True, text=True, check=True)
        (a.work/(mode+'.json.gz')).write_bytes(gzip.compress(result.stdout.encode(), mtime=0))
        (a.work/(mode+'.log')).write_text(result.stderr)
        sources[mode] = json.loads(result.stdout)
        if sources[mode]['status'] or len(sources[mode]['rows']) != len(points):
            raise ValueError('incomplete chemical source result')
    bulk_count = sum(state(sources['gas'], sources['gas']['rows'][0], masses, gas_fits, solid_fits)['bulk_counts_per_gram'].values())
    corrections = [correction(state(sources['gas'], g, masses, gas_fits, solid_fits),
                              state(sources['equilibrium'], e, masses, gas_fits, solid_fits))
                   for g, e in zip(sources['gas']['rows'], sources['equilibrium']['rows'], strict=True)]
    results = []; cursor = 0
    for center in centers:
        scans = []; t, pressure = center['temperature_K'], center['pressure_bar']*1e6
        for h in steps:
            stencil = dict(zip(OFFSETS, corrections[cursor:cursor+len(OFFSETS)], strict=True)); cursor += len(OFFSETS)
            f = {c: r['delta_g_over_T'] for c, r in stencil.items()}; middle = stencil[0, 0]
            dx = (f[-2, 0]-8*f[-1, 0]+8*f[1, 0]-f[2, 0])/(12*h)
            dy = (f[0, -2]-8*f[0, -1]+8*f[0, 1]-f[0, 2])/(12*h)
            dxx = (-f[2, 0]+16*f[1, 0]-30*f[0, 0]+16*f[-1, 0]-f[-2, 0])/(12*h*h)
            dyy = (-f[0, 2]+16*f[0, 1]-30*f[0, 0]+16*f[0, -1]-f[0, -2])/(12*h*h)
            dxy = (f[1, 1]-f[1, -1]-f[-1, 1]+f[-1, -1])/(4*h*h)
            heat = {c: r['delta_formation_enthalpy'] for c, r in stencil.items()}
            cp_independent = (heat[-2, 0]-8*heat[-1, 0]+8*heat[1, 0]-heat[2, 0])/(12*h*t)
            volume = t*dy/pressure; enthalpy = -t*dx; cp = -dx-dxx
            scales = {'gas_volume': 1/middle['gas_rho'], 'kT_per_bulk_atom': KB*t*bulk_count}
            scans.append({'logarithmic_step': h, 'central': middle,
                          'phase_inventory_changes_in_thermal_stencil': any(
                              stencil[c]['active_condensates'] != middle['active_condensates'] for c in [(1, 0), (-1, 0), (2, 0), (-2, 0)]),
                          'delta_entropy': -middle['delta_g_over_T']-dx,
                          'delta_enthalpy': enthalpy, 'delta_cp': cp,
                          'delta_cp_from_species_enthalpies': cp_independent,
                          'delta_specific_volume': volume,
                          'delta_dvolume_dT_at_P': (dy+dxy)/pressure,
                          'delta_dvolume_dP_at_T': t*(dyy-dy)/pressure**2,
                          'pressure_identity_error_in_gas_volume': (volume-middle['delta_specific_volume'])/scales['gas_volume'],
                          'enthalpy_identity_error_in_kT_per_bulk_atom': (enthalpy-middle['delta_formation_enthalpy'])/scales['kT_per_bulk_atom']})
        results.append({**center, 'scans': scans})
    if any(digest(path) != expected for path, expected in dependencies.items()):
        raise ValueError('chemical inputs changed during audit')
    report = {'scope': __doc__, 'script_sha256': digest(__file__),
              'sources': ['https://doi.org/10.1093/mnras/stad3515', 'https://github.com/NewStrangeWorlds/FastChem/tree/'+fc['commit']],
              'input_files_sha256': dependencies, 'states': results,
              'source_outputs_sha256': {name: digest(a.work/name) for name in ['audit_source.py', 'input.dat', 'gas.json.gz', 'equilibrium.json.gz']},
              'runtime_status': 'not installed; numerical and physical joins still require validation'}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps([{'T': r['temperature_K'], 'P_bar': r['pressure_bar'],
                       'delta_cp_by_step': [s['delta_cp'] for s in r['scans']],
                       'cp_from_species_by_step': [s['delta_cp_from_species_enthalpies'] for s in r['scans']],
                       'max_pressure_identity_error': max(abs(s['pressure_identity_error_in_gas_volume']) for s in r['scans']),
                       'max_enthalpy_identity_error': max(abs(s['enthalpy_identity_error_in_kT_per_bulk_atom']) for s in r['scans'])} for r in results], indent=2))


if __name__ == '__main__':
    main()
