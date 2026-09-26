#!/usr/bin/env python3
"""Compare uncorrected OP mixture means with a retained TOPS spectrum.

This diagnostic sets the charge and opacity of P, Cl, K and Ti to zero while
retaining their original reference mass. It adds none of OPserver's separate
mixture scattering or screening corrections, and no refractive correction.
It is not an accepted stellar opacity prescription.
"""
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np

from op_native_spectra_v2 import read_mesh, read_spectra, native_rosseland
from audit_tops_spectral_means import source

MU = 1.66053906660e-24
A0 = 5.29177210903e-9
KB = 1.380649e-16
KEV = 1.602176634e-9
DATA = Path('/tmp/ember-op-selected-spectra-v1/OP4STARS_1.3/mono')
WORK = Path('/tmp/ember-op-native-mixture-run-v1')
OUTPUT = Path('docs/results/op_native_mixture_comparison_v1.json')


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def interpolation(logx, logy, target):
    if len(logx) < 2 or not logx[0] <= target <= logx[-1]:
        raise ValueError('interpolation would leave the supported range')
    values = []
    for i in range(len(logx)):
        weight = 1.
        for j in range(len(logx)):
            if i != j:weight *= (target - logx[j]) / (logx[i] - logx[j])
        values.append(weight * logy[i])
    return float(sum(values))


def select_connected(indices, valid, target, number):
    if not np.all(np.diff(indices) > 0):raise ValueError('unordered source coordinates')
    # Find a contiguous sequence of valid states enclosing the requested point.
    segments = np.split(np.arange(len(indices)), np.flatnonzero(np.diff(valid.astype(int))) + 1)
    choices = [s for s in segments if len(s) >= number and valid[s[0]] and indices[s[0]] <= target <= indices[s[-1]]]
    if len(choices) != 1:raise ValueError('not enough contiguous positive source spectra')
    s = choices[0]; upper = int(np.searchsorted(indices[s], target, side='right'))
    first = min(max(upper - 1 - (number - 2)//2, 0), len(s) - number)
    return s[first:first+number]


def native_plane(it, weights, rho, dv, assets, inputs):
    sigma = None; charge = None; valid = None; jns = None; source_masks = []
    for z, w in weights.items():
        p = DATA / f'm{z:02d}.{it:03d}'
        r = read_spectra(p); inputs[str(p)] = digest(p)
        these = np.array([v['electron_index'] for v in r['states']])
        if jns is None:
            jns = these; sigma = np.zeros((len(jns), 10000)); charge = np.zeros(len(jns)); valid = np.ones(len(jns), dtype=bool)
        elif not np.array_equal(jns, these):raise ValueError('element state grids differ')
        for j, row in enumerate(r['states']):
            sigma[j] += w * row['cross_section_atomic']
            charge[j] += w * row['electrons_per_atom']
            if not row['positive_cross_section']:
                valid[j] = False; source_masks.append(dict(atomic_number=z, electron_index=int(jns[j])))
    if np.any(charge <= 0):raise ValueError('nonpositive mixture charge')
    density = MU * 10**(.25*jns) / charge
    if np.any(np.diff(density) <= 0):raise ValueError('nonmonotonic mixture density')
    valid &= np.all(sigma > 0, axis=1)
    means = np.full(len(jns), np.nan)
    means[valid] = A0*A0/MU * native_rosseland(sigma[valid], dv)
    rows = []
    for number in (2, 4):
        chosen = select_connected(np.log(density), valid, math.log(rho), number)
        kap = math.exp(interpolation(np.log(density[chosen]), np.log(means[chosen]), math.log(rho)))
        ne = math.exp(interpolation(np.log(density[chosen]), .25*jns[chosen]*math.log(10), math.log(rho)))
        rows.append(dict(density_points=number, electron_indices=jns[chosen].tolist(),
                         native_densities=density[chosen].tolist(), native_means_cm2_g=means[chosen].tolist(),
                         mean_at_requested_density_cm2_g=kap, electron_density_cm3=ne))
        for j in chosen:
            assets[f't{it}_n{jns[j]}'] = sigma[j]
    return dict(temperature_index=it, temperature_K=10**(.025*it), source_masks=source_masks, density_interpolations=rows)


def physical_case(label, its, weights, temperature, rho, dv, inputs):
    assets = {}; planes = [native_plane(it, weights, rho, dv, assets, inputs) for it in its]
    interpolations = []
    for nt in (2, 4):
        ts = [1,2] if nt == 2 else [0,1,2,3]
        for nr in (2,4):
            rows = [planes[i]['density_interpolations'][0 if nr == 2 else 1] for i in ts]
            logt = np.log([planes[i]['temperature_K'] for i in ts])
            k = math.exp(interpolation(logt, np.log([v['mean_at_requested_density_cm2_g'] for v in rows]), math.log(temperature)))
            ne = math.exp(interpolation(logt, np.log([v['electron_density_cm3'] for v in rows]), math.log(temperature)))
            interpolations.append(dict(temperature_points=nt, density_points=nr, rosseland_cm2_g=k, electron_density_cm3=ne))
    artifact = WORK / (label + '_native_mixture_spectra.npz')
    if artifact.exists():raise FileExistsError(artifact)
    np.savez_compressed(artifact, **assets)
    return dict(label=label, temperature_K=temperature, density_g_cm3=rho,
                included_atom_weights=weights, planes=planes, interpolations=interpolations,
                spectral_artifact=dict(path=str(artifact), bytes=artifact.stat().st_size, sha256=digest(artifact)))


def main():
    if OUTPUT.exists():raise FileExistsError(OUTPUT)
    validation = Path('docs/results/op_native_spectral_reader_v2.json')
    checked = json.loads(validation.read_text())
    if checked['failures']:raise ValueError('native reader validation failed')
    for p,h in checked['input_sha256'].items():
        if digest(p) != h:raise ValueError('validated OP input changed')
    inputs = {str(Path(p).resolve()): digest(p) for p in [__file__, 'scripts/op_native_spectra_v2.py',
              'scripts/audit_tops_spectral_means.py', validation]}
    dv, u = read_mesh(DATA/'m01.mesh'); inputs[str(DATA/'m01.mesh')] = digest(DATA/'m01.mesh')
    se = -np.expm1(-u)
    # An independently specified constant total opacity checks the stimulated
    # emission conversion and finite native frequency normalization.
    constant_opacity = 7.
    sigma_per_ref = constant_opacity * MU/(A0*A0) / se
    normalizing_weight = float(dv * np.sum(se))
    recovered = float(A0*A0/MU*native_rosseland(sigma_per_ref,dv))
    assert abs(recovered/(constant_opacity/normalizing_weight)-1) < 2e-15
    reference_report = Path('docs/results/tops_cool_full_spectra_v4.json')
    refs = json.loads(reference_report.read_text()); inputs[str(reference_report.resolve())] = digest(reference_report)
    ref = next(r for r in refs['records'] if r['X']==.1 and r['temperature_keV']==.008 and r['density_atomic_g_cm3']==.053367)
    root = Path(ref['full_spectrum_source']); spectra, _ = source(root)
    receipt = json.loads((root/'receipt.json').read_text())
    if digest(root/'source.txt') != receipt['sha256']:raise ValueError('TOPS source hash changed')
    for p in (root/'source.txt',root/'receipt.json',root/'request.json'):
        inputs[str(p)] = digest(p)
    text = (root/'source.txt').read_text().split('No. Fraction Mass Fraction  At. No.  Chem. Sym.  Mat ID.\n')[1].split('Temperature grid')[0]
    elements = [r.split() for r in text.splitlines() if len(r.split())==5]
    he = next(r for r in elements if r[3]=='He')
    mean_mass = 4.002602*float(he[0])/float(he[1])
    allowed = {int(p.name[1:3]) for p in DATA.glob('m*.mesh')}
    weights = {int(r[2]):float(r[0])/mean_mass for r in elements if int(r[2]) in allowed}
    tops_case = physical_case('tops_x010_t008000', [196,198,200,202], weights, .008*KEV/KB, .053367, dv, inputs)
    total = spectra[(.008,.053367)]
    tu = total[:,0]/.008
    if not tu[0] <= u[0] < u[-1] <= tu[-1]:raise ValueError('TOPS spectrum does not cover OP frequency interval')
    tops_k = np.exp(np.interp(np.log(u), np.log(tu), np.log(total[:,1])))
    tops_on_v = float(1/(dv*np.sum(se/tops_k)))
    tops_case.update(tops_full_uncut_cm2_g=ref['reference_uncut_rosseland_atomic_cm2_g'], tops_on_OP_mesh_cm2_g=tops_on_v,
                     tops_mesh_change=tops_on_v/ref['reference_uncut_rosseland_atomic_cm2_g']-1,
                     omitted_mass_fraction=sum(float(r[1]) for r in elements if int(r[2]) not in allowed),
                     comparison=[dict(**r, relative_to_TOPS_on_OP_mesh=r['rosseland_cm2_g']/tops_on_v-1) for r in tops_case['interpolations']])
    domain = Path('docs/results/op_native_stellar_domain_v2.json')
    star = json.loads(domain.read_text())['stellar_records'][493]; inputs[str(domain.resolve())] = digest(domain)
    metals_file=Path('include/ember/gs98_mixture.hpp');inputs[str(metals_file.resolve())]=digest(metals_file)
    metals=[tuple(map(float,r.split(','))) for r in re.findall(r'^\s*\{([^{}]+)\}, //',metals_file.read_text(),re.M)]
    weights={1:star['X'],2:star['Y3']/3+star['Y4']/4}
    weights.update({int(r[0]):star['Z']*r[3]/r[1] for r in metals if int(r[0]) in allowed})
    star_case=physical_case('stellar_zone493',[186,188,190,192],weights,star['temperature_K'],star['density_baryonic_g_cm3'],dv,inputs)
    star_case.update(density_convention='baryonic mass from tracked isotope counts', omitted_mass_fraction=star['missing_metal_mass_fraction'])
    for p,h in inputs.items():
        if digest(p)!=h:raise ValueError('input changed during comparison')
    result=dict(scope=__doc__,outcome='completed_uncorrected_mixture_comparison',accepted_for_stellar_opacity=False,
                constant_opacity_control=dict(target=constant_opacity,native_frequency_weight=normalizing_weight,recovered_finite_mesh=recovered),
                input_sha256=inputs,cases=[tops_case,star_case],
                limitations=['The OP values omit four trace-element spectra and all separate OPserver mixture scattering/screening corrections; differences cannot be assigned solely to atomic opacity physics.',
                             'Two- versus four-point interpolation measures sensitivity on the existing coarse grid, not independent source accuracy.',
                             'Native spectra end at u=20; comparison uses the same frequency mesh for TOPS and reports its difference from the retained broad-spectrum mean.',
                             'The negative helium source state is excluded as a whole; only contiguous positive source states can enter interpolation.',
                             'No bound/free absorption-scattering decomposition, refractive transport correction, new atmosphere or stellar continuation is accepted.'])
    OUTPUT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],tops_comparison=tops_case['comparison'],tops_mesh_change=tops_case['tops_mesh_change'],star=star_case['interpolations'])),flush=True)


if __name__=='__main__':main()
