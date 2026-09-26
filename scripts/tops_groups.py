"""Read TOPS frequency groups and evaluate a controlled transport approximation.

The group Rosseland mean is retained exactly when the refractive index is one.
The group Planck absorption estimates the scattering fraction in a narrow
group; this estimate is an approximation to test against complete spectra.
"""
import json
import math
from pathlib import Path
import re

import numpy as np
from numpy.polynomial.legendre import leggauss

from audit_tops_spectral_means import digest
from import_tops_composition import read
from audit_tops_electron_dispersion import electron_moments, index_squared, KEV, KB, HBAR, ME, RW


def source_groups(root):
    root = Path(root)
    receipt = json.loads((root/'receipt.json').read_text())
    if digest(root/'request.json') != receipt['request_sha256']:
        raise ValueError('group request changed')
    request = json.loads((root/'request.json').read_text())
    if request['datype'] != 'groups' or request['plasnu'] != 'off' or request['egrid'] != 'range':
        raise ValueError('expected uncut frequency groups on an explicit range')
    _, _, cells, excluded = read(root/'source.txt', receipt, dimensions=tuple(receipt['dimensions']))
    if excluded:
        raise ValueError('substituted states are not valid group sources')
    text = (root/'source.txt').read_text()
    n = int(request['ngpengs'])
    if request['espace'] == 'log':
        boundaries = np.geomspace(float(request['egplow']), float(request['egphigh']), n)
    elif request['espace'] == 'linear':
        boundaries = np.linspace(float(request['egplow']), float(request['egphigh']), n)
    else:
        raise ValueError('unknown group spacing')
    data = {}
    for part in text.split('Energy      Ross mg     Planck mg    for T, density =')[1:]:
        lines = part.splitlines(); key = tuple(map(float, lines[0].split()))
        rows = []
        for line in lines[1:]:
            values = line.split()
            if len(values) != 3 or not all(re.fullmatch(r'[0-9.E+-]+', x) for x in values):
                break
            rows.append(list(map(float, values)))
        a = np.array(rows)
        if (a.shape != (n-1, 3) or not np.isfinite(a).all() or np.any(a <= 0)
                or np.max(np.abs(a[:, 0]/boundaries[:-1]-1)) > 5e-5 or key in data):
            raise ValueError('incomplete or inconsistent frequency groups')
        data[key] = a
    grey = {}
    for part in text.split('Density     Ross opa    Planck opa  No. Free    Av Sq Free  T=  ')[1:]:
        lines = part.splitlines(); t = float(lines[0])
        for line in lines[1:]:
            values = line.split()
            if len(values) != 5 or not all(re.fullmatch(r'[0-9.E+-]+', x) for x in values):
                break
            rho, ross, planck, free, _ = map(float, values)
            grey[t, rho] = {'rosseland': ross, 'planck': planck, 'free_electrons_per_ion': free}
    if set(data) != set(cells) or set(grey) != set(cells):
        raise ValueError('group and grey source coordinates disagree')
    composition = text.split('No. Fraction Mass Fraction  At. No.  Chem. Sym.  Mat ID.\n')[1].split('Temperature grid')[0]
    helium = [v.split() for v in composition.splitlines() if len(v.split()) == 5 and v.split()[3] == 'He']
    if len(helium) != 1:
        raise ValueError('expected helium in the stellar source mixture')
    number, mass = map(float, helium[0][:2])
    mean_atomic_mass = 4.002602*number/mass
    for key in grey:
        grey[key]['electron_density_cm3'] = key[1]/(mean_atomic_mass*1.66053906660e-24)*grey[key]['free_electrons_per_ion']
    return boundaries, data, grey


def transport(boundaries, groups, temperature_keV, ne, order=32):
    e = electron_moments(temperature_keV*KEV/KB, ne)
    up = HBAR*math.sqrt(4*math.pi*4.80320471257e-10**2*ne/ME)/(temperature_keV*KEV)
    up *= math.sqrt(e['plasma_frequency_squared_ratio'])
    if boundaries[0]/temperature_keV > up or boundaries[-1]/temperature_keV < max(100., up+60.):
        raise ValueError('insufficient group coverage for this plasma cutoff')
    lo, hi = np.maximum(boundaries[:-1]/temperature_keV, up), boundaries[1:]/temperature_keV
    valid = hi > lo; lo, hi, a = lo[valid], hi[valid], groups[valid]
    nodes, weights = leggauss(order)
    u = (lo[:, None]+hi[:, None])/2+(hi-lo)[:, None]/2*nodes
    n = np.sqrt(index_squared(u, up, e['vstar_squared']))
    wr = u**4*np.exp(-u)/(-np.expm1(-u))**2
    fraction = np.clip(1-a[:, 2]/a[:, 1], 0., 1.)
    common = (hi-lo)[:, None]/2*weights*wr/a[:, 1, None]
    component = RW/np.sum(common*n**3/(1-(1-n)*fraction[:, None]))
    common_index = RW/np.sum(common*n**3)
    return {'rosseland_component_estimate': float(component),
            'rosseland_common_index': float(common_index), 'cutoff_u': up,
            'groups_above_cutoff': int(valid.sum()),
            'groups_with_planck_above_rosseland': int(np.sum(a[:, 2] > a[:, 1])),
            'electron': e}
