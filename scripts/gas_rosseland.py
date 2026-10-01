#!/usr/bin/env python3
"""Integrate H/He absorption, electron scattering and Rayleigh scattering.

Inputs are SYNSPEC absorption isotherms and neutral populations exported from
the same chemical calculation. This does not correct the source's dense-gas
chemistry, add grains, or model correlated scattering.
"""
import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from nongrey_opacity import read_table


def population_export(source):
    """Add diagnostic output to SYNSPEC54 without changing its equations."""
    marker = "c     Ember neutral populations for Rosseland scattering.\n"
    if marker in source:
        return source
    start = source.lower().index("      subroutine ougrid(abso)")
    end = source.lower().index("\n      end\n", start)
    anchor = re.search(r"(?im)^      if \(nfreq\.le\.3\) return[^\S\n]*\n", source[start:end+1])
    if anchor is None:
        raise ValueError("unrecognized SYNSPEC ougrid routine")
    pos = start + anchor.end()
    return source[:pos] + marker + (
        "      if(ipfreq.eq.0) write(184,*) temp(1),dens(1),elec(1),\n"
        "     *     popul(nfirst(ielh),1),anh2(1),rrr(1,1,2),\n"
        "     *     dens(1)/(wmm(1)*ytot(1))\n"
    ) + source[pos:]


def rayleigh_cross_sections(frequency):
    """H ground state, H2 and He I; SYNSPEC54 OPADD/TLUSTY208 RAYLEIGH."""
    nu = np.asarray(frequency)
    v = (np.minimum(nu, 2.463e15)/2.997925e18)**2
    hydrogen = (5.799e-13 + (1.422e-6 + 2.784*v)*v)*v*v
    v = (2.997925e18/np.minimum(nu, 5.150e15))**2
    helium = 5.484e-14/v**2 * (1 + (2.44e5 + 5.94e10/(v-2.90e5))/v)**2
    v = (2.997925e18/np.minimum(nu, 2.922e15))**2
    molecular = (8.14e-13 + 1.28e-6/v + 1.61/v**2)/v**2
    return hydrogen, molecular, helium


def rosseland(frequency, temperature, opacity):
    """Harmonic mean weighted by dB_nu/dT; frequency must increase."""
    nu = np.asarray(frequency)
    opacity = np.asarray(opacity)
    if (not np.isfinite(temperature) or temperature <= 0
            or nu.ndim != 1 or len(nu) < 2 or not np.isfinite(nu).all()
            or not np.all(np.diff(nu) > 0) or nu[0] <= 0
            or opacity.shape[0] != len(nu) or not np.all(opacity > 0)
            or not np.isfinite(opacity).all()):
        raise ValueError("positive opacity and increasing positive frequencies required")
    x = 6.62607015e-27*nu/(1.380649e-16*temperature)
    w = nu**4*np.exp(-x)/np.square(-np.expm1(-x))
    weights = w.reshape((len(nu),) + (1,)*(opacity.ndim-1))
    return np.trapezoid(w, nu)/np.trapezoid(weights/opacity, nu, axis=0)


def isotherm(table, populations):
    nf, nr, nt = table['shape']
    p = np.asarray(populations)
    if nt != 1 or p.shape != (nr, 7) or not np.isfinite(p).all() or np.any(p < 0):
        raise ValueError("one isotherm and seven population columns per density required")
    temperature = np.exp(table['log_temperature'][0])
    rho = np.exp(table['log_density'])
    ne = np.exp(table['log_electron_density'])
    # The source's density/electron iteration has finite accuracy. Record the
    # discrepancy; reject unrelated populations instead of interpolating them.
    errors = [float(np.max(abs(p[:,i]/ref-1)))
              for i, ref in enumerate([temperature, rho, ne])]
    if max(errors) > 1e-5:
        raise ValueError(f"population coordinates/electron density disagree: {errors}")
    nu = np.asarray(table['frequency'])[::-1]
    absorption = np.exp(np.asarray(table['log_opacity']).reshape(nf, nr)[::-1])
    h, h2, he = rayleigh_cross_sections(nu)
    scattering = (h[:,None]*p[None,:,3] + h2[:,None]*p[None,:,4]
                  + he[:,None]*p[None,:,5] + 6.6524587321e-25*ne[None,:])/rho
    # Column 6 is total particle density. It is not neutral atomic hydrogen.
    return rosseland(nu, temperature, absorption + scattering), errors


def write_opacity_table(rows, output):
    """Write a rectangular table or explicitly supported density prefixes."""
    xs, ts, rs = [sorted({key[i] for key in rows}) for i in range(3)]
    if min(len(ts), len(rs)) < 2 or any(not np.isfinite(v) or v <= 0 for v in rows.values()):
        raise ValueError("positive material values and at least two temperatures/densities required")
    counts, densities = {}, {}
    for x, t, r in rows:
        densities.setdefault((x,t), []).append(r)
    for x in xs:
        for t in ts:
            actual = sorted(densities.get((x,t), []))
            if len(actual) < 2 or actual != rs[:len(actual)]:
                raise ValueError("each isotherm must contain a contiguous density prefix")
            counts[x,t] = len(actual)
    ragged = any(n != len(rs) for n in counts.values())
    if ragged and (len(ts) < 4 or min(counts.values()) < 4):
        raise ValueError("masked tables require four temperatures and four densities per row")
    text = [f'{len(xs)} {len(ts)} {len(rs)} HHe_absorption_electron_Rayleigh',
            ' '.join(format(r, '.17g') for r in rs),
            ' '.join(format(np.log10(t), '.17g') for t in ts)]
    for x in xs:
        text.append(format(x, '.17g') + ' 0')
        for t in ts:
            values = ' '.join(format(np.log10(rows[x,t,r]), '.17g') for r in rs[:counts[x,t]])
            text.append((str(counts[x,t])+' ' if ragged else '') + values)
    if ragged: text.insert(0, 'EMBER_OPACITY_TABLE 2')
    output.write_text('\n'.join(text) + '\n')


def integrate(manifest, output):
    spec = json.loads(manifest.read_text())
    if spec.get('absorption_only') is not True:
        raise ValueError("inputs must declare scattering-free absorption spectra")
    rows, sources = {}, []
    for column in spec['columns']:
        paths = {}
        for kind in ['absorption', 'populations']:
            item = column[kind]
            path = manifest.parent/item['path']
            if hashlib.sha256(path.read_bytes()).hexdigest() != item['sha256']:
                raise ValueError(f"changed {kind} input: {path}")
            paths[kind] = path
        table = read_table(paths['absorption'])
        abundance = table['abundance']
        x, t = column['X'], column['T_K']
        actual_x = abundance[0][1]/(abundance[0][1] + 4*abundance[1][1])
        if (abs(x-actual_x) > 1e-8 or any(a[1] > 1e-30 for a in abundance[2:])
                or abs(t/np.exp(table['log_temperature'][0])-1) > 1e-8):
            raise ValueError("H/He abundance or temperature label mismatch")
        mean, errors = isotherm(table, np.loadtxt(paths['populations']))
        indices = column.get('density_indices', list(range(len(mean))))
        if (not isinstance(indices, list) or not indices
                or any(type(i) is not int or i < 0 or i >= len(mean) for i in indices)
                or len(set(indices)) != len(indices)):
            raise ValueError("density_indices must select distinct source rows")
        for i in indices:
            lr, value = table['log_density'][i], mean[i]
            key = x, t, round(lr/np.log(10), 10)
            if key in rows:
                raise ValueError("duplicate material node")
            rows[key] = float(value)
        sources.append(dict(column, population_relative_errors=errors))
    write_opacity_table(rows, output)
    return dict(nodes=len(rows), columns=sources,
                sha256=hashlib.sha256(output.read_bytes()).hexdigest())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    patch = commands.add_parser('export-populations')
    patch.add_argument('source', type=Path); patch.add_argument('output', type=Path)
    table = commands.add_parser('integrate')
    table.add_argument('manifest', type=Path); table.add_argument('output', type=Path)
    table.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'export-populations':
        args.output.write_text(population_export(args.source.read_text()))
    else:
        args.report.write_text(json.dumps(integrate(args.manifest, args.output), indent=2)+'\n')
