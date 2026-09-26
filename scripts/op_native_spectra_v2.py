#!/usr/bin/env python3
"""Read archived OP sequential Fortran records without executing OP software.

Native cross sections are in a0^2, on the distributed v mesh, without the
stimulated-emission factor. This module adds no mixture corrections and
performs no interpolation in thermodynamic state. Nonpositive source spectra
are returned with an explicit flag and cannot be passed to the mean integrator.
"""
from pathlib import Path
import struct

import numpy as np


class Records:
    def __init__(self, path, expected_first):
        self.body = Path(path).read_bytes()
        candidates = [endian for endian in ('<', '>') if len(self.body) >= 4
                      and struct.unpack(endian + 'i', self.body[:4])[0] == expected_first]
        if len(candidates) != 1:
            raise ValueError('unrecognized Fortran record markers')
        self.endian = candidates[0]; self.offset = 0

    def record(self):
        end = self.endian; p = self.offset; data = self.body
        if p + 8 > len(data):
            raise ValueError('missing Fortran record')
        size = struct.unpack_from(end + 'i', data, p)[0]
        if size < 0 or p + size + 8 > len(data):
            raise ValueError('truncated Fortran record')
        tail = struct.unpack_from(end + 'i', data, p + 4 + size)[0]
        if tail != size:
            raise ValueError('Fortran record marker mismatch')
        self.offset = p + size + 8
        return data[p + 4:p + 4 + size]

    def finished(self):
        if self.offset != len(self.body):
            raise ValueError('unexpected trailing records')


def read_mesh(path):
    r = Records(path, 40008)
    record = r.record()
    dv, ntot = struct.unpack_from(r.endian + 'fi', record)
    u = np.frombuffer(record, dtype=r.endian + 'f4', offset=8).astype(float)
    r.finished()
    if ntot != 10000 or len(u) != ntot or not dv > 0 or not np.isfinite(u).all() or np.any(np.diff(u) <= 0) or u[0] <= 0:
        raise ValueError('invalid frequency mesh')
    return dv, u


def read_spectra(path):
    r = Records(path, 44)
    header = struct.unpack(r.endian + 'iifffiifiii', r.record())
    z, it, mass, umin, umax, ncoarse, ntot, tolerance, j0, j1, dj = header
    if Path(path).name != f'm{z:02d}.{it:03d}' or ntot != 10000 or dj <= 0:
        raise ValueError('wrong source header')
    states = []
    for j in range(j0, j1 + 1, dj):
        row = r.record()
        jn, electrons, planck, rosseland, n0, n1 = struct.unpack_from(r.endian + 'ifffii', row)
        ions = np.frombuffer(row, dtype=r.endian + 'f4', offset=24).astype(float)
        if jn != j or n0 < -1 or n1 > 30 or len(ions) != n1 - n0 + 1 or not np.isfinite(ions).all() or np.any(ions < 0):
            raise ValueError('wrong density/ionization record')
        npacked, = struct.unpack(r.endian + 'i', r.record())
        packed = r.record()
        if npacked == 0:
            spectrum = np.frombuffer(packed, dtype=r.endian + 'f4').astype(float)
        elif 2 <= npacked <= ntot:
            pairs = np.frombuffer(packed, dtype=[('i', r.endian + 'i4'), ('f', r.endian + 'f4')])
            indices = pairs['i']; values = pairs['f'].astype(float)
            if len(pairs) != npacked or indices[0] != 1 or indices[-1] != ntot or np.any(np.diff(indices) <= 0):
                raise ValueError('invalid packed spectrum')
            spectrum = np.interp(np.arange(1, ntot + 1), indices, values)
        else:
            raise ValueError('invalid packed size')
        if len(spectrum) != ntot or not np.isfinite(spectrum).all():
            raise ValueError('nonpositive or invalid element cross section')
        states.append(dict(electron_index=jn, electrons_per_atom=electrons,
                           planck_atomic=planck, rosseland_atomic=rosseland,
                           ion_indices=np.arange(n0, n1 + 1), ion_fractions=ions,
                           packed_points=npacked, positive_cross_section=bool(np.all(spectrum > 0)), cross_section_atomic=spectrum))
    r.finished()
    return dict(atomic_number=z, temperature_index=it, atomic_mass=mass,
                u_min=umin, u_max=umax, packing_tolerance=tolerance, states=states)


def native_rosseland(cross_section, dv):
    """Finite distributed-mesh mean; not an infinite-frequency completion."""
    cross_section = np.asarray(cross_section)
    if not np.isfinite(cross_section).all() or np.any(cross_section <= 0):
        raise ValueError('nonpositive source spectrum is masked, not integrated')
    return 1 / (dv * np.sum(1 / cross_section, axis=-1))
