#!/usr/bin/env python3
"""Read and validate SYNSPEC's original binary absorption tables.

This is the rectangular density-grid format used by Ember's source recipe.
Fortran record markers are checked, including at EOF. Values are natural
logarithms; opacity storage order is wavelength, density, temperature.
"""
from array import array
import math
from pathlib import Path
import struct
import sys

# The original TLUSTY reader's temporary arrays support 100 temperature
# rows. Its compiled atmosphere capacity is checked separately at launch.
MAX_TEMPERATURE_ROWS = 100


def read_table(path):
    with Path(path).open("rb") as source:
        def record(size):
            marker = source.read(4)
            if len(marker) != 4 or struct.unpack("<i", marker)[0] != size:
                raise ValueError("invalid or truncated opacity record")
            value = source.read(size)
            if len(value) != size or source.read(4) != marker:
                raise ValueError("truncated opacity record or mismatched marker")
            return value

        abundance = [struct.unpack("<4sdd", record(20)) for _ in range(92)]
        ifmol, tmolim = struct.unpack("<id", record(12))
        flags = struct.unpack("<10i", record(40))
        nf, nt, nr = struct.unpack("<3i", record(12))
        if not (2 <= nf <= 100000 and 1 <= nt <= MAX_TEMPERATURE_ROWS and 2 <= nr <= 19):
            raise ValueError("unsupported opacity table dimensions")
        lt = struct.unpack(f"<{nt}d", record(8*nt))
        lr = struct.unpack(f"<{nr}d", record(8*nr))
        le = struct.unpack(f"<{nt*nr}d", record(8*nt*nr))
        frequency = array("d")
        opacity = array("f")
        for _ in range(nf):
            frequency.append(struct.unpack("<d", record(8))[0])
            for _ in range(nr):
                opacity.frombytes(record(4*nt))
        if sys.byteorder != "little":
            opacity.byteswap()
        if source.read(1):
            raise ValueError("trailing opacity data")
    for axis in [lt, lr]:
        if any(not math.isfinite(v) for v in axis) or any(a >= b for a,b in zip(axis,axis[1:])):
            raise ValueError("invalid opacity material axis")
    if any(not math.isfinite(v) or v <= 0 for v in frequency) or any(a <= b for a,b in zip(frequency,frequency[1:])):
        raise ValueError("invalid opacity frequency axis")
    if any(not math.isfinite(v) for v in opacity) or any(not math.isfinite(v) for v in le):
        raise ValueError("nonfinite source opacity or electron density")
    return {"abundance": abundance, "ifmol": ifmol, "tmolim": tmolim,
            "flags": flags, "shape": (nf,nr,nt), "log_temperature": lt,
            "log_density": lr, "log_electron_density": le,
            "frequency": frequency, "log_opacity": opacity}


def merge_isotherms(paths, output):
    """Assemble independent source isotherms without interpolating opacity."""
    tables = [read_table(p) for p in paths]
    if not 2 <= len(tables) <= MAX_TEMPERATURE_ROWS or any(t["shape"][2] != 1 for t in tables):
        raise ValueError("expected a complete set of source isotherms")
    first = tables[0]
    for table in tables[1:]:
        for key in ["abundance","ifmol","tmolim","flags","shape","log_density","frequency"]:
            if table[key] != first[key]:
                raise ValueError("incompatible source isotherms")
    temperatures = [t["log_temperature"][0] for t in tables]
    if any(a >= b for a,b in zip(temperatures,temperatures[1:])):
        raise ValueError("unordered source isotherms")
    for axis in [temperatures, first["log_density"]]:
        step=(axis[-1]-axis[0])/(len(axis)-1)
        if any(abs(v-axis[0]-i*step)>1e-10 for i,v in enumerate(axis)):
            raise ValueError("TLUSTY requires uniformly spaced logarithmic material axes")
    temporary = Path(str(output)+".part")
    with temporary.open("wb") as target:
        def record(fmt, *values):
            content = struct.pack("<"+fmt, *values)
            marker = struct.pack("<i",len(content))
            target.write(marker+content+marker)
        for row in first["abundance"]:record("4sdd",*row)
        record("id",first["ifmol"],first["tmolim"])
        record("10i",*first["flags"])
        nf,nr,_ = first["shape"]; nt = len(tables)
        record("3i",nf,nt,nr)
        record(f"{nt}d",*temperatures)
        record(f"{nr}d",*first["log_density"])
        record(f"{nt*nr}d",*(v for t in tables for v in t["log_electron_density"]))
        for k,f in enumerate(first["frequency"]):
            record("d",f)
            for j in range(nr):
                record(f"{nt}f",*(t["log_opacity"][k*nr+j] for t in tables))
    temporary.replace(output)


def validate_table(path, abundance, temperatures, densities):
    table = read_table(path)
    for (_, eos, absorption), expected in zip(table["abundance"], abundance):
        if not math.isclose(eos, expected, rel_tol=1e-12) or not math.isclose(absorption, expected, rel_tol=1e-12):
            raise ValueError("opacity source composition mismatch")
    for name, expected in [("log_temperature",temperatures),("log_density",densities)]:
        actual = table[name]
        if len(actual) != len(expected) or any(abs(a-math.log(b)) > 1e-10 for a,b in zip(actual,expected)):
            raise ValueError("opacity source support mismatch")
    if table["ifmol"] != 1 or table["tmolim"] != 10000 or table["flags"][-4:] != (1,1,1,1):
        raise ValueError("molecular opacity inputs missing")
    return table
