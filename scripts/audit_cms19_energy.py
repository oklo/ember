#!/usr/bin/env python3
"""Reproduce a CMS19 caloric-energy join defect from the original archive.

No imported values are repaired. Output is JSON so future EOS work can
compare a newer source against this exact diagnostic.
"""
import hashlib
import json
from pathlib import Path
import sys
import tarfile
from import_cms19 import SHA256

source = Path(sys.argv[1])
if hashlib.sha256(source.read_bytes()).hexdigest() != SHA256:
    raise ValueError('source archive differs from documented CMS19 version')
with tarfile.open(source) as archive:
    rows = [list(map(float, line.split())) for line in archive.extractfile(
        'DirTABLES-EOS2019/TABLE_HE_Trho_v1').read().decode().splitlines()
        if line.strip() and not line.startswith('#')]
    hrows = [list(map(float, line.split())) for line in archive.extractfile(
        'DirTABLES-EOS2019/TABLE_H_TP_v1').read().decode().splitlines()
        if line.strip() and not line.startswith('#')]
a, b = [next(r for r in rows if abs(r[0] - t) < 1e-12 and abs(r[2]) < 1e-12)
        for t in (5.95, 6.0)]
T0, T1 = 10**a[0], 10**b[0]
U0, U1 = 10**(a[3]+10), 10**(b[3]+10)
S1 = 10**(b[4]+10)
cp = S1 * b[7]
cv = cp + S1 * b[8] * (-b[5] / b[6])
h = next(r for r in hrows if abs(r[0] - 4.45) < 1e-12 and abs(r[1] - 1.35) < 1e-12)
# Maxwell: S_P = -(P/rho/T)*delta, and delta = -dlnrho/dlnT_P.
# Units of P and S both carry 1e10, which cancels in this ratio.
maxwell = 10**(h[1] - h[2] - h[0] - h[4]) * h[5] / h[8] - 1
print(json.dumps(dict(source_sha256=SHA256, species='He', density_g_cm3=1,
    temperature_K=[T0,T1], internal_energy_erg_g=[U0,U1],
    secant_dU_dT_erg_g_K=(U1-U0)/(T1-T0),
    listed_entropy_cv_erg_g_K=cv, listed_entropy_cp_erg_g_K=cp,
    energy_decreases_while_temperature_increases=U1<U0,
    hydrogen_source_maxwell=dict(logT_K=h[0], logP_GPa=h[1], logrho_g_cm3=h[2],
        logS_MJ_kg_K=h[4], dlnrho_dlnT_P=h[5], dlnS_dlnP_T=h[8],
        relative_defect=maxwell)), indent=2))
