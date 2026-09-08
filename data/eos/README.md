# Equation-of-state data

## FreeEOS 3.0 material Helmholtz potential

`freeeos300_hhe_x070_potential.dat` is ember's C2 biquintic representation
of a material free-energy potential from FreeEOS 3.0.0 EOS1. It is the
current fixed-composition equilibrium EOS. The numerical source mixture
is H=.7, He=.3 (metals-as-helium approximation); it is not a full-metal,
He3-capable EOS. Original direct numerical evaluations are retained in
`sources/freeeos300_hhe_x070_raw.json.gz`. Normal builds need no Fortran.

The 305×593 source grid has .0125-dex spacing in log T and
log Q=log[rho/(T/1e6)^1.5]. The potential retains 301×589 nodes after
trimming derivative stencils, with 106 masked nodes. Temperature/density
support is strict. Radiation is removed from the source before forming
the material potential and is added analytically at runtime.

See [FREEEOS.md](../../docs/FREEEOS.md) for equations, source provenance,
reproduction from the checksum-pinned external source archive, and the
small source-fit discontinuities regularized by this representation.
Sampled response discrepancies reach 1–2% near those joins; mathematical
consistency does not establish exact source reproduction or physical accuracy.

```
python3 scripts/import_freeeos_potential.py data/eos/sources/freeeos300_hhe_x070_raw.json.gz /tmp/freeeos-potential.dat
```

SHA-256:

- FreeEOS 3.0.0 source archive: `4ab1c15a51385a3eab3b08c6f3f240739c0105d92ec828d635ac95720edefb09`
- Raw evaluations: `26165748e1493438f4d62d1e04cb2eb67761a51d57c7abbf99cc0f1a01ff1e62`
- Potential: `34f5e83c898ab74273dd3d907b2ebac77c61976b8e1d7c5e1fda8530c69d823a`

All raw evaluations were reproduced byte for byte from a fresh source
build on the development host; rebuilding the potential is byte-identical.
FreeEOS source is GPL-2.0-or-later and remains outside ember's build; its
source is not relicensed or linked into ember. Cite Alan Irwin and the
FreeEOS release for these numerical calculations.

## Retained CMS19 H/He tables

`cms19_h_tp.dat` and `cms19_he_tp.dat` contain the original density and
entropy columns from `TABLE_H_TP_v1` and `TABLE_HE_TP_v1` in the authors'
[DirEOS2019.tar.gz distribution](https://perso.ens-lyon.fr/gilles.chabrier/DirEOS/DirEOS2019.tar.gz),
retrieved 2026-09-07. Cite Chabrier, Mazevet & Soubiran (2019), ApJ 872, 51,
[A new equation of state for dense hydrogen-helium mixtures](https://arxiv.org/abs/1902.01852).
The archive root is `DirTABLES-EOS2019/`. Its README identifies the pure
component files as February 2019 v1; the June 2021 archive update added a
mixture table. These external scientific data are not relicensed under
ember's MIT code licence.

The authors' June 2021 interacting-mixture README recommends the 2019
additive-volume tables for low-mass stars. We use pure H and He, never the
2021 effective-hydrogen table as pure H. The optional metals-as-helium
approximation must be explicitly selected in `Cms19Eos`.

### Reproduction

```
curl -fLO https://perso.ens-lyon.fr/gilles.chabrier/DirEOS/DirEOS2019.tar.gz
python3 scripts/import_cms19.py DirEOS2019.tar.gz data/eos
python3 scripts/audit_cms19_energy.py DirEOS2019.tar.gz
```

The importer and audit verify SHA-256
`736de2a0b26c02b897fbd906504bb6aa144f08a6d0f3cff36da62970054d149b`.
Reimport was checked byte for byte against both versioned files. Python's
standard library is sufficient; tables are read directly by C++ at runtime.

### Format and support

Each file has 121 temperatures by 441 pressures, or 53,361 pairs of original
values. The format is:

```
121 441 CMS19 H: <description>     # HE: for helium
<121 log10 T[K] values>
<441 log10 P[GPa] values>
<log10 rho[g/cm3]> <log10 S[MJ/kg/K]>  # pressure varies fastest
...
```

Axes run from log T=2 to 8 and log P[GPa]=-9 to 13, both in steps of .05.
The importer retains the original printed density/entropy strings, including
unphysical rectangular corners. It neither interpolates nor repairs them.
The runtime supports a much smaller fluid subset: log T=3.2..7.3, with
each entire 4×4 interpolation stencil inside density and phase masks.
Consequently an imported rectangle does **not** imply rectangular physical
coverage. See [CMS19 implementation and audit](../../docs/CMS19.md).

### Static calculations only

Energy columns are deliberately absent. At rho=1 g/cm³, the original He
table's internal energy decreases between log T=5.95 and 6.0 despite a
positive entropy-derived heat capacity. The source H table also violates
the pressure/entropy Maxwell identity by 26.5% at one audited dense
ionization point. These are reproducible diagnostics of this particular
source version, not a claim that every CMS19 state has such errors.

`Cms19Eos::has_internal_energy()` is false. Energy and its density
derivative are NaN; time-dependent structure and central boundary calls
reject this EOS. Density, entropy and their derivatives support an
**experimental static calculation**, with consistency errors exposed in
the equilibrium JSON. Do not manufacture internal energy or enforce one
identity by overwriting an independently derived response.
